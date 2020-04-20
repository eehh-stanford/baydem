library(baydem)
library(magrittr)

## Generate simulated data.
# Sample sizes of the simulations
sim_samp <- c(10, 100, 1000, 10000)

# The simulation distribution: a Gaussian mixture with ordering pi1, pi2, mu1,
# mu2, sig1, sig2
th_sim <-
  c(
    pi1 = 0.2,
    pi2 = 0.8,
    mu1 = 775,
    mu2 = 1000,
    sig1 = 35,
    sig2 = 45
  )

# Set the hyperparameters
hp <-
  list(
    # Class of fit (Gaussian mixture)
    fitType = "gaussmix",
    # Parameter for the dirichlet draw of the mixture probabilities
    alpha_d = 1,
    # The gamma distribution shape parameter for sigma
    alpha_s = 3,
    # The gamma distribution rate parameter for sigma, yielding a mode of 100 or 300, respectively
    alpha_r = list(
      (3 - 1) / 100,
      (3 - 1) / 300
    ),
    # Minimum calendar date (years BC/AD)
    taumin = 600,
    # Maximum calendar date (years BC/AD)
    taumax = 1300,
    # Spacing for the measurement matrix (years)
    dtau = 1,
    # Number of mixtures
    K = 2
  ) %>%
  purrr::cross()

# Generate the largest sample;
# We'll subset this in the function below, simulating
# the process of retrieving more radiocarbon dates
# Set the random number seed (seed from random.org)
set.seed(806372)
sim_dates <-
  tibble::tibble(
    date_AD = baydem::bd_sample_gauss_mix(
      # A really large number of samples from which to draw the
      # test datasets, in case someone wants to run more than 10,000
      N = 100000,
      th = th_sim,
      taumin = hp[[1]]$taumin,
      taumax = hp[[1]]$taumax
    )
  ) %>%
  dplyr::bind_cols(
    .,
    baydem::bd_draw_rc_meas_using_date(
      t_e = .$date_AD,

      # Load the calibration data frame by calling bd_load_calib_curve
      calibDf = bd_load_calib_curve("intcal13"),

      # For simulating radiocarbon measurements, a draw is made for the standard
      # deviation of the fraction modern from a uniform density on the interval 0.0021
      # to 0.0028. This is specified via the list errorSpec
      errorSpec = list(
        min = .0021,
        max = .0028
      ),
      isAD = T
    ) %>%
      tibble::as_tibble()
  )

sim_inference <-
  sim_samp %>%
  magrittr::set_names(., .) %>%
  as.list() %>%
  purrr::cross2(hp) %>%
  purrr::map(function(x) {
    x %<>%
      magrittr::set_names(c("n", "hp"))

    x$sim_dates <- sim_dates[1:x$n, ]

    x
  })

# If simulated data exists for runs, load it. If are missing for any run, generate
# new data for all runs.
# Doing inference involves three steps:
#
# (1) Generate the problem
# (2) Do the Bayesian sampling
# (3) Run some standard analyses
out_file <- here::here("analysis/data-derived/sim_inference.rds")

if (
  !file.exists(out_file) ||
    !identical(
      out_file %>%
        readr::read_rds() %>%
        purrr::map(magrittr::extract, c("n", "hp", "sim_dates")),
      sim_inference
    )
) {
  if (file.exists(out_file)) {
    saved_results <-
      out_file %>%
      readr::read_rds() %>%
      purrr::map(magrittr::extract, c("n", "hp", "sim_dates"))

    sim_inference %<>%
      setdiff(saved_results)
  }

  sim_inference %<>%
    purrr::map(
      function(x) {
        prob <-
          list(
            phi_m = x$sim_dates$phi_m,
            sig_m = x$sim_dates$sig_m,
            hp = x$hp,
            calibDf = bd_load_calib_curve("intcal13"),
            # Define the control parameters for the call to Stan. Use 4500 total MCMC
            # samples, of which 2000 are warmup samples. Since four chains are used, this
            # yields 4*(4500-2000) = 10,000 total samples.
            control = list(
              sampsPerChain = 4500,
              warmup = 2000
            )
          )

        soln <-
          baydem::bd_do_inference(prob)

        anal <-
          baydem::bd_analyze_soln(
            soln = soln,
            th_sim = th_sim
          )

        x$sim_output <-
          tibble::lst(
            prob,
            soln,
            anal
          )

        return(x)
      }
    )

  if (file.exists(out_file)) {
    saved_results <-
      out_file %>%
      readr::read_rds()

    # Discard whatever results were just calculated
    saved_results %<>%
      purrr::discard(saved_results %>%
        purrr::map(magrittr::extract, c("n", "hp", "sim_dates")) %>%
        magrittr::is_in(sim_inference %>%
          purrr::map(magrittr::extract, c("n", "hp", "sim_dates"))))

    sim_inference <- base::union(
      saved_results,
      sim_inference
    )
  }

  # Save the full result set
  sim_inference %>%
    readr::write_rds(out_file,
      compress = "gz"
    )
}



# Load the data.
sim_inference <-
  "analysis/data-derived/sim_inference.rds" %>%
  here::here() %>%
  readr::read_rds()

print("Key parameter values:")
sim_inference %>%
  purrr::map_dfr(function(x) {
    x[c("n", "hp")] %>%
      unlist(recursive = FALSE) %>%
      tibble::as_tibble()
  }) %>%
  dplyr::arrange(n) %T>%
  readr::write_csv(here::here("analysis/logs/sim_inference.csv"))




#### Make Simulated Plots ####
sim_inference %<>%
  purrr::keep(function(x) x$hp$alpha_r == ((3 - 1) / 300)) %>%
  purrr::keep(function(x) x$n != 10)

nplots <- length(sim_inference) + 1
# Generate a 4 x 1 graph figure summarizing the simulation results
pdf(here::here("analysis/figures/Fig1_sim_inference.pdf"), width = 5, height = 2.5 * nplots)

par(
  mfrow = c(nplots, 1),
  xaxs = "i", # No padding for x-axis
  yaxs = "i", # No padding for y-axis
  # outer margins with ordering bottom, left, top, right:
  oma = c(4, 2, 2, 2),
  # plot margins with ordering bottom, left, top, right:
  mar = c(2, 4, 0, 0)
  # Don't add data if it falls outside plot window
  # xpd = F
)

# (1) Calibration curve
par(mar = c(0, 4, 0, 0))
bd_vis_calib_curve(min(sim_inference[[1]]$sim_output$anal$tau),
  max(sim_inference[[1]]$sim_output$anal$tau),
  sim_inference[[1]]$sim_output$prob$calibDf,
  xlab = "",
  ylab = "Fraction Modern",
  xaxt = "n",
  invertCol = "gray80"
)
box()

# (2-nplots) Density plots
sim_inference %>%
  purrr::walk(function(x) {
    par(mar = c(0, 4, 0, 0))

    bd_make_blank_density_plot(x$sim_output$anal,
      ylim = c(0, 0.01),
      xlab = "",
      ylab = "Density",
      xaxt = "n",
      yaxt = "n"
    )

    bd_add_shaded_quantiles(x$sim_output$anal,
      col = "gray80"
    )

    bd_plot_summed_density(x$sim_output$anal,
      lwd = 2,
      add = T,
      col = "black"
    )

    bd_plot_50_percent_quantile(x$sim_output$anal,
      lwd = 2,
      add = T,
      col = "red"
    )

    bd_plot_known_sim_density(x$sim_output$anal,
      lwd = 2,
      add = T,
      col = "blue"
    )

    text(
      labels = paste0("n = ", x$n),
      x = 600,
      y = 0.009,
      pos = 4,
      cex = 2
    )

    axis(
      side = 2,
      at = c(0, 0.002, 0.004, 0.006, 0.008)
    )

    box()
  })

axis(side = 1)
mtext("Calendar Date [AD]", side = 1, line = 2.5, cex = 0.75)

dev.off()
