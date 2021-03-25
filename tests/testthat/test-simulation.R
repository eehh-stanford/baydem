# Run a simulation to ensure nothing major is broken in the core code
# Also, use the outputs of the simulation to do unit testing on other
# functions.

# Create a standalone function to run the simulation so that expect_error can be
# used to ensure the simulated data can be created.
run_simulation <- function() {

  # Set the random number seed (seed from random.org)
  set.seed(466985)

  # The simulation parameters: a Gaussian mixture with ordering pi1, pi2, mu1,
  # mu2, sig1, sig2
  th_sim <- c(.2, .8, 775, 1000, 35, 45)

  # Load the calibration data frame by calling bd_load_calib_curve
  calibDf <- baydem::bd_load_calib_curve("intcal20")

  # For simulating radiocarbon measurements, a draw is made for the standard
  # deviation of the fraction modern from a uniform density on the interval 0.0021
  # to 0.0028. This is specified via the list errorSpec
  errorSpec <- list(type = "unif_fm", min = .0021, max = .0028)

  # taumin and taumax are the minimum and maximum calendar date used in the
  # Bayesian inference. Dates older than taumin and younger than taumax are
  # assumed to have zero probability, and the parameterized density p(t|theta) is
  #  normalized to integrate to 1 on the interval taumin to taumax.
  taumin <- 600 # AD
  taumax <- 1300 # AD

  # Use 400 samples, of which 200 are warm-up (sometimes called burn-in)
  # More samples would be needed to test the functionality of the inference, but
  # this test just ensures nothing basic is broken in the code.
  mcSamp <- 400 # This includes warm-up
  mcWarmup <- 200 # Number of warm-up samples

  N <- 20 # number of samples
  sampDatesAD <- bd_sample_gauss_mix(N, th_sim, taumin, taumax)
  sampRcMeas <- bd_draw_rc_meas_using_date(sampDatesAD, calibDf, errorSpec, isAD = T)
  simData <- list(datesAD = sampDatesAD, rcMeas = sampRcMeas)

  # hyperparameters for inference
  hp <-
    list(
      # Class of fit (Gaussian mixture)
      fitType = "gaussmix",
      # Parameter for the dirichlet draw of the mixture probabilities
      alpha_d = 1,
      # The gamma distribution shape parameter for sigma
      alpha_s = 3,
      # The gamma distribution rate parameter for sigma, yielding a mode of 300
      alpha_r = (3 - 1) / 300,
      # Minimum calendar date (years BC/AD)
      taumin = taumin,
      # Maximum calendar date (years BC/AD)
      taumax = taumax,
      # Spacing for the measurement matrix (years)
      dtau = 5,
      # Number of mixtures
      K = 2
    )

  # The control settings to use for inference in the call to bd_do_inference.
  # Set the total number of samples (including warmup) and warmup samples to the
  # values specified above
  controlList <- list(
    sampsPerChain = mcSamp,
    warmup = mcWarmup,
    stanControl = list(adapt_delta = .99)
  )

  # Define a problem object to be input for the inference
  prob <- list(
    phi_m = simData$rcMeas$phi_m,
    sig_m = simData$rcMeas$sig_m,
    hp = hp,
    calibDf = calibDf,
    control = controlList
  )

  soln <- bd_do_inference(prob)
  anal <- bd_analyze_soln(soln)

  return(
    list(
      prob = prob,
      soln = soln,
      anal = anal,
      calibDf = calibDf,
      errorSpec = errorSpec,
      th_sim = th_sim
    )
  )
}

# Calling run_simulation should not raise an error. If it does, the test fails.
expect_error(
  simOutput <- run_simulation(),
  NA
)

# Check that building a density plot from individual functions does not raise
# an error. This checks:
#
# bd_make_blank_density_plot
# bd_plot_50_percent_quantile
# bd_add_shaded_quantiles
expect_error(
  bd_make_blank_density_plot(simOutput$anal),
  NA
)

expect_error(
  bd_plot_50_percent_quantile(simOutput$anal, add = T),
  NA
)

expect_error(
  bd_add_shaded_quantiles(simOutput$anal),
  NA
)

# Check that calling bd_draw_rc_meas_using_date does not raise an error,
# whether isAD is True or False
t_e_AD <- c(700, 705)
t_e_calBP <- 1950 - t_e_AD
expect_error(
  rcMeas1 <- bd_draw_rc_meas_using_date(t_e_AD, simOutput$calibDf, simOutput$errorSpec, isAD = T),
  NA
)

expect_error(
  rcMeas2 <- bd_draw_rc_meas_using_date(t_e_calBP, simOutput$calibDf, simOutput$errorSpec),
  NA
)

# Check that calling bd_calc_half_life_from_peak does not raise an error
expect_error(
  halfLife <- bd_calc_half_life_from_peak(simOutput$soln),
  NA
)

expect_error(
  halfLife2 <- bd_calc_half_life_from_peak(simOutput$soln, propChange = .25),
  NA
)

# Check that calling bd_calc_relative_density does not raise an error. Do two
# tests using all three ways of specifying the relative density to also check
# the three helper functions,
#
# bd_calc_point_density
# bd_calc_range_density
# bd_calc_peak_density
expect_error(
  relDens1 <- bd_calc_relative_density(simOutput$soln, "peak", 1100),
  NA
)

expect_error(
  relDens2 <- bd_calc_relative_density(simOutput$soln, 900, c(700, 750)),
  NA
)

# Check that calling bd_extract_param does not raise an error
expect_error(
  TH <- bd_extract_param(simOutput$soln$fit),
  NA
)

# Check that calling bd_calc_gauss_mix_pdf does not raise an error for all the
# valid calculation types. Also check the dimensions of the output. Use the
# 100-th sample of TH, which has dimensions 800 x 6. tau, which is created to do
# the test, has length 141.
tau <- seq(simOutput$prob$hp$taumin, simOutput$prob$hp$taumax, by = simOutput$prob$hp$dtau)
expect_error(
  pdfVect1 <- bd_calc_gauss_mix_pdf(TH[100, ], tau), # density is the default type
  NA
)

expect_equal(
  length(tau),
  length(pdfVect1)
)

expect_error(
  pdfVect2 <- bd_calc_gauss_mix_pdf(TH[100, ], tau, type = "cumulative"),
  NA
)

expect_equal(
  length(tau),
  length(pdfVect2)
)

expect_error(
  pdfVect3 <- bd_calc_gauss_mix_pdf(TH[100, ], tau, type = "derivative"),
  NA
)

expect_equal(
  length(tau),
  length(pdfVect3)
)

expect_error(
  pdfVect4 <- bd_calc_gauss_mix_pdf(TH[100, ], tau, type = "rate"),
  NA
)

expect_equal(
  length(tau),
  length(pdfVect4)
)

# Check that calling bd_summarize_trunc_gauss_mix_sample does not raise an
# error. Use the 100-th sample from TH.
expect_error(
  summ <- bd_summarize_trunc_gauss_mix_sample(TH[100, ], simOutput$prob$hp$taumin, simOutput$prob$hp$taumax),
  NA
)

# Check that calling bd_calc_gauss_mix_pdf_mat does not raise an error for all the
# valid calculation types. Also check the dimensions of the output.
expect_error(
  pdfMat1 <- bd_calc_gauss_mix_pdf_mat(TH, tau), # density is the default type
  NA
)

expect_equal(
  c(nrow(TH), length(tau)),
  dim(pdfMat1)
)

expect_error(
  pdfMat2 <- bd_calc_gauss_mix_pdf_mat(TH, tau, type = "cumulative"),
  NA
)

expect_equal(
  c(nrow(TH), length(tau)),
  dim(pdfMat2)
)

expect_error(
  pdfMat3 <- bd_calc_gauss_mix_pdf_mat(TH, tau, type = "derivative"),
  NA
)

expect_equal(
  c(nrow(TH), length(tau)),
  dim(pdfMat3)
)

expect_error(
  pdfMat4 <- bd_calc_gauss_mix_pdf_mat(TH, tau, type = "rate"),
  NA
)

expect_equal(
  c(nrow(TH), length(tau)),
  dim(pdfMat4)
)

# Check the functioning of bd_calc_meas_matrix.

# (1) Calculate the measurement matrix using tau, already defined above, which
#     is regularly spaced. Do this for both using and not using calibration
#     uncertainty. Also check the dimensions of the output

expect_error(
  measMat1a <- bd_calc_meas_matrix(tau, simOutput$prob$phi_m, simOutput$prob$sig_m, simOutput$calibDf, addCalibUnc = T),
  NA
)

expect_equal(
  c(length(simOutput$prob$phi_m), length(tau)),
  dim(measMat1a)
)

expect_error(
  measMat1b <- bd_calc_meas_matrix(tau, simOutput$prob$phi_m, simOutput$prob$sig_m, simOutput$calibDf, addCalibUnc = F),
  NA
)

expect_equal(
  c(length(simOutput$prob$phi_m), length(tau)),
  dim(measMat1b)
)

# (2) Calculate the measurement matrix using tau_irreg, newly defined, which
#     is irregularly spaced. Do this for both using and not using calibration
#     uncertainty. Also check the dimensions of the output

tau_irreg <- c(800, 805, 810, 825)
# an error is expected if useTrapez is not set to TRUE
expect_error(
  measMat2a <- bd_calc_meas_matrix(tau_irreg, simOutput$prob$phi_m, simOutput$prob$sig_m, simOutput$calibDf),
  "tau is irregularly spaced but useTrapez is FALSE"
)

expect_error(
  measMat2a <- bd_calc_meas_matrix(tau_irreg, simOutput$prob$phi_m, simOutput$prob$sig_m, simOutput$calibDf, addCalibUnc = T, useTrapez = T),
  NA
)

expect_equal(
  c(length(simOutput$prob$phi_m), length(tau_irreg)),
  dim(measMat2a)
)

expect_error(
  measMat2b <- bd_calc_meas_matrix(tau_irreg, simOutput$prob$phi_m, simOutput$prob$sig_m, simOutput$calibDf, addCalibUnc = F, useTrapez = T),
  NA
)

expect_equal(
  c(length(simOutput$prob$phi_m), length(tau_irreg)),
  dim(measMat2b)
)

# Check the functioning of bd_calc_quantiles. Check the dimensions for a
# a call in which probs is the default, c(.025, .5, .975), and for which it is
# explicitly set.
fMat <- bd_calc_gauss_mix_pdf_mat(TH, tau)

expect_error(
  Qdens1 <- bd_calc_quantiles(fMat),
  NA
)

expect_equal(
  dim(Qdens1),
  c(3, ncol(fMat))
)

expect_error(
  Qdens2 <- bd_calc_quantiles(fMat, c(.025, .05, .5, .90, .975)),
  NA
)
expect_equal(
  dim(Qdens2),
  c(5, ncol(fMat))
)

# Check that bd_sample_gauss_mix does not throw an error. Also check the lengths
# of the outputs and that no errors are whether or not truncation is used.
expect_error(
  gmSamp1 <- bd_sample_gauss_mix(50, simOutput$th_sim),
  NA
)

expect_equal(
  length(gmSamp1),
  50
)

expect_error(
  gmSamp2 <- bd_sample_gauss_mix(50, simOutput$th_sim, simOutput$prob$hp$taumin, simOutput$prob$hp$taumax),
  NA
)

expect_equal(
  length(gmSamp2),
  50
)

# Check that bd_plot_summed_density does not throw an error.
expect_error(
  bd_plot_summed_density(simOutput$anal),
  NA
)

# Check that bd_plot_known_sim_density does not throw an error.
expect_error(
  bd_plot_known_sim_density(simOutput$anal),
  NA
)

