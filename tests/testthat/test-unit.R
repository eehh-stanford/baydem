# Unit tests for baydem functions. These tests are fast but do not check
# function. Thus, they can be quickly run as code is actively developed. The
# function tests are located in test-function.R. They take much longer to run
# and should primarily be run before pushing updates back to github.

# ------------------------------------------------------------------------------
# (1) Do unit tests for functions that generate simulated data
#
# Coverage: sample_gauss_mix
#           draw_rc_meas_using_date
#           simulate_rc_data
# ------------------------------------------------------------------------------

# (1a) test sample_gauss_mix
N <- 75
th_sim <-
  c(
    pi1 = 0.2,
    pi2 = 0.8,
    mu1 = 775,
    mu2 = 1000,
    sig1 = 35,
    sig2 = 45
  )

expect_error(
  sample <- sample_gauss_mix(N,
                             th_sim),
  NA
)

expect_equal(
  length(sample),
  N
)

tau_min <- 600
tau_max <- 1300
expect_error(
  sample <- sample_gauss_mix(N,
                             th_sim,
                             tau_min,
                             tau_max),
  NA
)

expect_equal(
  length(sample),
  N
)

expect_error(
  sample_gauss_mix(N, 1:15),
  "The maximum number of supported mixture components is 4"
)

# (1b) test draw_rc_meas_using_date
calib_df <- load_calib_curve("intcal20")
# test with is_AD FALSE
expect_error(
  rc_meas <- draw_rc_meas_using_date(sample,
                                     calib_df,
                                     list(type="unif_fm",
                                             min=.0021,
                                             max=.0028),
                                     is_AD=F),
  NA
)

expect_equal(
  names(rc_meas),
  c("phi_m","sig_m","trc_m","sig_trc_m")
)

for(vector_name in names(rc_meas)) {
  expect_equal(
    length(rc_meas[[vector_name]]),
    N
  )
}

expect_error(
  rc_meas <- draw_rc_meas_using_date(sample,
                                     calib_df,
                                     list(type="not_valid",
                                             min=.0021,
                                             max=.0028),
                                     is_AD=F),
  "Unrecognized error type = not_valid"
)

# test with is_AD TRUE
expect_error(
  rc_meas <- draw_rc_meas_using_date(sample,
                                     calib_df,
                                     list(type="unif_fm",
                                             min=.0021,
                                             max=.0028),
                                     is_AD=T),
  NA
)

expect_equal(
  names(rc_meas),
  c("phi_m","sig_m","trc_m","sig_trc_m")
)

for(vector_name in names(rc_meas)) {
  expect_equal(
    length(rc_meas[[vector_name]]),
    N
  )
}

expect_error(
  rc_meas <- draw_rc_meas_using_date(sample,
                                     calib_df,
                                     list(type="not_valid",
                                             min=.0021,
                                             max=.0028),
                                     is_AD=T),
  "Unrecognized error type = not_valid"
)

# (1c) test simulate_rc_data

# Check simulation of a non-truncated Gaussian mixture with no random number
# seed set
sim_spec <- list(model_spec=
                   list(density_type = "gauss_mix",
                        th=th_sim,
                        error_spec=list(type="unif_fm",min=.0021,max=.0028),
                        is_AD=T),
                 N=N,
                 calib_curve="intcal20")

expect_error(
  sim <- simulate_rc_data(sim_spec),
  NA
)

expect_equal(
  names(sim),
  c("sim_spec","data")
)

expect_equal(
  names(sim$data),
  c("dates","rc_meas")
)

expect_equal(
  length(sim$data$dates),
  N
)

# Check simulation of a truncated Gaussian mixture with no random number seed
# set
sim_spec <- list(model_spec=
                   list(density_type = "trunc_gauss_mix",
                        th=c(th_sim,tau_min,tau_max),
                        error_spec=list(type="unif_fm",min=.0021,max=.0028),
                        is_AD=T),
                 N=N,
                 calib_curve="intcal20")

expect_error(
  sim <- simulate_rc_data(sim_spec),
  NA
)

expect_equal(
  names(sim),
  c("sim_spec","data")
)

expect_equal(
  names(sim$data),
  c("dates","rc_meas")
)

expect_equal(
  length(sim$data$dates),
  N
)

# Check simulation of a truncated Gaussian mixture with a random number seed
# set
sim_spec <- list(model_spec=
                   list(density_type = "trunc_gauss_mix",
                        th=c(th_sim,tau_min,tau_max),
                        error_spec=list(type="unif_fm",min=.0021,max=.0028),
                        is_AD=T),
                 N=N,
                 calib_curve="intcal20",
                 seed=1002)

expect_error(
  sim <- simulate_rc_data(sim_spec),
  NA
)

expect_equal(
  names(sim),
  c("sim_spec","data")
)

expect_equal(
  names(sim$data),
  c("dates","rc_meas")
)

expect_equal(
  length(sim$data$dates),
  N
)

# ------------------------------------------------------------------------------
# (2) Do unit tests for functions related to the maximum likelihood fitting
#
# coverage: is_th_reduced_valid
#           calc_trunc_gauss_mix_neg_log_lik
#           temper_trunc_gauss_mix
# ------------------------------------------------------------------------------

# (2a) is_th_reduced_valid
expect_equal(
  is_th_reduced_valid(c(0.5,10,20,2,3)),
  TRUE
)

# Fails for pi
expect_equal(
  is_th_reduced_valid(c(1.5,10,20,2,3)),
  FALSE
)

# Fails for mu
expect_equal(
  is_th_reduced_valid(c(0.5,20,10,2,3)),
  FALSE
)

# Fails for s
expect_equal(
  is_th_reduced_valid(c(0.5,10,20,2,3),s_min=10),
  FALSE
)

# (2b) calc_trunc_gauss_mix_neg_log_lik
tau <- seq(tau_min,tau_max,by=5)
M <- calc_meas_matrix(tau,
                      sim$data$rc_meas$phi_m,
                      sim$data$rc_meas$sig_m,
                      calib_df)
expect_error(
  neg_log_lik <- calc_trunc_gauss_mix_neg_log_lik(th_sim[-1],M,tau),
  NA
)

expect_equal(
  is.na(neg_log_lik),
  FALSE
)
expect_equal(
  is.finite(neg_log_lik),
  TRUE
)

expect_equal(
  calc_trunc_gauss_mix_neg_log_lik(th_sim[-1],M,tau,s_min=100),
  Inf
)

# (2c) fit_trunc_gauss_mix
# Check that fit_trunc_gauss_mix does not throw an error (and also that the
# output has the right fields). Do this for both not using and using multiple
# cores. Also check the behavior and reproducibility when setting the random
# number seed(s).

expect_error(
  max_lik_fit <-
    fit_trunc_gauss_mix(
                           2,
                           sim$data$rc_meas$phi_m,
                           sim$data$rc_meas$sig_m,
                           600,
                           1300,
                           1,
                           calib_df,
                           num_restarts=3),
  NA
)

expect_equal(
  names(max_lik_fit),
  c("th","neg_log_lik","tau","f","bic","aic","hjkb_fit_list",
    "input_seed","base_seed","seed_vect")
)

expect_error(
  max_lik_fit_a <-
    fit_trunc_gauss_mix(
                           2,
                           sim$data$rc_meas$phi_m,
                           sim$data$rc_meas$sig_m,
                           600,
                           1300,
                           1,
                           calib_df,
                           num_restarts=3,
                           input_seed=10),
  NA
)

expect_error(
  max_lik_fit_b <-
    fit_trunc_gauss_mix(
                           2,
                           sim$data$rc_meas$phi_m,
                           sim$data$rc_meas$sig_m,
                           600,
                           1300,
                           1,
                           calib_df,
                           num_restarts=3,
                           input_seed=10),
  NA
)

expect_equal(
  max_lik_fit_a,
  max_lik_fit_b
)

expect_error(
  max_lik_fit_a <-
    fit_trunc_gauss_mix(
                           2,
                           sim$data$rc_meas$phi_m,
                           sim$data$rc_meas$sig_m,
                           600,
                           1300,
                           1,
                           calib_df,
                           num_restarts=3,
                           input_seed=10:13),
  NA
)

expect_error(
  max_lik_fit_b <-
    fit_trunc_gauss_mix(
                           2,
                           sim$data$rc_meas$phi_m,
                           sim$data$rc_meas$sig_m,
                           600,
                           1300,
                           1,
                           calib_df,
                           num_restarts=3,
                           input_seed=10:13),
  NA
)

expect_equal(
  max_lik_fit_a,
  max_lik_fit_b
)

expect_error(
  max_lik_fit <-
    fit_trunc_gauss_mix(
                        2,
                        sim$data$rc_meas$phi_m,
                        sim$data$rc_meas$sig_m,
                        600,
                        1300,
                        1,
                        calib_df,
                        num_restarts=3,
                        input_seed=10:12),
  paste0("input_seed must be NA, a single integer, or a vector of ",
         "integers with length 1 + num_restarts"),
  fixed = TRUE
)

# Multiple cores
expect_error(
  max_lik_fit <-
    fit_trunc_gauss_mix(2,
                        sim$data$rc_meas$phi_m,
                        sim$data$rc_meas$sig_m,
                        600,
                        1300,
                        1,
                        calib_df,
                        num_restarts=3,
                        num_cores=2),
  NA
)

expect_equal(
  names(max_lik_fit),
  c("th","neg_log_lik","tau","f","bic","aic","hjkb_fit_list",
    "input_seed","base_seed","seed_vect")
)

expect_error(
  max_lik_fit_a <-
    fit_trunc_gauss_mix(
                           2,
                           sim$data$rc_meas$phi_m,
                           sim$data$rc_meas$sig_m,
                           600,
                           1300,
                           1,
                           calib_df,
                           num_restarts=3,
                           num_cores=2,
                           input_seed=10),
  NA
)

expect_error(
  max_lik_fit_b <-
    fit_trunc_gauss_mix(
                           2,
                           sim$data$rc_meas$phi_m,
                           sim$data$rc_meas$sig_m,
                           600,
                           1300,
                           1,
                           calib_df,
                           num_restarts=3,
                           num_cores=2,
                           input_seed=10),
  NA
)

expect_equal(
  max_lik_fit_a,
  max_lik_fit_b
)

expect_error(
  max_lik_fit_a <-
    fit_trunc_gauss_mix(
                           2,
                           sim$data$rc_meas$phi_m,
                           sim$data$rc_meas$sig_m,
                           600,
                           1300,
                           1,
                           calib_df,
                           num_restarts=3,
                           num_cores=2,
                           input_seed=10:13),
  NA
)

expect_error(
  max_lik_fit_b <-
    fit_trunc_gauss_mix(
                           2,
                           sim$data$rc_meas$phi_m,
                           sim$data$rc_meas$sig_m,
                           600,
                           1300,
                           1,
                           calib_df,
                           num_restarts=3,
                           num_cores=2,
                           input_seed=10:13),
  NA
)

expect_equal(
  max_lik_fit_a,
  max_lik_fit_b
)

expect_error(
  max_lik_fit <-
    fit_trunc_gauss_mix(
                        2,
                        sim$data$rc_meas$phi_m,
                        sim$data$rc_meas$sig_m,
                        600,
                        1300,
                        1,
                        calib_df,
                        num_restarts=3,
                        num_cores=2,
                        input_seed=10:12),
  paste0("input_seed must be NA, a single integer, or a vector of ",
         "integers with length 1 + num_restarts"),
  fixed = TRUE
)

# ------------------------------------------------------------------------------
# (3) Do unit tests for functions related to the Bayesian inference
#
# coverage: do_bayesian_inference
#           summarize_bayesian_inference
# ------------------------------------------------------------------------------


# (3a) do_bayesian_inference
# Define the density model
density_model <- list(type="trunc_gauss_mix",
                      tau_min=600,
                      tau_max=1300,
                      K=2)

# Define the hyperparameters
hp <-
  list(
    # Parameter for the dirichlet draw of the mixture probabilities
    alpha_d = 1,
    # The gamma distribution shape parameter for sigma
    alpha_s = 3,
    # The gamma distribution rate parameter for sigma, yielding a mode of 300
    alpha_r = (3 - 1) / 300,
    # Spacing for the measurement matrix (years)
    dtau = 5
  )

input_control <- list(samps_per_chain=2000,warmup=1000,
                      stan_control=list(adapt_delta=.99))
expect_error(
  soln <- do_bayesian_inference(
    sim$data$rc_meas,
    density_model,
    hp,
    calib_df,
    control=input_control
  ),
  NA
)

expect_is(
  soln,
  "bd_bayesian_soln"
)

expect_equal(
  names(soln),
  c("fit","control","init_seed","stan_seed","optional_inputs")
)

expect_is(
  soln$fit,
  "stanfit"
)

expect_equal(
  names(soln$control),
  c("num_chains","samps_per_chain","warmup","stan_control")
)

expect_equal(
  soln$optional_inputs,
  list(th0=NA,init_seed=NA,stan_seed=NA,control=input_control)
)

# (3b) summarize_bayesian_inference
expect_error(
  summary <- summarize_bayesian_inference(soln,
                                          rc_meas,
                                          density_model,
                                          calib_df),
  NA
)

expect_is(
  summary,
  "bd_bayesian_summary"
)

expect_equal(
  names(summary),
  c("tau","f_spdf","Qdens","Qrate","probs",
    "rate_prop","rate_ind","growth_state","summ_list")
)

expect_equal(
  dim(summary$Qdens),
  c(3,length(seq(density_model$tau_min,density_model$tau_max,by=5)))
)

# ------------------------------------------------------------------------------
# (4) Do unit tests for plotting functions
#
# coverage: make_blank_density_plot
#           plot_50_percent_quantile
#           add_shaded_quantiles
#           plot_summed_density
#           plot_known_sim_density
# ------------------------------------------------------------------------------

# (4a) make_blank_density_plot
expect_error(
  make_blank_density_plot(summary),
  NA
)

# (4b) plot_50_percent_quantile
expect_error(
  plot_50_percent_quantile(summary, add = T),
  NA
)

# (4c) add_shaded_quantiles
expect_error(
  add_shaded_quantiles(summary),
  NA
)

# (4d) plot_known_sim_density
#      plot_summed_density
expect_error(
  plot_known_sim_density(summary, add=T),
  NA
)

expect_error(
  plot_summed_density(summary, add=T),
  NA
)

expect_error(
  plot_known_sim_density(summary),
  NA
)

expect_error(
  plot_summed_density(summary),
  NA
)

# (4e) vis_calib_curve
expect_error(
  vis_calib_curve(600, 700, calib_df),
  NA
)

# ------------------------------------------------------------------------------
# (5) Do unit tests for density calculations and related "helper functions"
#
# coverage: extract_param
#           calc_gauss_mix_pdf
#           calc_gauss_mix_pdf_mat
#           calc_quantiles
#           calc_half_life_from_peak
#           calc_relative_density
#             unpack_spec           [indirect]
#             calc_point_density [indirect]
#             calc_range_density [indirect]
#             calc_peak_density  [indirect]
#           summarize_trunc_gauss_mix_sample
# ------------------------------------------------------------------------------

# (5a) extract_param
total_samples <- 4000
expect_error(
  TH <- extract_param(soln$fit),
  NA
)

expect_equal(
  dim(TH),
  c(total_samples,6)
)

# (5b) calc_gauss_mix_pdf
tau <- seq(density_model$tau_min,density_model$tau_max,by=hp$dtau)

# density is the default type
expect_error(
  pdf_vect <- calc_gauss_mix_pdf(TH[50, ], tau),
  NA
)

expect_equal(
  length(tau),
  length(pdf_vect)
)

expect_equal(
  any(is.na(pdf_vect)),
  FALSE
)

expect_error(
  cdf_vect <- calc_gauss_mix_pdf(TH[50, ], tau, type = "cumulative"),
  NA
)

expect_equal(
  length(tau),
  length(cdf_vect)
)

expect_equal(
  any(is.na(cdf_vect)),
  FALSE
)

expect_error(
  deriv_vect <- calc_gauss_mix_pdf(TH[50, ], tau, type = "derivative"),
  NA
)

expect_equal(
  length(tau),
  length(deriv_vect)
)

expect_equal(
  any(is.na(deriv_vect)),
  FALSE
)

expect_error(
  rate_vect <- calc_gauss_mix_pdf(TH[50, ], tau, type = "rate"),
  NA
)

expect_equal(
  length(tau),
  length(rate_vect)
)

expect_equal(
  any(is.na(rate_vect)),
  FALSE
)

# (5c) calc_gauss_mix_pdf_mat

# density is the default type
expect_error(
  pdf_mat <- calc_gauss_mix_pdf_mat(TH, tau),
  NA
)

expect_equal(
  c(nrow(TH), length(tau)),
  dim(pdf_mat)
)

expect_equal(
  any(is.na(pdf_mat)),
  FALSE
)

expect_error(
  cdf_mat <- calc_gauss_mix_pdf_mat(TH, tau, type = "cumulative"),
  NA
)

expect_equal(
  c(nrow(TH), length(tau)),
  dim(cdf_mat)
)

expect_equal(
  any(is.na(cdf_mat)),
  FALSE
)

expect_error(
  deriv_mat <- calc_gauss_mix_pdf_mat(TH, tau, type = "derivative"),
  NA
)

expect_equal(
  c(nrow(TH), length(tau)),
  dim(deriv_mat)
)

expect_equal(
  any(is.na(deriv_mat)),
  FALSE
)

expect_error(
  rate_mat <- calc_gauss_mix_pdf_mat(TH, tau, type = "rate"),
  NA
)

expect_equal(
  c(nrow(TH), length(tau)),
  dim(rate_mat)
)

expect_equal(
  any(is.na(rate_mat)),
  FALSE
)

# (5d) calc_quantiles

expect_error(
  Qdens1 <- calc_quantiles(pdf_mat),
  NA
)

expect_equal(
  dim(Qdens1),
  c(3, ncol(pdf_mat))
)

expect_error(
  Qdens2 <- calc_quantiles(pdf_mat, c(.025, .05, .5, .90, .975)),
  NA
)

expect_equal(
  dim(Qdens2),
  c(5, ncol(pdf_mat))
)

# (5e) calc_half_life_from_peak
expect_error(
  half_life1 <- calc_half_life_from_peak(soln,
                                         density_model,
                                         rc_meas=rc_meas,
                                         calib_df=calib_df),
  NA
)

expect_equal(
  length(half_life1),
  total_samples
)

expect_equal(
  any(is.na(half_life1)),
  FALSE
)

expect_error(
  half_life2 <- calc_half_life_from_peak(soln,
                                         density_model,
                                         rc_meas=rc_meas,
                                         calib_df=calib_df,
                                         prop_chang=.25),
  NA
)

expect_equal(
  length(half_life2),
  total_samples
)

expect_equal(
  any(is.na(half_life2)),
  FALSE
)

expect_equal(
  all(half_life1 < half_life2),
  TRUE
)

# (5f) calc_relative_density
#        unpack_spec           [indirect]
#        calc_point_density [indirect]
#        calc_range_density [indirect]
#        calc_peak_density  [indirect]
# Check two calls to calc_relative_density to check all four helper
# functions (which are only checked indirectly).

expect_error(
  rel_dens1 <- calc_relative_density(soln,
                                     density_model,
                                     "peak",
                                     1100,
                                     rc_meas=rc_meas,
                                     calib_df=calib_df),
  NA
)

expect_equal(
  length(rel_dens1),
  total_samples
)

expect_equal(
  any(is.na(rel_dens1)),
  FALSE
)

expect_error(
  rel_dens2 <- calc_relative_density(soln,
                                     density_model,
                                     900,
                                     c(700, 750),
                                     rc_meas=rc_meas,
                                     calib_df=calib_df),
  NA
)

expect_equal(
  length(rel_dens2),
  total_samples
)

expect_equal(
  any(is.na(rel_dens2)),
  FALSE
)

# (5g) summarize_trunc_gauss_mix_sample
expect_error(
  summ <- summarize_trunc_gauss_mix_sample(TH[50, ],
                                           density_model$tau_min,
                                           density_model$tau_max),
  NA
)

expect_equal(
  names(summ),
  c("periods","ind_peak","t_peak","f_peak","pattern")
)

# ------------------------------------------------------------------------------
# (6) Do unit tests for the measurement matrix calculation
#
# coverage: calc_meas_matrix
#           calc_trapez_weights
# ------------------------------------------------------------------------------


# (1) Calculate the measurement matrix using tau, already defined above, which
#     is regularly spaced. Do this for both using and not using calibration
#     uncertainty. Also check the dimensions of the output

# (6a) calc_meas_matrix
#        tau regularly spaced
#        with calibration uncertainty
expect_error(
  M6a <- calc_meas_matrix(tau,
                          rc_meas$phi_m,
                          rc_meas$sig_m,
                          calib_df,
                          add_calib_unc=T),
  NA
)

expect_equal(
  dim(M6a),
  c(length(rc_meas$phi_m), length(tau)),
)

expect_equal(
  any(is.na((M6a))),
  FALSE
)

# (6b) calc_meas_matrix
#        tau regularly spaced
#        without calibration uncertainty
expect_error(
  M6b <- calc_meas_matrix(tau,
                          rc_meas$phi_m,
                          rc_meas$sig_m,
                          calib_df,
                          add_calib_unc=F),
  NA
)

expect_equal(
  dim(M6b),
  c(length(rc_meas$phi_m), length(tau)),
)

expect_equal(
  any(is.na((M6b))),
  FALSE
)

# (6c) calc_meas_matrix
#        tau irrregularly spaced
#        with calibration uncertainty
tau_irreg <- c(800, 805, 810, 825)

# an error is expected if useTrapez is not set to TRUE
expect_error(
  calc_meas_matrix(tau_irreg,
                   rc_meas$phi_m,
                   rc_meas$sig_m,
                   calib_df,
                   add_calib_unc=T),
  "tau is irregularly spaced but useTrapez is FALSE"
)

expect_error(
  M6c <- calc_meas_matrix(
    tau_irreg,
    rc_meas$phi_m,
    rc_meas$sig_m,
    calib_df,
    add_calib_unc=T,
    use_trapez=T),
  NA
)

expect_equal(
  dim(M6c),
  c(length(rc_meas$phi_m),length(tau_irreg)),
)

expect_equal(
  any(is.na((M6c))),
  FALSE
)

# (6d) calc_meas_matrix
#        tau irrregularly spaced
#        without calibration uncertainty

# an error is expected if useTrapez is not set to TRUE
expect_error(
  calc_meas_matrix(tau_irreg,
                   rc_meas$phi_m,
                   rc_meas$sig_m,
                   calib_df,
                   add_calib_unc=F),
  "tau is irregularly spaced but useTrapez is FALSE"
)

expect_error(
  M6d <- calc_meas_matrix(
    tau_irreg,
    rc_meas$phi_m,
    rc_meas$sig_m,
    calib_df,
    add_calib_unc=F,
    use_trapez=T),
  NA
)

expect_equal(
  dim(M6d),
  c(length(rc_meas$phi_m),length(tau_irreg)),
)

expect_equal(
  any(is.na((M6d))),
  FALSE
)

# (6e) calc_trapez_weights
expect_equal(
  calc_trapez_weights(c(-1.5, 2, 3, 4, 7)),
  c(1.75, 2.25, 1, 2, 1.5)
)

# ------------------------------------------------------------------------------
# (7) Do unit tests for identifiability functions
#
# coverage: assess_calib_curve_equif
#             phi2tau [indirect]
#           calc_calib_curve_equif_dates
#           calc_calib_curve_frac_modern
# ------------------------------------------------------------------------------

# (7a) assess_calib_curve_equif
#        phi2tau [indirect]
#      calc_calib_curve_equif_dates
expect_error(
  equif_result1 <- assess_calib_curve_equif(calib_df),
  NA
)

expect_equal(
  names(equif_result1),
  c("inv_span_list","can_invert")
)

expect_equal(
  length(equif_result1$can_invert),
  nrow(calib_df)
)

expect_error(
  equif_dates <- calc_calib_curve_equif_dates(calib_df),
  NA
)

expect_error(
  equif_result2 <- assess_calib_curve_equif(calib_df, equif_dates),
  NA
)

expect_equal(
  names(equif_result2),
  c("inv_span_list","can_invert")
)

expect_equal(
  length(equif_result2$can_invert),
  nrow(calib_df)
)

expect_equal(
  equif_result1,
  equif_result2
)

# (7b) calc_calib_curve_frac_modern

expect_error(
  phi <- calc_calib_curve_frac_modern(calib_df),
  NA
)

expect_equal(
  phi,
  exp(-calib_df$uncal_year_BP / 8033)
)

tau2 <- c(600, 602, 805.89)
expect_error(
  phi2 <- calc_calib_curve_frac_modern(calib_df, tau2),
  NA
)

expect_equal(
  length(phi2),
  length(tau2)
)

# ------------------------------------------------------------------------------
# (8) Do unit tests for truncated exponential model
#
# coverage: assess_calib_curve_equif
# ------------------------------------------------------------------------------
expect_error(
  exp_samp1 <- sample_trunc_exp(50, 0.01, 600, 1300),
  NA
)

expect_equal(
  length(exp_samp1),
  50
)

expect_error(
  exp_samp2 <- sample_trunc_exp(50, -0.01, 600, 1300),
  NA
)

expect_equal(
  length(exp_samp2),
  50
)

# ------------------------------------------------------------------------------
# (9) Do unit tests for data io (input/output) functions
#
# coverage: import_rc_data    (9a)
#           set_rc_meas       (9b)
#           set_sim           (9c)
#           calc_tau_range    (9d)
#           set_density_model (9e)
#           set_calib_curve   (9f)
#           do_max_lik_fits   (9g)
# ------------------------------------------------------------------------------

# (9a) import_rc_data

# Specifying phi_m and sig_m with named columns
file_name <- system.file("extdata",
                         "fraction_modern_named_columns.csv",
                         package = "baydem")

expect_error(
  rc_meas <- import_rc_data(file_name,
                            phi_m_col="fraction_modern",
                            sig_m_col="fraction_modern_error"),
  NA
)

phi_m <- c(.80,.75,.74,.82)
sig_m <- c(.001,.002,.002,.001)
rc_meas_direct <- list(phi_m=phi_m,
                       sig_m=sig_m,
                       trc_m=-8033 * log(phi_m),
                       sig_trc_m=8033 * sig_m / phi_m)

expect_equal(
  rc_meas,
  rc_meas_direct
)

# Specifying an invalid set of inputs
expect_error(
  rc_meas <- import_rc_data(file_name,
                            trc_m_col="this_should_fail",
                            phi_m_col="fraction_modern",
                            sig_m_col="fraction_modern_error"),
  "Unsupported input pattern for specifying data columns"
)

# Specifying phi_m and sig_m with column numbers
expect_error(
  rc_meas <- import_rc_data(file_name,
                            phi_m_col=1,
                            sig_m_col=2),
  NA
)

expect_equal(
  rc_meas,
  rc_meas_direct
)

# Specifying trc_m and sig_trc_m with named columns
file_name <- system.file("extdata",
                         "rc_years_BP_named_columns.csv",
                         package = "baydem")

expect_error(
  rc_meas <- import_rc_data(file_name,
                            trc_m_col="rc_years_BP",
                            sig_trc_m_col="rc_years_BP_error"),
  NA
)

trc_m     <- c(900,800,850)
sig_trc_m <- c( 20,400, 70)
phi_m <- exp(-trc_m/8033)
sig_m <- phi_m * sig_trc_m / 8033
rc_meas_direct <- list(phi_m=phi_m,
                       sig_m=sig_m,
                       trc_m=trc_m,
                       sig_trc_m=sig_trc_m)

expect_equal(
  rc_meas,
  rc_meas_direct
)

# Specifying trc_m and sig_trc_m with column numbers and skipping the column
# names but specying header=F (to test the functioning of ...)
expect_error(
  rc_meas <- import_rc_data(file_name,
                            trc_m_col=1,
                            sig_trc_m_col=2,
                            skip=1,
                            header=F),
  NA
)

expect_equal(
  rc_meas,
  rc_meas_direct
)

# Specifying no columns, for which it is assumed that the first and second
# columns are, respectively, trc_m and sig_trc_m.

expect_error(
  rc_meas <- import_rc_data(file_name),
  NA
)

expect_equal(
  rc_meas,
  rc_meas_direct
)

# (9b) set_rc_meas

# Call data_dir to get the temporary directory to use in testing
data_dir <- tempdir()

rc_meas <- sim$data$rc_meas

analysis_name <- "test_analysis"
path_to_analysis_file <- file.path(data_dir,paste0(analysis_name,".rds"))

# If the analysis file exists, delete it
if (file.exists(path_to_analysis_file)) {
  file.remove(path_to_analysis_file)
}

expect_error(
  set_rc_meas(data_dir,analysis_name,rc_meas),
  NA
)

# Make sure the file was created and the stored variable is correct
expect_equal(
  file.exists(path_to_analysis_file),
  T
)

expect_equal(
  readRDS(path_to_analysis_file),
  list(rc_meas=rc_meas)
)

# Make sure that we cannot set rc_meas for an analysis that already exists
expect_error(
  set_rc_meas(data_dir,analysis_name,rc_meas),
  "A save file for analysis_name already exists in data_dir"
)

# (9c) set_sim

# Call data_dir to get the temporary directory to use in testing

sim_analysis_name <- "sim"
path_to_sim_analysis_file <-
  file.path(data_dir,paste0(sim_analysis_name,".rds"))

# If the analysis file exists, delete it
if (file.exists(path_to_sim_analysis_file)) {
  file.remove(path_to_sim_analysis_file)
}

expect_error(
  set_sim(data_dir,sim_analysis_name,sim),
  NA
)

# Make sure the file was created and the stored variable is correct
expect_equal(
  file.exists(path_to_sim_analysis_file),
  T
)

expect_equal(
  readRDS(path_to_sim_analysis_file),
  list(rc_meas=sim$data$rc_meas,calib_df=calib_df,sim=sim)
)

# Make sure that we cannot set a simulation for an analysis that already exists
expect_error(
  set_sim(data_dir,sim_analysis_name,sim),
  "A save file for analysis_name already exists in data_dir"
)

# (9d) calc_tau_range

# The tau range should be the same across test runs since all relative number
# seeds are set above
expect_error(
  tau_range <- calc_tau_range(rc_meas),
  NA
)

expect_equal(
  tau_range,
  list(tau_min=603,tau_max=1265)
)

# Using dtau=1 should have no effect
expect_error(
  tau_range <- calc_tau_range(rc_meas,dtau=1),
  NA
)

expect_equal(
  tau_range,
  list(tau_min=603,tau_max=1265)
)

# A warning should be thrown if dtau is not an integer
expect_warning(
  tau_range <- calc_tau_range(rc_meas,dtau=1.5),
  "dtau is being ignored because it is not an integer"
)

expect_equal(
  tau_range,
  list(tau_min=603,tau_max=1265)
)

# Check that the range is extended if dtau=5 and dtau=7
expect_error(
  tau_range <- calc_tau_range(rc_meas,dtau=5),
  NA
)

expect_equal(
  tau_range,
  list(tau_min=600,tau_max=1265)
)

expect_error(
  tau_range <- calc_tau_range(rc_meas,dtau=7),
  NA
)

expect_equal(
  tau_range,
  list(tau_min=602,tau_max=1267)
)
# (9e) set_density_model

# First check the error messages. The order they are checked is not the same
# order they would be encountered in reading through the source code.

empty_density_model <- list()
bad_analysis_name <- "not_yet_defined"
expect_error(
  set_density_model(data_dir,bad_analysis_name,density_model),
  "A save file for analysis_name does not exist in data_dir"
)

expect_error(
  set_density_model(data_dir,analysis_name,empty_density_model),
  "type must be a named field in the list density_model"
)

invalid_density_model <- list(type="invalid_type")
expect_error(
  set_density_model(data_dir,analysis_name,invalid_density_model),
  "Unsupported type for density_model"
)

# no tau_min
invalid_density_model <- list(type="trunc_gauss_mix",
                              tau_max=1300,
                              K=4)
expect_error(
  set_density_model(data_dir,analysis_name,invalid_density_model),
  "tau_min must be a field in density_model for trunc_gauss_mix"
)

# no tau_max
invalid_density_model <- list(type="trunc_gauss_mix",
                              tau_min=600,
                              K=4)
expect_error(
  set_density_model(data_dir,analysis_name,invalid_density_model),
  "tau_max must be a field in density_model for trunc_gauss_mix"
)

# no K
invalid_density_model <- list(type="trunc_gauss_mix",
                              tau_min=600,
                              tau_max=1300)
expect_error(
  set_density_model(data_dir,analysis_name,invalid_density_model),
  "K must be a field in density_model for trunc_gauss_mix"
)

# tau_max less than tau_min
invalid_density_model <- list(type="trunc_gauss_mix",
                              tau_min=1300,
                              tau_max=600,
                              K=4)
expect_error(
  set_density_model(data_dir,analysis_name,invalid_density_model),
  "tau_min must be less than tau_max"
)

# Make a correct call to update the saved .rds file, check the result,  and make
# sure that an error is thrown if an attempt is made to reset the density_model.
density_model <- list(type="trunc_gauss_mix",
                              tau_min=600,
                              tau_max=1300,
                              K=4)
expect_error(
  set_density_model(data_dir,analysis_name,density_model),
  NA
)

expect_equal(
  readRDS(path_to_analysis_file),
  list(rc_meas=rc_meas,density_model=density_model)
)

expect_error(
  set_density_model(data_dir,analysis_name,density_model),
  "A density model has already been defined for this analysis"
)

# (9f) set_calib_curve
expect_error(
  set_calib_curve(data_dir,bad_analysis_name),
  "A save file for analysis_name does not exist in data_dir"
)

expect_error(
  set_calib_curve(data_dir,analysis_name,"intcal20"),
  NA
)

analysis <- readRDS(path_to_analysis_file)
expect_equal(
  names(analysis),
  c("rc_meas","density_model","calib_df")
)

calib_df <- load_calib_curve("intcal20")
expect_equal(
  analysis$calib_df,
  calib_df
)

# Directly delete calib_df and ensure from the .rds save file and check that it
# can also be set with a dataframe as input

analysis["calib_df"] <- NULL
saveRDS(analysis,path_to_analysis_file)
expect_error(
  set_calib_curve(data_dir,analysis_name,calib_df),
  NA
)

analysis <- readRDS(path_to_analysis_file)
expect_equal(
  names(analysis),
  c("rc_meas","density_model","calib_df")
)

calib_df <- load_calib_curve("intcal20")
expect_equal(
  analysis$calib_df,
  calib_df
)

# (9g) do_max_lik_fits
expect_error(
  do_max_lik_fits(data_dir,bad_analysis_name),
  "A save file for analysis_name does not exist in data_dir"
)

bad_analysis_name <- "bad"
path_to_bad_analysis_file <- file.path(data_dir,
                                       paste0(bad_analysis_name,".rds"))
saveRDS(list(),path_to_bad_analysis_file)
expect_error(
  do_max_lik_fits(data_dir,bad_analysis_name),
  "Radiocarbon measurements have not specified for this analysis"
)

saveRDS(list(rc_meas=rc_meas),path_to_bad_analysis_file)
expect_error(
  do_max_lik_fits(data_dir,bad_analysis_name),
  "A density model has not been specified for this analysis"
)

saveRDS(list(rc_meas=rc_meas,density_model=density_model),
        path_to_bad_analysis_file)
expect_error(
  do_max_lik_fits(data_dir,bad_analysis_name),
  "A calibration curve has not been specified for this analysis"
)
invalid_density_model <- list(type="invalid_type")
saveRDS(list(rc_meas=rc_meas,
             density_model=invalid_density_model,
             calib_df=calib_df),
        path_to_bad_analysis_file)
expect_error(
  do_max_lik_fits(data_dir,bad_analysis_name),
  "Unsupported type for density_model"
)

file.remove(path_to_bad_analysis_file)

expect_error(
  do_max_lik_fits(data_dir,analysis_name,num_restarts=4,maxfeval=200),
  NA
)

expect_error(
  do_max_lik_fits(data_dir,analysis_name),
  "Maximum likelihood fits have already been defined for this analysis"
)

# Check reproducbility with and without multiple cores
analysis0 <- readRDS(path_to_analysis_file)
input_seed <- analysis0$max_lik_fits_seeds$base_seed
analysis <- analysis0
# Undo the maximum likelihood fitting
analysis$max_lik_fits       <- NULL
analysis$max_lik_fits_seeds <- NULL
saveRDS(analysis,path_to_analysis_file)

expect_error(
  do_max_lik_fits(data_dir,analysis_name,num_restarts=4,maxfeval=200,
                  input_seed=input_seed),
  NA
)

# analysis0 and new_analysis should be identical except that
# analysis0$max_lik_fits_seeds$input_seed is NA. Check that it is indeed NA,
# set it to base_seed, and make sure that the resulting lists are identical.
expect_equal(
  is.na(analysis0$max_lik_fits_seeds$input_seed),
  TRUE
)

analysis0_modified <- analysis0
analysis0_modified$max_lik_fits_seeds$input_seed <- input_seed
expect_equal(
  readRDS(path_to_analysis_file),
  analysis0_modified
)

# Rewrite the file without the maximum likelihood fits, then check the
# reproducibility when using multiple cores.
saveRDS(analysis,path_to_analysis_file)

expect_error(
  do_max_lik_fits(data_dir,analysis_name,num_restarts=4,maxfeval=200,
                  input_seed=input_seed,num_cores=2),
  NA
)

expect_equal(
  readRDS(path_to_analysis_file),
  analysis0_modified
)

file.remove(path_to_analysis_file)