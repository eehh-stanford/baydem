# Unit tests for baydem functions. These tests are fast but do not check
# function. Thus, they can be quickly run as code is actively developed. The
# function tests are located in test-function.R. They take much longer to run
# and should primarily be run before pushing updates back to github.

# ------------------------------------------------------------------------------
# (1) Do unit tests for functions that generate simulated data
#
# Coverage: bd_sample_gauss_mix
#           bd_draw_rc_meas_using_date
#           bd_simulate_rc_data
# ------------------------------------------------------------------------------

# (1a) test bd_sample_gauss_mix
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
  sample <- bd_sample_gauss_mix(N,
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
  sample <- bd_sample_gauss_mix(N,
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
  bd_sample_gauss_mix(N,1:15),
  "The maximum number of supported mixture components is 4"
)

# (1b) test bd_draw_rc_meas_using_date
calib_df <- bd_load_calib_curve("intcal20")
# test with isAD FALSE
expect_error(
  rc_meas <- bd_draw_rc_meas_using_date(sample,
                                        calib_df,
                                        list(type="unif_fm",
                                             min=.0021,
                                             max=.0028),
                                        isAD=F),
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
  rc_meas <- bd_draw_rc_meas_using_date(sample,
                                        calib_df,
                                        list(type="not_valid",
                                             min=.0021,
                                             max=.0028),
                                        isAD=F),
  "Unrecognized error type = not_valid"
)

# test with isAD TRUE
expect_error(
  rc_meas <- bd_draw_rc_meas_using_date(sample,
                                        calib_df,
                                        list(type="unif_fm",
                                             min=.0021,
                                             max=.0028),
                                        isAD=T),
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
  rc_meas <- bd_draw_rc_meas_using_date(sample,
                                        calib_df,
                                        list(type="not_valid",
                                             min=.0021,
                                             max=.0028),
                                        isAD=T),
  "Unrecognized error type = not_valid"
)

# (1c) test bd_simulate_rc_data

# Check simulation of a non-truncated Gaussian mixture with no random number
# seed set
sim_spec <- list(model_spec=
                   list(density_type = "gauss_mix",
                        th=th_sim,
                        error_spec=list(type="unif_fm",min=.0021,max=.0028),
                        isAD=T),
                 N=N,
                 calib_curve="intcal20")

expect_error(
  sim <- bd_simulate_rc_data(sim_spec),
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
                        isAD=T),
                 N=N,
                 calib_curve="intcal20")

expect_error(
  sim <- bd_simulate_rc_data(sim_spec),
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
                        isAD=T),
                 N=N,
                 calib_curve="intcal20",
                 seed=1002)

expect_error(
  sim <- bd_simulate_rc_data(sim_spec),
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
  is_th_reduced_valid(c(0.5,10,20,2,3),sig_min=10),
  FALSE
)

# (2b) calc_trunc_gauss_mix_neg_log_lik
tau <- seq(tau_min,tau_max,by=5)
M <- bd_calc_meas_matrix(tau,
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
  calc_trunc_gauss_mix_neg_log_lik(th_sim[-1],M,tau,sig_min=100),
  Inf
)

# (2c) temper_trunc_gauss_mix
# Check that temper_trunc_gauss_mix does not throw an error (and also that the
# output has the right fields). Do this for both not using and using multiple
# cores.
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
  c("th","neg_log_lik","tau","f","bic","aic","hjkb_fit_list")
)

expect_error(
  max_lik_fit <-
    fit_trunc_gauss_mix(2,
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
  c("th","neg_log_lik","tau","f","bic","aic","hjkb_fit_list")
)

# ------------------------------------------------------------------------------
# (3) Do unit tests for functions related to the Bayesian inference
#
# coverage: bd_do_inference
#           bd_analyze_soln
# ------------------------------------------------------------------------------

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
    taumin = tau_min,
    # Maximum calendar date (years BC/AD)
    taumax = tau_max,
    # Spacing for the measurement matrix (years)
    dtau = 5,
    # Number of mixtures
    K = 2
  )

# Define a problem object to be input for the inference
prob <- list(
  phi_m = sim$data$rc_meas$phi_m,
  sig_m = sim$data$rc_meas$sig_m,
  hp = hp,
  calibDf = calib_df,
  control = list(
    sampsPerChain = 400,
    warmup = 200,
    stanControl = list(adapt_delta = .99)
  )
)

# (3a) bd_do_inference
expect_error(
  soln <- bd_do_inference(prob),
  NA
)

expect_equal(
  names(soln),
  c("prob","fit","control")
)

expect_equal(
  soln$prob,
  prob
)

expect_is(
  soln$fit,
  "stanfit"
)

expect_equal(
  names(soln$control),
  c("numChains","sampsPerChain","warmup","stanControl","initList")
)

# (3b) bd_analyze_soln
expect_error(
  anal <- bd_analyze_soln(soln),
  NA
)

expect_equal(
  names(anal),
  c("tau","f_spdf","Qdens","Qrate","probs",
    "rateProp","rateInd","growthState","dtau","summList")
)

# ------------------------------------------------------------------------------
# (4) Do unit tests for plotting functions
#
# coverage: bd_make_blank_density_plot
#           bd_plot_50_percent_quantile
#           bd_add_shaded_quantiles
#           bd_plot_summed_density
#           bd_plot_known_sim_density
# ------------------------------------------------------------------------------

# (4a) bd_make_blank_density_plot
expect_error(
  bd_make_blank_density_plot(anal),
  NA
)

# (4b) bd_plot_50_percent_quantile
expect_error(
  bd_plot_50_percent_quantile(anal, add = T),
  NA
)

# (4c) bd_add_shaded_quantiles
expect_error(
  bd_add_shaded_quantiles(anal),
  NA
)

# (4d) bd_plot_known_sim_density
#      bd_plot_summed_density
expect_error(
  bd_plot_known_sim_density(anal,add=T),
  NA
)

expect_error(
  bd_plot_summed_density(anal,add=T),
  NA
)

expect_error(
  bd_plot_known_sim_density(anal),
  NA
)

expect_error(
  bd_plot_summed_density(anal),
  NA
)

# (4e) bd_vis_calib_curve
expect_error(
  bd_vis_calib_curve(600, 700, calib_df),
  NA
)

# ------------------------------------------------------------------------------
# (5) Do unit tests for density calculations and related "helper functions"
#
# coverage: bd_extract_param
#           bd_calc_gauss_mix_pdf
#           bd_calc_gauss_mix_pdf_mat
#           bd_calc_quantiles
#           bd_calc_half_life_from_peak
#           bd_calc_relative_density
#             unpack_spec           [indirect]
#             bd_calc_point_density [indirect]
#             bd_calc_range_density [indirect]
#             bd_calc_peak_density  [indirect]
#           bd_summarize_trunc_gauss_mix_sample
# ------------------------------------------------------------------------------

# (5a) bd_extract_param
total_samples <- 800
expect_error(
  TH <- bd_extract_param(soln$fit),
  NA
)

expect_equal(
  dim(TH),
  c(total_samples,6)
)

# (5b) bd_calc_gauss_mix_pdf
tau <- seq(hp$taumin,hp$taumax,by=hp$dtau)

# density is the default type
expect_error(
  pdf_vect <- bd_calc_gauss_mix_pdf(TH[50, ],tau),
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
  cdf_vect <- bd_calc_gauss_mix_pdf(TH[50, ], tau, type = "cumulative"),
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
  deriv_vect <- bd_calc_gauss_mix_pdf(TH[50, ], tau, type = "derivative"),
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
  rate_vect <- bd_calc_gauss_mix_pdf(TH[50, ], tau, type = "rate"),
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

# (5c) bd_calc_gauss_mix_pdf_mat

# density is the default type
expect_error(
  pdf_mat <- bd_calc_gauss_mix_pdf_mat(TH,tau),
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
  cdf_mat <- bd_calc_gauss_mix_pdf_mat(TH, tau, type = "cumulative"),
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
  deriv_mat <- bd_calc_gauss_mix_pdf_mat(TH, tau, type = "derivative"),
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
  rate_mat <- bd_calc_gauss_mix_pdf_mat(TH, tau, type = "rate"),
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

# (5d) bd_calc_quantiles

expect_error(
  Qdens1 <- bd_calc_quantiles(pdf_mat),
  NA
)

expect_equal(
  dim(Qdens1),
  c(3, ncol(pdf_mat))
)

expect_error(
  Qdens2 <- bd_calc_quantiles(pdf_mat, c(.025, .05, .5, .90, .975)),
  NA
)
expect_equal(
  dim(Qdens2),
  c(5, ncol(pdf_mat))
)

# (5e) bd_calc_half_life_from_peak
expect_error(
  half_life1 <- bd_calc_half_life_from_peak(soln),
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
  half_life2 <- bd_calc_half_life_from_peak(soln, propChange = .25),
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

# (5f) bd_calc_relative_density
#        unpack_spec           [indirect]
#        bd_calc_point_density [indirect]
#        bd_calc_range_density [indirect]
#        bd_calc_peak_density  [indirect]
# Check two calls to bd_calc_relative_density to check all four helper
# functions (which are only checked indirectly).

expect_error(
  rel_dens1 <- bd_calc_relative_density(soln, "peak", 1100),
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
  rel_dens2 <- bd_calc_relative_density(soln, 900, c(700, 750)),
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

# (5g) bd_summarize_trunc_gauss_mix_sample
expect_error(
  summ <- bd_summarize_trunc_gauss_mix_sample(TH[50, ], hp$taumin, hp$taumax),
  NA
)

expect_equal(
  names(summ),
  c("periods","indPeak","tpeak","fpeak","pattern")
)

# ------------------------------------------------------------------------------
# (6) Do unit tests for the measurement matrix calculation
#
# coverage: bd_calc_meas_matrix
#           bd_calc_trapez_weights
# ------------------------------------------------------------------------------


# (1) Calculate the measurement matrix using tau, already defined above, which
#     is regularly spaced. Do this for both using and not using calibration
#     uncertainty. Also check the dimensions of the output

# (6a) bd_calc_meas_matrix
#        tau regularly spaced
#        with calibration uncertainty
expect_error(
  M6a <- bd_calc_meas_matrix(tau,prob$phi_m,prob$sig_m,calib_df,addCalibUnc=T),
  NA
)

expect_equal(
  dim(M6a),
  c(length(prob$phi_m), length(tau)),
)

expect_equal(
  any(is.na((M6a))),
  FALSE
)

# (6b) bd_calc_meas_matrix
#        tau regularly spaced
#        without calibration uncertainty
expect_error(
  M6b <- bd_calc_meas_matrix(tau,prob$phi_m,prob$sig_m,calib_df,addCalibUnc=F),
  NA
)

expect_equal(
  dim(M6b),
  c(length(prob$phi_m), length(tau)),
)

expect_equal(
  any(is.na((M6b))),
  FALSE
)

# (6c) bd_calc_meas_matrix
#        tau irrregularly spaced
#        with calibration uncertainty
tau_irreg <- c(800, 805, 810, 825)

# an error is expected if useTrapez is not set to TRUE
expect_error(
  bd_calc_meas_matrix(tau_irreg,prob$phi_m,prob$sig_m,calib_df,addCalibUnc=T),
  "tau is irregularly spaced but useTrapez is FALSE"
)

expect_error(
  M6c <- bd_calc_meas_matrix(
    tau_irreg,
    prob$phi_m,
    prob$sig_m,
    calib_df,
    addCalibUnc=T,
    useTrapez=T),
  NA
)

expect_equal(
  dim(M6c),
  c(length(prob$phi_m), length(tau_irreg)),
)

expect_equal(
  any(is.na((M6c))),
  FALSE
)

# (6d) bd_calc_meas_matrix
#        tau irrregularly spaced
#        without calibration uncertainty

# an error is expected if useTrapez is not set to TRUE
expect_error(
  bd_calc_meas_matrix(tau_irreg,prob$phi_m,prob$sig_m,calib_df,addCalibUnc=F),
  "tau is irregularly spaced but useTrapez is FALSE"
)

expect_error(
  M6d <- bd_calc_meas_matrix(
    tau_irreg,
    prob$phi_m,
    prob$sig_m,
    calib_df,
    addCalibUnc=F,
    useTrapez=T),
  NA
)

expect_equal(
  dim(M6d),
  c(length(prob$phi_m), length(tau_irreg)),
)

expect_equal(
  any(is.na((M6d))),
  FALSE
)

# (6e) bd_calc_trapez_weights
expect_equal(
  bd_calc_trapez_weights(c(-1.5, 2, 3, 4, 7)),
  c(1.75, 2.25, 1, 2, 1.5)
)

# ------------------------------------------------------------------------------
# (7) Do unit tests for identifiability functions
#
# coverage: bd_assess_calib_curve_equif
#             bd_phi2tau [indirect]
#           bd_calc_calib_curve_equif_dates
#           bd_calc_calib_curve_frac_modern
# ------------------------------------------------------------------------------

# (7a) bd_assess_calib_curve_equif
#        bd_phi2tau [indirect]
#      bd_calc_calib_curve_equif_dates
expect_error(
  equif_result1 <- bd_assess_calib_curve_equif(calib_df),
  NA
)

expect_equal(
  names(equif_result1),
  c("invSpanList","canInvert")
)

expect_equal(
  length(equif_result1$canInvert),
  nrow(calib_df)
)

expect_error(
  equif_dates <- bd_calc_calib_curve_equif_dates(calib_df),
  NA
)

expect_error(
  equif_result2 <- bd_assess_calib_curve_equif(calib_df,equif_dates),
  NA
)

expect_equal(
  names(equif_result2),
  c("invSpanList","canInvert")
)

expect_equal(
  length(equif_result2$canInvert),
  nrow(calib_df)
)

expect_equal(
  equif_result1,
  equif_result2
)

# (7b) bd_calc_calib_curve_frac_modern

expect_error(
  phi <- bd_calc_calib_curve_frac_modern(calib_df),
  NA
)

expect_equal(
  phi,
  exp(-calib_df$uncalYearBP / 8033)
)

tau2 <- c(600, 602, 805.89)
expect_error(
  phi2 <- bd_calc_calib_curve_frac_modern(calib_df,tau2),
  NA
)

expect_equal(
  length(phi2),
  length(tau2)
)

# ------------------------------------------------------------------------------
# (8) Do unit tests for truncated exponential model
#
# coverage: bd_assess_calib_curve_equif
# ------------------------------------------------------------------------------
expect_error(
  exp_samp1 <- bd_sample_trunc_exp(50, 0.01, 600, 1300),
  NA
)

expect_equal(
  length(exp_samp1),
  50
)

expect_error(
  exp_samp2 <- bd_sample_trunc_exp(50, -0.01, 600, 1300),
  NA
)

expect_equal(
  length(exp_samp2),
  50
)
