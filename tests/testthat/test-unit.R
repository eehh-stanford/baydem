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
# (2) do unit tests for functions related to the maximum likelihood fitting
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
print(names(sim))
expect_error(
  max_lik_fit <-
    temper_trunc_gauss_mix(
                           2,
                           sim$data$rc_meas$phi_m,
                           sim$data$rc_meas$sig_m,
                           600,
                           1300,
                           1,
                           calib_df,
                           samps_per_cyc=5,
                           num_cyc=10),
  NA
)

expect_equal(
  names(max_lik_fit),
  c("th","neg_log_lik","tau","f","bic","aic")
)

expect_error(
  max_lik_fit <-
    temper_trunc_gauss_mix(2,
                           sim$data$rc_meas$phi_m,
                           sim$data$rc_meas$sig_m,
                           600,
                           1300,
                           1,
                           calib_df,
                           samps_per_cyc=5,
                           num_cyc=10,
                           num_cores=2),
  NA
)

expect_equal(
  names(max_lik_fit),
  c("th","neg_log_lik","tau","f","bic","aic")
)