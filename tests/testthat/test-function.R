# Functional tests for baydem functions. These tests are slow and should
# primarily be run before pushing updates back to github.

# Use multiple cores to speed up processing. Always reserve two cores that are
# not used.
cores_to_use <- parallel::detectCores() - 2

# ------------------------------------------------------------------------------
# (1) Do functional tests for maximum likelihood fitting
#
# ------------------------------------------------------------------------------

# Simulate a (truncated) Gaussian mixture with N=20000 samples and ensure that
# the original simulation parameter vector is recovered with the maximum
# likelihood fit. Do this num_sims=4 times with a tolerance on the relative
# error of 2%.

N <- 20000
maximum_relative_error <- rep(.02,6)

th_sim <-
  c(
    pi1 = 0.4,
    pi2 = 0.6,
    mu1 = 1000,
    mu2 = 1250,
    s1 = 45,
    s2 = 35
  )

tau_min <- 700
tau_max <- 1500
dtau <- 1
tau <- seq(tau_min,tau_max,by=dtau)

calib_df <- load_calib_curve("intcal20")

# Explicitly create a vector of random seeds for the simulations
set.seed(135066)
num_sims <- 4
seed_vect <- sample.int(1000000,4)

for(n_s in 1:num_sims) {
  sim_spec <- list(model_spec=
                   list(density_type = "trunc_gauss_mix",
                        th=c(th_sim,tau_min,tau_max),
                        error_spec=list(type="unif_fm",min=.0021,max=.0028),
                        is_AD=T),
                   N=N,
                   calib_curve="intcal20",
                   seed=seed_vect[n_s])
  sim <- simulate_rc_data(sim_spec)
  expect_error(
    max_lik_fit <- fit_trunc_gauss_mix(2,
                                   sim$data$rc_meas$phi_m,
                                   sim$data$rc_meas$sig_m,
                                   tau_min,
                                   tau_max,
                                   1,
                                   calib_df,
                                   num_restarts=20,
                                   maxfeval=10000,
                                   num_cores=cores_to_use),
    NA
  )
  relative_error <- (max_lik_fit$th - th_sim)/th_sim
  expect_equal(
    all(abs(relative_error) <= maximum_relative_error),
    TRUE
  )
  # Make sure memory is free by doing garbage collection
  dummy <- gc(verbose=F)
}
