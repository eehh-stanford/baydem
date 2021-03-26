# Functional tests for baydem functions. These tests are slow and should
# primarily be run before pushing updates back to github.

# ------------------------------------------------------------------------------
# (1) Do functional tests for maximum likelihood fitting
#
# ------------------------------------------------------------------------------

# Simulate a (truncated) Gaussian mixture with N=100000 samples and ensure that
# the original simulation parameter vector is recovered with the maximum
# likelihood fit.

N <- 100000

th_sim <-
  c(
    pi1 = 0.4,
    pi2 = 0.6,
    mu1 = 775,
    mu2 = 1000,
    sig1 = 35,
    sig2 = 45
  )


tau_min <- 600
tau_max <- 1300

calib_df <- bd_load_calib_curve("intcal20")

sim_spec <- list(model_spec=
                   list(density_type = "trunc_gauss_mix",
                        th=c(th_sim,tau_min,tau_max),
                        error_spec=list(type="unif_fm",min=.0021,max=.0028),
                        isAD=T),
                 N=N,
                 calib_curve="intcal20",
                 seed=135066)

sim <- bd_simulate_rc_data(sim_spec)

# Set the random number seed so that the tempering is reproducible. Also, do
# not use multiple cores.
set.seed(250070)

max_lik_fit <- temper_trunc_gauss_mix(2,
                                      sim$data$rc_meas$phi_m,
                                      sim$data$rc_meas$sig_m,
                                      tau_min,
                                      tau_max,
                                      1,
                                      calib_df,
                                      samps_per_cyc=20,
                                      num_cyc=1000)

relative_error <- (max_lik_fit$th - th_sim)/th_sim

# Require the relative error to be within 1% for the mixture proportions and
# mean and 5% for the standard deviations (which are harder to estimate). s1 is
# especially hard to estimate due to the shape of the radiocarbon calibration
# curve.
maximum_error <- c(rep(.01,4),rep(.05,2))

expect_equal(
  all(abs(relative_error) <= maximum_error),
  TRUE
)
