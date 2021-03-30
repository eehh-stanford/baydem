#' Create simulated radiocarbon data
#'
#' The input is a multi-tiered list specifying all aspects of the simulation. In
#' particular, sim_spec (for simulation specification) consists of three
#' required fields and one optional field: model_spec (model specification),
#' N (number of samples), calib_curve (calibration curve), and seed (an optional
#' random number seed). model_spec, in turn, is a list than consists of three
#' required fields: density_type (the type of parametric model to use for the
#' target density), th (the parameter vector for the target density), isAD (a
#' boolean variable indicating whether the target density is for years AD or
#' BP), and error_spec (a specification for how to model the measurement
#' errors). Currently, two model types are supported for the target density,
#' gauss_mix = a Gaussian mixture and trunc_gauss_mix = a truncated Gaussian
#' mixture. Currently, one specfiaction of the error_spec (a list) is supported:
#' unif_fm = a uniform draw for the error of the fraction modern value.
#' For unif_fm, error_spec must contain the fields min and max (the minimum and
#' maximum values for the uniform draw).
#'
#' The output is a list with the fields dates and rc_meas. Dates is a length N
#' vector of dates (on non-simulated data, cannot be directly known). rc_meas
#' is a list consisting of four length N vectors: phi_m, the fracture modern,
#' sig_m, the uncertainty of the fraction modern, trc_m, trc_m, the measurement
#' in uncalibrated radiocarbon years, and sig_trc_m, the error for trc_m. The
#' following hierarchical summary may be easier to digest:
#'
#' sim_spec            A full specification of the simulation
#'     model_spec          A specification for the model used to generate the
#'                          data
#'         density_type    The model type to use for the target density (e.g., a
#'                          a truncated Gaussian mixture)
#'         th              The parameter vector for the target density
#'         isAD            A boolean indicating whether the parameter vector
#'                          assumes AD or BP (AD is 1950 - BP, which allows
#'                          negative dates)
#'         error_spec      A specification for the measurement errors
#'     N               The number of random samples to make
#'     calib_curve     The calibration curve (currently, must be a named intcal
#'                      curve, but support for arbitrary curves could be added
#'                      in the future)
#'     seed            An optional random number seed to use to ensure
#'                      reproducibility
#'
#' data                The simulated data
#'     dates               A vector of original dates for the samples
#'     rc_meas             A list with the radiocarbon measurements
#'         phi_m               A vector of simulated fraction moderns
#'         sig_m               A vector of simulated standard deviations for
#'                              phi_m
#'         trc_m               A vector of simulated uncalibrated radiocarbon
#'                              years (BP)
#'         sig_trc_m           A vector of simulated uncertainties for trc_m
#'
#' @export
simulate_rc_data <- function(sim_spec) {
  # Right now, two model types are supported for the target density:
  # (a) A non-truncated Gaussian mixture (gauss_mix)
  # (b) A     truncated Gaussian mixture (trunc_gauss_mix)

  # (1) If necessary, set the random number seed
  if( "seed" %in% names(sim_spec)) {
    set.seed(sim_spec$seed)
  }

  # (2) Make draws for the dates
  if (sim_spec$model_spec$density_type == "gauss_mix") {
    dates <- baydem::sample_gauss_mix(sim_spec$N,
                                         sim_spec$model_spec$th)
  } else if (sim_spec$model_spec$density_type == "trunc_gauss_mix") {
    K <- (length(sim_spec$model_spec$th) - 2)/3
    dates <- baydem::sample_gauss_mix(sim_spec$N,
                                         sim_spec$model_spec$th[1:(3*K)],
                                         taumin=sim_spec$model_spec$th[3*K+1],
                                         taumax=sim_spec$model_spec$th[3*K+2])
  } else {
    error_message <- paste0("Unsupported density_type = ",
                            sim_spec$model_spec$density_type)
    stop(error_message)
  }

  # (2) Make draws for the radiocarbon measurements using the already sampled
  #     dates
  # TODO: consider adding an error check of the input calibration curve
  # TODO: consider allowing arbitrary calibration curves to be used
  calib_df <- load_calib_curve(sim_spec$calib_curve)
  rc_meas <-
    draw_rc_meas_using_date(dates,
                            calib_df,
                            sim_spec$model_spec$error_spec,
                            sim_spec$model_spec$isAD)
  # rc_meas contains phi_m, sig_m, trc_m, and sig_trc_m, though phi_m and sig_m
  # uniquely determine trc_m and sig_trc_m.
  return(list(sim_spec=sim_spec,data=list(dates=dates,rc_meas=rc_meas)))
}

#' @title Sample from a possibly truncated Gaussian mixture
#'
#' @details `N` is the number of samples to draw and `th` specifies the Gaussian
#' mixture with the ordering (pi_1,...,pi_K,mu_1,...,mu_K,sig_1,...,sig_K),
#' where `K` is the number of mixture components. Optionally, the samples are
#' drawn on the truncated interval `taumin` to `taumax`. Because of limitations
#' in the package `distr`, the maximum number of mixture components suppored is
#' K=4.
#'
#' @param N Number of samples
#' @param th Parameterization vector for the Gaussian mixture
#' @param taumin (Optional) Lower bound for samples
#' @param taumax (Optional) Upper bound for samples
#'
#' @return N samples from the Gaussian mixture
#'
#' @export
sample_gauss_mix <- function(N, th, taumin = NA, taumax = NA) {
  K <- length(th) / 3 # Number of  mixtures

  # (Somewhat awkwardly), define the mixing distribution for K = 1
  # to 4 directly. This is necessary because truncate expects an
  # AbscontDistribution but creating a mixture distribution via vectors
  # creates a UnivarMixingDistribution, which cannot be cast to an
  # AbscontDistribution. Perhaps this shortcoming of distr will be addressed
  # in the future. In the meantime, the number of mixtures is limited to 4
  if (K == 1) { # Allow K = 1 (i.e., a Gaussian, not a mixture) as a special case
    normMix <- distr::Norm(mean = th[2], sd = th[3])
  } else if (K == 2) {
    normMix <- distr::UnivarMixingDistribution(distr::Norm(mean = th[3], sd = th[5]),
      distr::Norm(mean = th[4], sd = th[6]),
      mixCoeff = th[1:2]
    )
  } else if (K == 3) {
    normMix <- distr::UnivarMixingDistribution(distr::Norm(mean = th[4], sd = th[7]),
      distr::Norm(mean = th[5], sd = th[8]),
      distr::Norm(mean = th[6], sd = th[9]),
      mixCoeff = th[1:3]
    )
  } else if (K == 4) {
    normMix <- distr::UnivarMixingDistribution(distr::Norm(mean = th[5], sd = th[9]),
      distr::Norm(mean = th[6], sd = th[10]),
      distr::Norm(mean = th[7], sd = th[11]),
      distr::Norm(mean = th[8], sd = th[12]),
      mixCoeff = th[1:4]
    )
  } else {
    stop("The maximum number of supported mixture components is 4")
  }

  # Use distr to sample from a two-component, truncated Gaussian mixture
  if (!is.na(taumin) && !is.na(taumax)) {
    normMixTrunc <- distr::Truncate(normMix, taumin, taumax)
    samp <- distr::r(normMixTrunc)(N)
  } else {
    samp <- distr::r(normMix)(N)
  }
  return(samp)
}

#' Simulate radiocarbon measurements from input calendar dates
#'
#' Given input calendar dates t_e, simulate the radiocarbon measurement process
#' to yield a fraction modern and uncertainty. error_spec is a list with named
#' type and additional fields specifying the parameters. Currently only type
#' unif_fm (for uniform fraction modern) is supported, which makes a uniform
#' draw on the interval min to max for the uncertainty of the fraction modern
#' value, where min and max are fields in error_spec giving the boundaries of
#' the uniform draw. By default, t_e is assumed to be calBP, but can be AD (more
#' more precisely, 1950 - BP value) if isAD is True.
#'
#' @param t_e A vector of event dates (calBP by default)
#' @param calibDf The calibration dataframe, with columns yearBP, uncalYearBP,
#'   and uncalYearBPError
#' @param errorSpec A list with entries type and parameter values specifying how
#'  to add uncertainty to the simualted measurement. Currently, only type =
#'  "unif_fm" is supported, for which the parameters min and max must be in
#'  errorSpec.
#' @param isAD (Optional; default F) Whether t_e is calBP or AD

#' @export
draw_rc_meas_using_date <- function(t_e, calibDf, errorSpec, isAD = F) {

  # Check that the errorSpec is of a know type
  if (errorSpec$type != "unif_fm") {
    error_message <- paste0("Unrecognized error type = ",errorSpec$type)
    stop(error_message)
  }

  # If input calendar dates are AD, convert to calBP
  if (isAD) {
    t_e <- 1950 - t_e
  }

  N <- length(t_e)

  tau_curve <- rev(calibDf$yearBP)
  mu_k_curve <- exp(-rev(calibDf$uncalYearBP) / 8033)
  sig_k_curve <- mu_k_curve * rev(calibDf$uncalYearBPError) / 8033

  # Interpolate curves at t_e to yield mu_k
  mu_k <- stats::approx(tau_curve, mu_k_curve, t_e)
  mu_k <- mu_k$y

  # Interpolate curves at t_e to yield sig_k
  sig_k <- stats::approx(tau_curve, sig_k_curve, t_e)
  sig_k <- sig_k$y

  # Sample the measurement errors
  sig_m <- stats::runif(N, errorSpec$min, errorSpec$max)
  sig_tot <- sqrt(sig_m^2 + sig_k^2)

  # The measured "ratios"
  phi_m <- stats::rnorm(N, mu_k, sig_tot)

  # Calculate radiocarbon years (uncal) measurement and error
  sig_trc_m <- 8033 * sig_m / phi_m
  trc_m <- -8033 * log(phi_m)
  return(list(phi_m=phi_m,sig_m=sig_m,trc_m=trc_m,sig_trc_m=sig_trc_m))
}
