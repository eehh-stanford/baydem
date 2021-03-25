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

bd_simulate_rc_data <- function(sim_spec) {
  # Right now, two model types are supported for the target density:
  # (a) A non-truncated Gaussian mixture (gauss_mix)
  # (b) A     truncated Gaussian mixture (trunc_gauss_mix)

  # (1) If necessary, set the random number seed
  if( "seed" %in% names(sim_spec)) {
    set.seed(sim_spec$speed)
  }

  # (2) Make draws for the dates
  if (sim_spec$model_spec$density_type == "gauss_mix") {
    dates <- baydem::bd_sample_gauss_mix(sim_spec$N,
                                         sim_spec$model_spec$th)
  } else if (sim_spec$model_spec$density_type == "trunc_gauss_mix") {
    K <- (length(sim_spec$model_spec$th) - 2)/3
    dates <- baydem::bd_sample_gauss_mix(sim_spec$N,
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
  calib_df <- bd_load_calib_curve(sim_spec$calib_curve)
  rc_meas <-
    bd_draw_rc_meas_using_date(dates,
                               calib_df,
                               sim_spec$model_spec$error_spec,
                               sim_spec$model_spec$isAD)
  # rc_meas contains phi_m, sig_m, trc_m, and sig_trc_m, though phi_m and sig_m
  # uniquely determine trc_m and sig_trc_m.
  return(list(sim_spec=sim_spec,data=list(dates=dates,rc_meas=rc_meas)))
}