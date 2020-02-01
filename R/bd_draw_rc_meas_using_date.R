#' @title Simulate radiocarbon measurements from input calendar dates
#'
#' @description Given input calendar dates t_e, simulate the radiocarbon
#'              measurement process to yield a fraction modern and uncertainty.
#'              errorSpec is a list with named components min and max. The
#'              uncertainty, sig_m, is drawn uniformily from this range. By
#'              default, t_e is assumed to be calBP, but can be AD if isAD is
#'              True.
#'
#' @param t_e A vector of event dates (calBP by default)
#' @param calibDf The calibration dataframe, with columns yearBP, uncalYearBP, and uncalYearBPError
#' @param errorSpec A list with entries min and max specifying the interval for a uniform draw
#' @param isAD (Optional; default F) Whether t_e is calBP or AD

#' @export
bd_draw_rc_meas_using_date <- function(t_e, calibDf, errorSpec, isAD = F) {
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
  return(list(phi_m = phi_m, sig_m = sig_m, trc_m = trc_m, sig_trc_m = sig_trc_m))
}
