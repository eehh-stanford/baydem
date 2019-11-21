#' @title Simulate a radiocabron measurement from input calendar dates
#'
#' @description Given an input calendar dates, simulate the radiocarbon
#'              measurement process to yield a fraction modern and uncertainty.
#'              errorSpec is a list with named components min and max. The
#'              uncertainty, sig_m, is drawn uniformily from this range.
#'
#' @export

bd_draw_rc_meas_using_date <- function(y_e, calibDf, errorSpec, isAD = F) {
  # If input calendar dates are AD, convert to calBP
  if (isAD) {
    y_e <- 1950 - y_e
  }

  N <- length(y_e)

  y_curve <- rev(calibDf$yearBP)
  mu_k_curve <- exp(-rev(calibDf$uncalYearBP) / 8033)
  sig_k_curve <- mu_k_curve * rev(calibDf$uncalYearBPError) / 8033

  # Interpolate curves at y_e to yield mu_k
  mu_k <- stats::approx(y_curve, mu_k_curve, y_e)
  mu_k <- mu_k$y

  # Interpolate curves at y_e to yield sig_k
  sig_k <- stats::approx(y_curve, sig_k_curve, y_e)
  sig_k <- sig_k$y

  # Sample the measurement errors
  sig_m <- stats::runif(N, errorSpec$min, errorSpec$max)
  sig_tot <- sqrt(sig_m^2 + sig_k^2)

  # The measured "ratios"
  phi_m <- stats::rnorm(N, mu_k, sig_tot)

  # Calculate radiocarbon years (uncal) measurement and error
  sig_yrc_m <- 8033 * sig_m / phi_m
  yrc_m <- -8033 * log(phi_m)
  return(list(phi_m = phi_m, sig_m = sig_m, yrc_m = yrc_m, sig_yrc_m = sig_yrc_m))
}
