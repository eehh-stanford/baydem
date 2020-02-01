#' @title Calculation the calibration curve fraction taumodern for input calendar dates
#'
#' @details
#' The input calibration data frame has three columns: yearBP, uncalYearBP, and
#' uncalYearBPError. For each input calendar date, tau, use the calibration curve
#' information to estimate the fraction modern of the calibration curve. If
#' tau is not specified, calibDf$uncalYearBP is used as the dates for the
#' calculation. By default, dates are assumed to be AD, but this can be changed
#' using the optional input isBP, which is FALSE by default.
#'
#' @param calibDf The calibration data frame, with columns yearBP, uncalYearBP, and uncalYearBPError
#' @param tau The calendar dates tau (if not input, calibDf$uncalYearBP is used)
#' @param isBP (default FALSE) Whether the input dates are before present (BP), as opposed to AD
#'
#' @return The vector of calibration curve fraction modern values
#' @export
bd_calc_calib_curve_frac_modern <- 
  function(calibDf, 
           tau = NA, 
           isBP = FALSE) {
  phi_curve <- exp(-calibDf$uncalYearBP / 8033)
  if (all(is.na(tau))) {
    return(phi_curve)
  }

  tau_curve <- 1950 - calibDf$yearBP
  if (isBP) {
    tau <- 1950 - tau # Convert to AD
  }

  phi <- stats::approx(tau_curve, phi_curve, tau)
  phi <- phi$y
  return(phi)
}
