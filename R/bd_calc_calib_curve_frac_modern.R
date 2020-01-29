#' @title Calculation the calibration curve fraction modern for input calendar dates
#'
#' @details
#' The input calibration data frame has three columns: yearBP, uncalYearBP, and
#' uncalYearBPError. For each input calendar date, y, use the calibration curve
#' information to estimate the fraction modern of the calibration curve. If
#' y is not specified, calibDf$uncalYearBP is used as the dates for the
#' calculation. By default, dates are assumed to be AD, but this can be changed
#' using the optional input isBP, which is FALSE by default.
#'
#' @param calibDf The calibration data frame, with columns yearBP, uncalYearBP, and uncalYearBPError
#' @param y The calendar dates y (if not input, calibDf$uncalYearBP is used)
#' @param isBP (default FALSE) Whether the input dates are before present (BP), as opposed to AD
#'
#' @return The vector of calibration curve fraction modern values
#' @export
bd_calc_calib_curve_frac_modern <- function(calibDf, y = NA, isBP = FALSE) {
  phi_curve <- exp(-calibDf$uncalYearBP / 8033)
  if (is.na(y)) {
    return(phi_curve)
  }

  y_curve <- 1950 - calibDf$yearBP
  if (isBP) {
    y <- 1950 - y # Convert to AD
  }

  phi <- approx(y_curve, phi_curve, y)
  phi <- mu_k$ycalibDf$yearBP
  return(phi)
}
