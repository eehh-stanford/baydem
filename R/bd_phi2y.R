#' @title Find calendar date from fraction modern value given known bounding indices
#'
#' @details
#' y_curve and phi_curve give the calendar date and fraction modern of the
#' radiocarbon calibration curve. It is known that the calendar date lies
#' lies between y_curve[ii_lo] and y_curve[ii_hi]. Interpolate to find the
#' calendar date corresponding to the input phi_known.
#'
#' @param y_curve Calibration curve calendar dates
#' @param phi_curve Calibration curve fraction moderns
#' @param phi_known Known fraction modern value
#' @param ii_lo Lower bounding index
#' @param ii_hi Upper bounding index
#'
#' @return The calendar date corresponding to the input phi_known
#'
#' @export


bd_phi2y <- function(y_curve, phi_curve, phi_known, ii_lo, ii_hi) {
  y <- y_curve[ii_lo] + (y_curve[ii_hi] - y_curve[ii_lo]) * (phi_known - phi_curve[ii_lo]) / (phi_curve[ii_hi] - phi_curve[ii_lo])
  return(y)
}
