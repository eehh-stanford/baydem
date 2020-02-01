#' @title Find calendar date from fraction modern value given known bounding indices
#'
#' @details
#' tau_curve and phi_curve give the calendar date and fraction modern of the
#' radiocarbon calibration curve. It is known that the calendar date lies
#' lies between tau_curve[ii_lo] and tau_curve[ii_hi]. Interpolate to find the
#' calendar date corresponding to the input phi_known. This function is called
#' by bd_assess_calib_curve_equif.
#'
#' @param tau_curve Calibration curve calendar dates
#' @param phi_curve Calibration curve fraction moderns
#' @param phi_known Known fraction modern value
#' @param ii_lo Lower bounding index
#' @param ii_hi Upper bounding index
#'
#' @return The calendar date corresponding to the input phi_known
#'
#' @export


bd_phi2tau <- function(tau_curve, phi_curve, phi_known, ii_lo, ii_hi) {
  tau <- tau_curve[ii_lo] + (tau_curve[ii_hi] - tau_curve[ii_lo]) * (phi_known - phi_curve[ii_lo]) / (phi_curve[ii_hi] - phi_curve[ii_lo])
  return(tau)
}
