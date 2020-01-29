#' @title Assess the equifinality or time spans in the input radiocarbon calibration curve
#'
#' @details
#' The input calibration data frame has three columns: yearBP, uncalYearBP, and
#' uncalYearBPError. A time span is equifinal if, for each date in the span,
#' there is at least one other date in the radiocarbon calibration curve with
#' the same fraction modern value. Conversely, a time span is not equifinal if
#' this is not true (all dates in the time span correspond to a unique fraction
#' modern value). Thus, the calibration curve is divided into alternating
#' equifinal and non-equifinal spans. This function identifies these regions
#' for the input calibration data frame, calibDf.
#'
#' Although calibDf uses year BP, all calculations and returned data use AD.
#'
#' A list is returned with two named variables: canInvert and invSpanList (inv =
#' invert, which is conceptually identical to non-equifinal).
#'
#' canInvert is a boolean that indicates whether each entry (row) in the
#' calibration data frame, calibDf, is inside or outside an invertible region.
#'
#' invSpanList is a list summarizing information about all the invertible (non-
#' equifinal) time spans. It contains the following named variables:
#'
#' ind       -- Indices (rows) of calibDf following within the invertible time
#'              span
#' y_left    -- Calendar date (AD) of the left (earlier) boundary of the
#'              invertible time span
#' phi_left  -- Fraction modern value of the left (earlier) boundary of the
#'              invertible time span
#' y_right   -- Calendar date (AD) of the right (later) boundary of the
#'              invertible time span
#' phi_right -- Fraction modern value of the right (later) boundary of the
#'              invertible time span
#' ii_prev   -- The index (row) of the calibDf earlier than the invertible
#'              region with closest fraction modern value
#' ii_next   -- The index (row) of the calibDf later than the invertible
#'
#' @param calibDf The calibration data frame, with columns yearBP, uncalYearBP, and uncalYearBPError
#' @param equiList (Optional) The result of a call to bd_calc_calib_curve_equif_dates
#'
#' @return A list of with named variables canInvert and invSpanlist (see details)
#'
#' @export
bd_assess_calib_curve_equif <- function(calibDf, equiList = NA) {
  if (all(is.na(equiList))) {
    equiList <- baydem::bd_calc_calib_curve_equif_dates(calibDf)
  }

  canInvert <- rep(T, length(calibDf$yearBP))
  y_curve <- 1950 - calibDf$yearBP
  phi_curve <- baydem::bd_calc_calib_curve_frac_modern(calibDf)
  for (equiEntry in equiList) {
    canInvert[equiEntry$indBase] <- F
  }
  clustList <- evd::clusters(canInvert, .5)
  invSpanList <- list()
  for (cc in 1:length(clustList)) {
    clust <- clustList[[cc]]
    lo <- as.numeric(names(clust)[1])
    hi <- as.numeric(names(clust)[length(clust)])
    # Find left boundary
    if (lo == 1) {
      y_left <- y_curve[1]
      phi_left <- phi_curve[1]
      ii_prev <- 1
    } else {
      ii_prev <- which.min(abs(phi_curve[1:(lo - 1)] - phi_curve[lo]))
      y_left <- baydem::bd_phi2y(y_curve, phi_curve, phi_curve[ii_prev], lo - 1, lo)
      phi_left <- phi_curve[ii_prev]
    }

    # Find right boundary
    if (hi == length(y_curve)) {
      y_right <- y_curve[length(y_curve)]
      phi_right <- phi_curve[length(y_curve)]
      ii_next <- length(y_curve)
    } else {
      ii_next <- hi + which.min(abs(phi_curve[(hi + 1):length(phi_curve)] - phi_curve[hi]))
      y_right <- baydem::bd_phi2y(y_curve, phi_curve, phi_curve[ii_next], hi, hi + 1)
      phi_right <- phi_curve[ii_next]
    }
    invSpan <- list(ind = lo:hi, y_left = y_left, phi_left = phi_left, y_right = y_right, phi_right = phi_right, ii_prev = ii_prev, ii_next = ii_next)
    invSpanList[[cc]] <- invSpan
  }
  return(list(invSpanList = invSpanList, canInvert = canInvert))
}
