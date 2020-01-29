#' @title Visualize the calibration curve with equifinal and non-equifinal time spans
#'
#' @details
#' Vizualize the input calibration curve, calibDf, on the interval ymin to
#' ymax. This involves plotting the curve itself and shading invertible
#' (non-equifinal) and non-invertible (equifinal) regions.
#'
#' @param ymin The minimum calendar date for plotting (AD)
#' @param ymax The maximum calendar date for plotting (AD)
#' @param calibDf The calibration data frame, with columns yearBP, uncalYearBP, and uncalYearBPError
#' @param invertCol [Default gray90] The color for shading invertible regions
#' @param pointCol [Default black] The color for calibration curve points
#' @param ... Additional inputs to plots
#'
#' @export
bd_vis_calib_curve <- function(ymin, ymax, calibDf, invertCol = "gray90", pointCol = "black", pointPch = 19, ...) {
  y_curve <- 1950 - calibDf$yearBP
  phi_curve <- exp(-calibDf$uncalYearBP / 8033)
  ind <- (y_curve >= ymin) & (y_curve <= ymax)
  yVect <- y_curve[ind]
  phiVect <- phi_curve[ind]
  phiMin <- min(phiVect)
  phiMax <- max(phiVect)

  # Create an empty plot
  plot(1, type = "n", xlim = c(ymin, ymax), ylim = c(phiMin, phiMax), ...)
  # plot(1, type="n",xlim=c(ymin,ymax),ylim=c(phiMin,phiMax),xlab='Calendar Date [AD]',ylab='Fraction Modern',...)

  equiInfo <- baydem::bd_assess_calib_curve_equif(calibDf)
  canInvert <- equiInfo$canInvert
  invSpanList <- equiInfo$invSpanList
  for (ii in 1:length(invSpanList)) {
    invSpan <- invSpanList[[ii]]
    if (dplyr::between(invSpan$y_left, ymin, ymax) || dplyr::between(invSpan$y_right, ymin, ymax)) {
      rect(invSpan$y_left, phiMin, invSpan$y_right, phiMax, border = NA, col = invertCol)
    }
  }

  # Draw calibration curve points
  points(yVect, phiVect, col = pointCol, pch = pointPch)
}
