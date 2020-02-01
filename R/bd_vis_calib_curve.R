#' @title Visualize the calibration curve with equifinal and non-equifinal time spans
#'
#' @details
#' Vizualize the input calibration curve, calibDf, on the interval taumin to
#' taumax. This involves plotting the curve itself and shading invertible
#' (non-equifinal) and non-invertible (equifinal) regions.
#'
#' @param taumin The minimum calendar date for plotting (AD)
#' @param taumax The maximum calendar date for plotting (AD)
#' @param calibDf The calibration data frame, with columns yearBP, uncalYearBP, and uncalYearBPError
#' @param invertCol (default: `gray90`) The color for shading invertible regions
#' @param pointCol (default: `black`) The color for calibration curve points
#' @param pointPch (default: `19`) The symbol for calibration curve points
#' @param ... Additional inputs to plots
#'
#' @export
bd_vis_calib_curve <-
  function(taumin,
           taumax,
           calibDf,
           invertCol = "gray90",
           pointCol = "black",
           pointPch = 19,
           ...) {
    tau_curve <- 1950 - calibDf$yearBP
    phi_curve <- exp(-calibDf$uncalYearBP / 8033)
    ind <- (tau_curve >= taumin) & (tau_curve <= taumax)
    tauVect <- tau_curve[ind]
    phiVect <- phi_curve[ind]
    phiMin <- min(phiVect)
    phiMax <- max(phiVect)

    # Create an empty plot
    graphics::plot(1, type = "n", xlim = c(taumin, taumax), ylim = c(phiMin, phiMax), ...)

    equiInfo <- baydem::bd_assess_calib_curve_equif(calibDf)
    canInvert <- equiInfo$canInvert
    invSpanList <- equiInfo$invSpanList
    for (ii in 1:length(invSpanList)) {
      invSpan <- invSpanList[[ii]]
      if (dplyr::between(invSpan$tau_left, taumin, taumax) || dplyr::between(invSpan$tau_right, taumin, taumax)) {
        graphics::rect(invSpan$tau_left, phiMin, invSpan$tau_right, phiMax, border = NA, col = invertCol)
      }
    }

    # Draw calibration curve points
    graphics::points(tauVect, phiVect, col = pointCol, pch = pointPch)
  }
