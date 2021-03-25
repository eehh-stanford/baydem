#' @title Calculate the measurement matrix given N observations and G grid points
#'
#' @description tau is a vector of calendar dates indexed from g = 1,2,...,G.
#'              phi_m is a vector of n = 1,2,...,N radiocarbon determinations
#'              (fraction modern) with associated uncertainties sig_m.
#'              Calculate the measurement matrix M, which  has dimensions N x G
#'              and for which each element is M_ig = p(phi_m,i|tau_g). The total
#'              uncertainty of the measurement comes from measurement error
#'              (SIG_M, calculated using the measurement error for each
#'              measurement) and the calibration curve error (SIG_c, calculated
#'              using the uncertainty for the calibration curve at each grid
#'              point). These uncertainties (and the associated measurements)
#'              should already be "projected" to 1950 equivalents. By default
#'              the measurement matrix, M, is multiplied by dtau so that M * f
#'              is an approximation to the integral over p(t|th) * p(phi|t)
#'              using a Riemann sum with the density calculated at the points
#'              tau and the width for each point being dtau. Alternatively,
#'              the trapezoidal rule can be used (for details see
#'              bd_calc_trapez_weights). If the spacing of tau is irregular,
#'              the trapezoidal rule must be used. An error is thrown if tau is
#'              irregularly spaced and the the trapezoidal rule is not used
#'
#' @param tau A vector of calendar dates indexed by g
#' @param phi_m A vector of fraction moderns indexed by i
#' @param sig_m A vector of standard deviations for phi_m indexed by i
#' @param calibDf Calibration curve (see bd_load_calib_curve)
#' @param addCalibUnc (default TRUE) Whether to add calibration uncertainty
#' @param useTrapez (default FALSE) Whether to use the trapezoidal rule for integration
#'
#' @export

bd_calc_meas_matrix <- function(tau, phi_m, sig_m, calibDf, addCalibUnc = T, useTrapez = F) {
  # First, check the consistency of the spacing in tau and the value of
  # useTrapez, which must be TRUE if tau is irregularly spaced.
  irreg <- length(unique(diff(tau))) != 1
  if (irreg && !useTrapez) {
    stop("tau is irregularly spaced but useTrapez is FALSE")
  }

  # tau is in AD
  tau_BP <- 1950 - tau

  # extract the calibration curve variables and convert to fraction modern
  tau_curve <- rev(calibDf$yearBP)
  mu_c_curve <- exp(-rev(calibDf$uncalYearBP) / 8033)
  sig_c_curve <- rev(calibDf$uncalYearBPError) * mu_c_curve / 8033

  # Interpolate curves at tau_BP to yield mu_c and sig_c
  mu_c <- stats::approx(tau_curve, mu_c_curve, tau_BP)
  mu_c <- mu_c$y
  sig_c <- stats::approx(tau_curve, sig_c_curve, tau_BP)
  sig_c <- sig_c$y

  PHI_m <- replicate(length(tau_BP), phi_m)
  SIG_m <- replicate(length(tau_BP), sig_m)

  MU_c <- t(replicate(length(phi_m), mu_c))
  if (addCalibUnc) {
    SIG_c <- t(replicate(length(sig_m), sig_c))
    SIG_sq <- SIG_m^2 + SIG_c^2
  } else {
    SIG_sq <- SIG_m^2
  }

  M <- exp(-(PHI_m - MU_c)^2 / (SIG_sq) / 2) / sqrt(SIG_sq) / sqrt(2 * pi)

  # If necessary, add the integration widths
  if (!useTrapez) {
    M <- M * (tau[2] - tau[1])
  } else {
    G <- length(tau)
    dtauVect <- bd_calc_trapez_weights(tau)
    dtauMat <- t(replicate(length(phi_m), dtauVect))
    M <- M * dtauMat
  }

  return(M)
}
