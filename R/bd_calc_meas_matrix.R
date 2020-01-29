#' @title Calculate the measurement matrix given N observations and G grid points
#'
#' @description ygrid is a vector of calendar dates indexed from g = 1,2,...,G.
#'              phi_m is a vector of n = 1,2,...,N radiocarbon determinations
#'              (fraction modern) with associated uncertainties sig_m.
#'              Calculate the measurement matrix M, which  has dimensions N x G
#'              and for which each element is M_ng = p(phi_m,n|y_g). The total
#'              uncertainty of the measurement comes from measurement error
#'              (SIG_M, calculated using the measurement error for each
#'              measurement) and the calibration curve error (SIK_k, calculated
#'              using the uncertainty for the calibration curve at each grid
#'              point). These uncertainties (and the associated measurements)
#'              should already be "projected" to 1950 equivalents (e.g., in
#'              bd_convert_date_samp_to_c14_samp). By default, the grid spacing
#'              of ygrid is accounted for in each element of the measurement
#'              matrix, M_ng, so that the matrix multiplication M * f is an
#'              approximation to the integral over p(y|th) * p(phi|y); the
#'              trapezoidal rule is assumed for this, so that the spacing to
#'              use for the the integration width is
#'              dy_g = (y_(g+1) - y_(g-1)) / 2, where the conventions y_0 = y_1
#'              and y_(G+1) = y_G are used: M_ng -> M_ng * dy_g
#'
#' @param ygrid A vector of calendar dates indexed by g
#' @param phi_m A vector of fraction moderns indexed by n
#' @param sig_m A vector of standard deviations for phi_m by n
#' @param calibDf Calibration curve (see bd_load_calib_curve)
#' @param useSpacing (default TRUE) Whether to account for the spacing of ygrid (see description)
#' @param addCalibUnc (default TRUE) Whether to add calibration uncertainty
#'
#' @export

bd_calc_meas_matrix <- function(ygrid, phi_m, sig_m, calibDf, useSpacing = T, addCalibUnc = T) {
  # ygrid is in AD
  ygrid_BP <- 1950 - ygrid

  # extract the calibration curve variables and convert to fraction modern
  y_curve <- rev(calibDf$yearBP)
  mu_k_curve <- exp(-rev(calibDf$uncalYearBP) / 8033)
  sig_k_curve <- rev(calibDf$uncalYearBPError) * mu_k_curve / 8033

  # Interpolate curves at ygrid_BP to yield mu_k and sig_k
  mu_k <- stats::approx(y_curve, mu_k_curve, ygrid_BP)
  mu_k <- mu_k$y
  sig_k <- stats::approx(y_curve, sig_k_curve, ygrid_BP)
  sig_k <- sig_k$y

  PHI_m <- replicate(length(ygrid_BP), phi_m)
  SIG_m <- replicate(length(ygrid_BP), sig_m)

  MU_k <- t(replicate(length(phi_m), mu_k))
  if (addCalibUnc) {
    SIG_k <- t(replicate(length(sig_m), sig_k))
    SIG_sq <- SIG_m^2 + SIG_k^2
  } else {
    SIG_sq <- SIG_m^2
  }

  M <- exp(-(PHI_m - MU_k)^2 / (SIG_sq) / 2) / sqrt(SIG_sq) / sqrt(2 * pi)

  # If necessary, add the integration widths
  if (useSpacing) {
    G <- length(ygrid)
    dyVect <- rep(NA, length(ygrid))
    indCent <- 2:(G - 1)
    dyVect[indCent] <- (ygrid[indCent + 1] - ygrid[indCent - 1]) / 2
    dyVect[1] <- (ygrid[2] - ygrid[1]) / 2
    dyVect[G] <- (ygrid[G] - ygrid[G - 1]) / 2
    dyMat <- t(replicate(length(phi_m), dyVect))
    M <- M * dyMat
  }

  return(M)
}
