#' @title Calculation the calibration curve fraction taumodern for input
#' calendar dates
#'
#' @details
#' The input calibration data frame has three columns: year_BP, uncal_year_BP,
#' and uncal_year_BP_error. For each input calendar date, tau, use the
#' calibration curve information to estimate the fraction modern of the
#' calibration curve. If tau is not specified, calib_df$uncal_year_BP is used as
#' the dates for the calculation. By default, dates are assumed to be AD, but
#' this can be changed using the optional input is_BP, which is FALSE by
#' default.
#'
#' @param calib_df The calibration data frame, with columns year_BP,
#'   uncal_year_BP, and uncal_year_BP_error
#' @param tau The calendar dates tau (if not input, calib_df$uncal_year_BP is
#'   used)
#' @param is_BP (default FALSE) Whether the input dates are before present (BP),
#'   as opposed to AD
#'
#' @return The vector of calibration curve fraction modern values
#' @export
calc_calib_curve_frac_modern <-
  function(calib_df,
           tau = NA,
           is_BP = FALSE) {
    phi_curve <- exp(-calib_df$uncal_year_BP / 8033)
    if (all(is.na(tau))) {
      return(phi_curve)
    }

    tau_curve <- 1950 - calib_df$year_BP
    if (is_BP) {
      tau <- 1950 - tau # Convert to AD
    }

    phi <- stats::approx(tau_curve, phi_curve, tau)
    phi <- phi$y
    return(phi)
  }


#' @title Calculate the measurement matrix given N observations and G grid
#' points
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
#'              calc_trapez_weights). If the spacing of tau is irregular,
#'              the trapezoidal rule must be used. An error is thrown if tau is
#'              irregularly spaced and the the trapezoidal rule is not used
#'
#' @param tau A vector of calendar dates indexed by g
#' @param phi_m A vector of fraction moderns indexed by i
#' @param sig_m A vector of standard deviations for phi_m indexed by i
#' @param calib_df Calibration curve (see load_calib_curve)
#' @param add_calib_unc (default TRUE) Whether to add calibration uncertainty
#' @param use_trapez (default FALSE) Whether to use the trapezoidal rule for
#'   integration
#'
#' @export

calc_meas_matrix <- function(tau,
                             phi_m,
                             sig_m,
                             calib_df,
                             add_calib_unc = T,
                             use_trapez = F) {
  # First, check the consistency of the spacing in tau and the value of
  # use_trapez, which must be TRUE if tau is irregularly spaced.
  irreg <- length(unique(diff(tau))) != 1
  if (irreg && !use_trapez) {
    stop("tau is irregularly spaced but useTrapez is FALSE")
  }

  # tau is in AD
  tau_BP <- 1950 - tau

  # extract the calibration curve variables and convert to fraction modern
  tau_curve <- rev(calib_df$year_BP)
  mu_c_curve <- exp(-rev(calib_df$uncal_year_BP) / 8033)
  sig_c_curve <- rev(calib_df$uncal_year_BP_error) * mu_c_curve / 8033

  # Interpolate curves at tau_BP to yield mu_c and sig_c
  mu_c <- stats::approx(tau_curve, mu_c_curve, tau_BP)
  mu_c <- mu_c$y
  sig_c <- stats::approx(tau_curve, sig_c_curve, tau_BP)
  sig_c <- sig_c$y

  PHI_m <- replicate(length(tau_BP), phi_m)
  SIG_m <- replicate(length(tau_BP), sig_m)

  MU_c <- t(replicate(length(phi_m), mu_c))
  if (add_calib_unc) {
    SIG_c <- t(replicate(length(sig_m), sig_c))
    SIG_sq <- SIG_m^2 + SIG_c^2
  } else {
    SIG_sq <- SIG_m^2
  }

  M <- exp(-(PHI_m - MU_c)^2 / (SIG_sq) / 2) / sqrt(SIG_sq) / sqrt(2 * pi)

  # If necessary, add the integration widths
  if (!use_trapez) {
    M <- M * (tau[2] - tau[1])
  } else {
    G <- length(tau)
    dtau_vect <- calc_trapez_weights(tau)
    dtau_mat <- t(replicate(length(phi_m), dtau_vect))
    M <- M * dtau_mat
  }

  return(M)
}

#' @title Calculate integration weights assuming the midpoint rule
#'
#' @description tau is a vector of locations where a function to be integrated
#'              is evaluated (the function values, f, are not input). Let
#'              tau_g be the locations of integration and f_g the corresponding
#'              function values for g=1, 2, ... G. Assuming trapezoidal
#'              integration at the midpoints between elements of tau, the
#'              weights to use for integration are
#'              dtau_g = (tau_(g+1) - tau_(g-1)) / 2, where the conventions
#'              tau_0 = tau_1 and tau_(G+1) = tau_G are used.
#'
#' @param tau A vector of locations where the function is sampled, possibly
#'   irregularly
#'
#' @return A vector of integration weights the same length as tau

#' @export
calc_trapez_weights <- function(tau) {
  G <- length(tau)
  weight_vect <- rep(NA, length(tau))
  ind_cent <- 2:(G - 1)
  weight_vect[ind_cent] <- (tau[ind_cent + 1] - tau[ind_cent - 1]) / 2
  weight_vect[1] <- (tau[2] - tau[1]) / 2
  weight_vect[G] <- (tau[G] - tau[G - 1]) / 2
  return(weight_vect)
}

#if (getRversion() >= "2.15.1") {
#  utils::globalVariables(c(
#    "CAL BP",
#    "14C age",
#    "Error"
#  ))
#}

#' @title Load Calibration Curve
#'
#' @description Parse and return the radiocarbon calibration curve stored in
#'   data
#'
#' @param calib_curve Name of calibration curve
#'
#' @return The calibration dataframe, with columns year_BP, uncal_year_BP, and
#'   uncal_year_BP_error
#'
#' @importFrom magrittr %>%
#' @export

load_calib_curve <- function(calib_curve) {
  if (!(calib_curve %in%
    c("intcal20", "marine20", "shcal20", "intcal13", "marine13", "shcal13"))) {
    stop(paste("Unknown calibration curve name:", calib_curve))
  }

  if (calib_curve == "intcal20") {
    calib_curve <- baydem::intcal20
  } else if (calib_curve == "marine20") {
    calib_curve <- baydem::marine20
  } else if (calib_curve == "shcal20") {
    calib_curve <- baydem::shcal20
  } else if (calib_curve == "intcal13") {
    calib_curve <- baydem::intcal13
  } else if (calib_curve == "marine13") {
    calib_curve <- baydem::marine13
  } else if (calib_curve == "shcal13") {
    calib_curve <- baydem::shcal13
  }

  calib_df <- calib_curve %>%
    dplyr::select(
      `CAL BP`,
      `14C age`,
      Error
    ) %>%
    dplyr::rename(
      year_BP = `CAL BP`,
      uncal_year_BP = `14C age`,
      uncal_year_BP_error = `Error`
    )

  return(calib_df)
}

#' @title Find calendar date from fraction modern value given known bounding
#' indices
#'
#' @details
#' tau_curve and phi_curve give the calendar date and fraction modern of the
#' radiocarbon calibration curve. It is known that the calendar date lies
#' lies between `tau_curve[ii_lo]` and `tau_curve[ii_hi]`. Interpolate to find
#' the calendar date corresponding to the input phi_known. This function is
#' called by assess_calib_curve_equif.
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
phi2tau <- function(tau_curve, phi_curve, phi_known, ii_lo, ii_hi) {
  tau <- tau_curve[ii_lo] +
    (tau_curve[ii_hi] - tau_curve[ii_lo]) *
      (phi_known - phi_curve[ii_lo]) / (phi_curve[ii_hi] - phi_curve[ii_lo])
  return(tau)
}
