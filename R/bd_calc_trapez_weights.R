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
#' @param tau A vector of locations where the function is sampled, possibly irregularly
#'
#' @return A vector of integration weights the same length as tau

#' @export
bd_calc_trapez_weights <- function(tau) {
  G <- length(tau)
  weightVect <- rep(NA, length(tau))
  indCent <- 2:(G - 1)
  weightVect[indCent] <- (tau[indCent + 1] - tau[indCent - 1]) / 2
  weightVect[1] <- (tau[2] - tau[1]) / 2
  weightVect[G] <- (tau[G] - tau[G - 1]) / 2
  return(weightVect)
}
