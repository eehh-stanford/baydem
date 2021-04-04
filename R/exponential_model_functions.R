#' @title
#' Sample from a truncated exponential distribution
#'
#' @description
#' Make N samples from a truncated exponential distribution with the density
#' function r0*exp(r0*t). The rate r0 can be positive (growth) or negative
#' (decay). The density function is normalized to integrate to 1 on the interval
#' tau_min to tau_max.
#'
#' @param N The number of samples to make
#' @param r0 The rate parameter
#' @param tau_min The lower limit for truncation
#' @param tau_max The upper limit for truncation
#'
#' @return A vector of N samples between tau_min and tau_max
#'
#' @export
sample_trunc_exp <- function(N, r0, tau_min, tau_max) {
  # Call RGeode::rexptr for the sampling. rexptr cannot be used direclty
  # because it does not support r0 being zero or negative. The input to rexptr
  # is -r0.
  if (r0 > 0) {
    return(tau_max - RGeode::rexptr(N, r0, range = c(0, tau_max - tau_min)))
  } else if (r0 == 0) {
    return(stats::runif(N, tau_min, tau_max))
  } else {
    return(RGeode::rexptr(N, -r0, range = c(tau_min, tau_max)))
  }
}
