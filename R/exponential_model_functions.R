#' @title Sample from a truncated exponential distribution
#'
#' @description Make N samples from a truncated exponential distribution with the density function r0*exp(r0*t). The rate r0 can be positive (growth) or negative (decay). The density function is normalized to integrate to 1 on the interval taumin to taumax.
#'
#' @param N The number of samples to make
#' @param r0 The rate parameter
#' @param taumin The lower limit for truncation
#' @param taumax The upper limit for truncation
#'
#' @return A vector of N samples between taumin and taumax
#'
#' @author Michael Holton Price <MichaelHoltonPrice@@gmail.com>
#'
#' @export

sample_trunc_exp <- function(N, r0, taumin, taumax) {
  # Call RGeode::rexptr for the sampling. rexptr cannot be used direclty
  # because it does not support r0 being zero or negative. The input to rexptr
  # is -r0.
  if (r0 > 0) {
    return(taumax - RGeode::rexptr(N, r0, range = c(0, taumax - taumin)))
  } else if (r0 == 0) {
    return(stats::runif(N, taumin, taumax))
  } else {
    return(RGeode::rexptr(N, -r0, range = c(taumin, taumax)))
  }
}
