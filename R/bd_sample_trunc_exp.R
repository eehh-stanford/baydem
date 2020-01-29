#' @title Sample from a truncated exponential distribution
#'
#' @description Make N samples from a truncated exponential distribution with the density function r0*exp(r0*y). The rate r0 can be positive (growth) or negative (decay). The density function is normalized to integrate to 1 on the interval ymin to ymax.
#'
#' @param N The number of samples to make
#' @param r0 The rate parameter
#' @param ymin The lower limit for truncation
#' @param ymax The upper limit for truncation
#'
#' @return A vector of N samples between ymin and ymax
#'
#' @author Michael Holton Price <MichaelHoltonPrice@@gmail.com>
#'
#' @export

bd_sample_trunc_exp <- function(N, r0, ymin, ymax) {
  # Call RGeode::rexptr for the sampling. rexptr cannot be used direclty
  # because it does not support r0 being zero or negative. The input to rexptr
  # is -r0.
  if (r0 > 0) {
    return(ymax - RGeode::rexptr(N, r0, range = c(0, ymax - ymin)))
  } else if (r0 == 0) {
    return(runif(N, ymin, ymax))
  } else {
    return(RGeode::rexptr(N, -r0, range = c(ymin, ymax)))
  }
}
