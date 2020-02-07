#' @title Sample from a possibly truncated Gaussian mixture
#'
#' @details `N` is the number of samples to draw and `th` specifies the Gaussian
#' mixture with the ordering [pi_1,...,pi_K,mu_1,...,mu_K,sig_1,...,sig_K],
#' where `K` is the number of mixture components. Optionally, the samples are
#' drawn on the truncated interval `taumin` to `taumax`. Because of limitations
#' in the package `distr`, the maximum number of mixture components suppored is
#' K=4.
#'
#' @param N Number of samples
#' @param th Parameterization vector for the Gaussian mixture
#' @param taumin (Optional) Lower bound for samples
#' @param taumax (Optional) Upper bound for samples
#'
#' @return N samples from the Gaussian mixture
#'
#' @export
bd_sample_gauss_mix <- function(N, th, taumin = NA, taumax = NA) {
  K <- length(th) / 3 # Number of  mixtures

  # (Somewhat awkwardly), define the mixing distribution directly for K = 1
  # to 4 directly. This is necessary because truncate expects an
  # AbscontDistribution but creating a mixture distribution via vectors
  # creates a UnivarMixingDistribution, which cannot be cast to an
  # AbscontDistribution. Perhaps this shortcoming of distr will be addressed
  # in the future. In the meantime, the number of mixtures is limited to 4
  if (K == 1) { # Allow K = 1 (i.e., a Gaussian, not a mixture) as a special case
    normMix <- distr::Norm(mean = th[2], sd = th[3])
  } else if (K == 2) {
    normMix <- distr::UnivarMixingDistribution(distr::Norm(mean = th[3], sd = th[5]),
      distr::Norm(mean = th[4], sd = th[6]),
      mixCoeff = th[1:2]
    )
  } else if (K == 3) {
    normMix <- distr::UnivarMixingDistribution(distr::Norm(mean = th[4], sd = th[7]),
      distr::Norm(mean = th[5], sd = th[8]), ,
      distr::Norm(mean = th[6], sd = th[9]),
      mixCoeff = th[1:3]
    )
  } else if (K == 4) {
    normMix <- distr::UnivarMixingDistribution(distr::Norm(mean = th[5], sd = th[9]),
      distr::Norm(mean = th[6], sd = th[10]), ,
      distr::Norm(mean = th[7], sd = th[11]),
      distr::Norm(mean = th[8], sd = th[12]),
      mixCoeff = th[1:4]
    )
  } else {
    stop("The maximum number of supported mixture components is 4")
  }

  # Use distr to sample from a two-component, truncated Gaussian mixture
  if (!is.na(taumin) && !is.na(taumax)) {
    normMixTrunc <- distr::Truncate(normMix, taumin, taumax)
    samp <- distr::r(normMixTrunc)(N)
  } else {
    samp <- distr::r(normMix)(N)
  }
  return(samp)
}
