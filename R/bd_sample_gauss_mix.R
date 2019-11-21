# @keywords
# Description
#   Sample from a Gaussian mixture. This provides calendar dates of radiocarbon
#   samples from the demographic model specified by Gaussian mixture.
#
# Example calls(s)
#
#   samp <- bd_sample_gauss_mix(N,th)
#
# Input(s)
#   Name    Type           Description
#   N       integer        The number of samples
#   th      vector-like    A vector-like object of length 3*K with the
#                          following entries:
#                          pik   -- [K entries] Weight of the k-th mixture
#                          muk   -- [K entries] Mean of the k-th mixture
#                          sigk  -- [K entries] Standard deviation of the k-th
#                                               mixture
#
# Output(s)
#   Name    Type           Description
#   samp    vector         The samples (length = N)

bd_sample_gauss_mix <- function(N, th, ymin = NA, ymax = NA) {
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
                                               mixCoeff = th[1:2])
  } else if (K == 3) {
    normMix <- distr::UnivarMixingDistribution(distr::Norm(mean = th[4], sd = th[7]),
                                               distr::Norm(mean = th[5], sd = th[8]), ,
                                               distr::Norm(mean = th[6], sd = th[9]),
                                               mixCoeff = th[1:3])
  } else if (K == 4) {
    normMix <- distr::UnivarMixingDistribution(distr::Norm(mean = th[5], sd = th[9]),
                                               distr::Norm(mean = th[6], sd = th[10]), ,
                                               distr::Norm(mean = th[7], sd = th[11]),
                                               distr::Norm(mean = th[8], sd = th[12]),
                                               mixCoeff = th[1:4])
  } else {
    stop("The maximum number of supported mixture components is 4")
  }

  # Use distr to sample from a two-component, truncated Gaussian mixture
  if (!is.na(ymin) && !is.na(ymax)) {
    normMixTrunc <- distr::Truncate(normMix, ymin, ymax)
    samp <- distr::r(normMixTrunc)(N)
  } else {
    samp <- distr::r(normMix)(N)
  }
  return(samp)
}
