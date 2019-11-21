# @keywords
# Description
# Sample from the prior, given the type of fit and associated hyperparameters
# for that type of fit
#
# Example calls(s)
#
#   samp <- bd_sample_prior(hp,N)
#
# Input(s)
#   Name    Type           Description
#   hp      list           The hyperparameters
#   N       integer        Number of samples
#
# Output(s)
#   Name    Type           Description
#   samps   list           N samples from the prior for the given fit type.
#                          samps is a list of lists, which is the format
#                          expected by stan for initializing chains. Details for
#                          each fit type:
#
#                gaussmix  sig [K x 1]  -- vector of standard deviations
#                          mu  [K x 1]  -- vector of means
#                          pi  [K x 1]  -- vector of mixture proportions
#' @export

bd_sample_prior <- function(hp, N) {
  samps <- list()
  if (hp$fitType == "gaussmix") {
    for (n in 1:N) {
      samps[[n]] <- list(fitType = hp$fitType,
                         sig = abs(stats::rgamma(hp$K,
                                                 shape = hp$sigAlpha,
                                                 rate = hp$sigBeta)),
                         mu = sort(stats::runif(hp$K, hp$ymin, hp$ymax)),
                         pi = as.vector(gtools::rdirichlet(1, rep(hp$dirichParam, hp$K))))
    }
  } else {
    stop(paste("Unrecognized fit type:", hp$fitType))
  }
  return(samps)
}
