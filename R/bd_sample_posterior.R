# Description
# Get a random sample of the posteriors from the stan fit.
#
# Example calls(s)
#
#   samp <- bd_sample_posterior(fit,N)
#
# Input(s)
#   Name    Type           Description
#   fit     rstan          The rstan fit object
#   hp      list           The hyperparameters
#   N       integer        Number of samples
#
# Output(s)
#   Name    Type           Description
#   samps   list           N samples from the posterior for the given fit type.
#                          samps is a list of lists, which is the format
#                          expected by stan for initializing chains. Details for
#                          each fit type:
#
#                gaussmix  sig [K x 1]  -- vector of standard deviations
#                          mu  [K x 1]  -- vector of means
#                          pi  [K x 1]  -- vector of mixture proportions
#                          ymin         -- minimum calendar date
#                          ymax         -- maximum calendar date

bd_sample_posterior <- function(fit, N) {
  samps <- bd_extract_param(fit)
  numSamp <- length(samps)
  ind <- sample.int(numSamp, size = N)
  return(samps[ind])
}
