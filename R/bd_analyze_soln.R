#' @title Analyze the result of a call to \code{bd_do_inference}
#'
#' @description \code{bd_do_inference} calls Stan to do Bayesian inference by generating
#' a sample of parameters from the posterior of theta (or th).
#' This function analyzes the result of that inference.
#' In particular, it calculates the quantiles of the density function and growth rate.
#'
#' @details `soln` is the result of a call to \code{bd_do_inference}.
#' It contains both the resulting samples and the parameters used in the inference,
#' such as the hyperparameters (see \code{bd_do_inference} for further details).
#' The primary thing \code{bd_analyze_soln} does is calculate quantiles of both
#' the parameterized density and growth rate. For example,
#' for a calendar date y each sample yields a density and growth rate.
#' The quantile is the value of the density or growth rate such that
#' a given proportion of samples are smaller than that value.
#' The probabilities used to calculate these quantiles are `probs = c(lev, 0.5, 1-lev)`,
#' where `lev` is the level (0.025 by default, so that 95% of the observations
#' lie between the first and last quantile bands).
#'
#' In addition, \code{bd_analyze_soln} identifies calendar dates for which
#' the growth rate quantiles defined by `lev` and `1 - lev` do not contain zero.
#' This indicates significant positive or negative growth for the density curve.
#' The output vector `growthState` codes calendar dates by growth state as 'negative',
#' 'zero', and 'positive'. For the Gaussian mixture parameterization of the density,
#' the rate is not typically meaningful near the calendar date boundaries where it
#' increases linearly as the calendar date goes to positive or negative infinity.
#' The parameter `rateProp` provides control on how calendar dates are classified by
#' growth rate near these boundaries. In particular, the calendar dates with a cumulative
#' density (50% quantile) below `rateProp` (for the lower boundary) or above `1 - rateProp`
#' (for the upper boundary) are classified as 'missing' in `growthState`.
#' By default, `rateProp` is NA and no calendar dates are classified as missing.
#'
#' By default, a summary is done for each sample by calling
#' bd_summarize_sample. This is not done of doSummary is FALSE
#'
#' @param soln The solution, a list-like object of class bd_soln (see \code{bd_do_inference})
#' @param tau (optional) The calendar dates at which to evaluate densities.
#' If tau is not input, tau is built from the hyperparameters.
#' @param th_sim (optional) The known parameters used to create simulation data
#' @param lev (default: 0.025) The level to use for the quantile bands
#' @param rateProp (optional) The cumulative density needed to define rate growth bands
#' @param doSummary [Default TRUE] Whether to summarize each sample by calling bd_summarize_sample
#'
#' @return A list with information on the quantiles of the density function and growth rate (and sample summaries)
#'
#' @export
bd_analyze_soln <- function(soln, tau = NA, th_sim = NA, lev = 0.025, rateProp = NA, doSummary = T) {
  if (all(is.na(tau))) {
    tau <- seq(soln$prob$hp$taumin, soln$prob$hp$taumax, by = soln$prob$hp$dtau)
  }

  probs <- c(lev, 0.5, 1 - lev) # The probabilities to use for quantiles

  # Determine tau spacing, dtau, and ensure that tau is evenly spaced
  dtau <- unique(diff(tau))
  if (length(dtau) > 1) {
    stop("tau should by uniformily spaced")
  }

  # Extract the samples of theta in the variable TH. TH is matrix like object,
  # possibly of a specific class (e.g., gaussmix) with dimensions
  # numSamp x numParam, where numSamp is the number of samples and numParam is
  # the number of parameters (length of th).
  TH <- bd_extract_param(soln$fit)

  numMix <- ncol(TH) / 3 # This assumes a gassian mixture fit. Future updates may generalize this
  numSamp <- nrow(TH)
  numGrid <- length(tau)

  # Calculate the pdf matrix, which is the density of the parametric model for
  # theta for each sample and each grid point. fMat has dimensions
  # N x G, where N is the number of samples in TH and G is the length of the vector tau.
  # Because bd_calc_gauss_mix_pdf_mat is called with taumin and taumax, the density
  # is normalized to integrate to 1 on the interval taumin to taumax.
  fMat <- bd_calc_gauss_mix_pdf_mat(TH, tau, ymin = soln$prob$hp$taumin, ymax = soln$prob$hp$taumax)

  # Calculate the rate for each sample and grid point (f' / f, where f is density)
  rateMat <- bd_calc_gauss_mix_pdf_mat(TH, tau, ymin = soln$prob$hp$taumin, ymax = soln$prob$hp$taumax, type = "rate")

  # Calculate the quantiles of the normalized density matrix
  Qdens <- bd_calc_quantiles(fMat, probs)

  # Normalized 50% densities (not normalized to integrate to 1)
  f50 <- Qdens[2, ] # The second row gives the 50% quantiles

  # Restrict to indices with enough probability mass (if necessary)
  if (!is.na(rateProp)) {
    rateInd <- which(cumsum(f50 * dtau) > rateProp & rev(cumsum(rev(f50) * dtau)) > rateProp)
  } else {
    rateInd <- 1:length(f50)
  }

  # Identify regions with growth rates that differ from zero per the input quantile level (lev)
  # growthState0 is -1 for significant negative growth, 1 for significant positive growth, and 0 otherwise
  Qrate <- bd_calc_quantiles(rateMat[, rateInd], probs)
  growthState0 <- rep("zero", length(rateInd)) # growthState0 indices in rateInd
  growthState0[Qrate[2, ] > 0 & Qrate[1, ] > 0] <- "positive"
  growthState0[Qrate[2, ] < 0 & Qrate[3, ] < 0] <- "negative"
  growthState <- rep("missing", length(tau))
  growthState[rateInd] <- growthState0 # growthState for all indices


  # Calculate the measurement matrix
  M <- bd_calc_meas_matrix(tau, soln$prob$phi_m, soln$prob$sig_m, soln$prob$calibDf)

  # Calculate and normalize the summed probability density vector
  f_spdf <- colSums(M)
  f_spdf <- f_spdf / sum(f_spdf) / dtau

  out <- list(
    tau = tau,
    f_spdf = f_spdf,
    Qdens = Qdens,
    Qrate = Qrate,
    probs = probs,
    rateProp = rateProp,
    rateInd = rateInd,
    growthState = growthState,
    dtau = dtau
  )
  class(out) <- "bd_analysis"

  if (doSummary) {
    summList <- list()
    for (n in 1:numSamp) {
      th <- TH[n, ]
      summList[[n]] <- bd_summarize_trunc_gauss_mix_sample(th, soln$prob$hp$taumin, soln$prob$hp$taumax)
    }
    out$summList <- summList
  }

  haveSim <- !all(is.na(th_sim))
  if (haveSim) {
    f_sim <- bd_calc_gauss_mix_pdf(th_sim, tau, ymin = soln$prob$hp$taumin, ymax = soln$prob$hp$taumax)
    rate_sim <- bd_calc_gauss_mix_pdf(th_sim, tau, ymin = soln$prob$hp$taumin, ymax = soln$prob$hp$taumax, type = "rate")
    out$f_sim <- f_sim
    out$rate_sim <- rate_sim
  }
  return(out)
}
