#' @title Identify growth periods and the peak value for a truncated Gaussian mixture
#'
#' @description
#' The input vector th parameterizes a Gaussian mixture, and taumin / taumax
#' give the limits of truncation. Summarize the sample by identifying growth /
#' decay periods and the peak value using the following procedure.
#'
#' (1) Calculate the derivative, f'(t), at the points t = seq(taumin,taumax,len=N),
#'     where N is 1000 by default.
#'
#' (2) Identify points where f'(t) changes sign, then numerically estimate the
#'     crossing point between the two t values where there was a sign change.
#'
#' (3) Create a vector of critical points, tcrit, which includes taumin / taumax
#'     as well as the crossing points found in the preceding step.
#'
#' (4) Calculate the density at the critical points to identify the peak value,
#'     fpeak, and corresponding calendar date, tpeak, as well as the index of
#'     of the peak in tcrit, indPeak.
#'
#' (5) For each time period (the length(tpeak)-1 durations defined by ypeak)
#'     determine the sign of the density function, f(t), and create a character
#'     vector, slope, that has the value 'pos' if f(t) is positive and 'neg' if
#'     f(t) is negative.
#'
#' (6) Finally, create a character vector, pattern, that appends the index of
#'     the peak in tcrit (converted to a character) to the character vector
#'     slope. This defines a unique pattern of the sample that takes into
#'     account periods of growth / decline and the relative location of the
#'     peak.
#'
#' @param  th The Gaussian mixture parameterization
#' @param  taumin The lower truncation value
#' @param  taumin The upper truncation value
#' @param  N (Default 1000) The number of points use for identifying slope changes
#'
#' @return A list consisting of tlo / thi (specifying the time periods), indPeak, tpeak, fpeak, and pattern (see Description)
#'
#' @export
bd_summarize_trunc_gauss_mix_sample <- function(th, taumin, taumax, N = 1000) {
  # (1) Calculate the derivative of the density
  K <- length(th) / 3 # Number of mixtures
  t <- seq(taumin, taumax, len = N)
  fprime <- bd_calc_gauss_mix_pdf(th, t, taumin, taumax, type = "derivative")

  # (2) Identify locations in t where the derivative changes sign. This happens
  #     if fprime[n] * fprime[n+1] is less than zero. Then, numerically
  #     estimate the exact t-value of the crossing.
  ind <- which(fprime[1:(length(fprime) - 1)] * fprime[2:length(fprime)] < 0)
  M <- length(ind) # Number of cross-overs

  # Vectors for t / f values of crossings
  tcross <- rep(NA, M)
  fcross <- rep(NA, M)

  if (M > 0) {
    # Objective function to maximize
    rootFun <- function(t) {
      return(bd_calc_gauss_mix_pdf(th, t, taumin, taumax, type = "derivative"))
    }

    # Iterate over crossings
    for (m in 1:M) {
      root <- uniroot(rootFun, lower = t[ind[m]], upper = t[ind[m] + 1])
      tcross[m] <- root$root
      fcross[m] <- bd_calc_gauss_mix_pdf(th, tcross[m], taumin, taumax)
    }
  }

  # (3-4) Create the vector of critical points, calculate densities, and
  #       identify peak
  tcrit <- c(taumin, tcross, taumax)
  fcrit <- c(bd_calc_gauss_mix_pdf(th, taumin, taumin, taumax), fcross, bd_calc_gauss_mix_pdf(th, taumin, taumin, taumax))
  indPeak <- which.max(fcrit)
  tpeak <- tcrit[indPeak]
  fpeak <- fcrit[indPeak]

  # (5) Create tlo, thi, and slope
  numPer <- length(tcrit) - 1 # Number of periods
  tlo <- tcrit[1:numPer]
  thi <- tcrit[2:(numPer + 1)]
  df <- diff(fcrit)
  slope <- rep("pos", numPer)
  slope[df < 0] <- "neg"

  # (6) Create the pattern (then return the result)
  pattern <- c(slope, as.character(indPeak))

  return(list(periods = data.frame(tlo = tlo, thi = thi, slope = slope), indPeak = indPeak, tpeak = tpeak, fpeak = fpeak, pattern = pattern))
}
