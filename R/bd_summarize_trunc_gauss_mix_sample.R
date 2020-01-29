#' @title Identify growth periods and the peak value for a truncated Gaussian mixture
#'
#' @description
#' The input vector th parameterizes a Gaussian mixture, and ymin / ymax give
#' the limits of truncation. Summarize the sample by identifying growth / decay
#' periods and the peak value using the following procedure.
#'
#' (1) Calculate the derivative, f'(y), at the points y = seq(ymin,ymax,len=N),
#'     where N is 1000 by default.
#'
#' (2) Identify points where f'(y) changes sign, then numerically estimate the
#'     crossing point between the two y values where there was a sign change.
#'
#' (3) Create a vector of critical points, ycrit, which includes ymin / ymax
#'     as well as the crossing points found in the preceding step.
#'
#' (4) Calculate the density at the critical points to identify the peak value,
#'     fpeak, and corresponding calendar date, ypeak, as well as the index of
#'     of the peak in ycrit, indPeak.
#'
#' (5) For each time period (the length(ypeak)-1 durations defined by ypeak)
#'     determine the sign of the density function, f(y), and create a character
#'     vector, slope, that has the value 'pos' if f(y) is positive and 'neg' if
#'     f(y) is negative.
#'
#' (6) Finally, create a character vector, pattern, that appends the index of
#'     the peak in ycrit (converted to a character) to the character vector
#'     slope. This defines a unique pattern of the sample that takes into
#'     account periods of growth / decline and the relative location of the
#'     peak.
#'
#' @param  th The Gaussian mixture parameterization
#' @param  ymin The lower truncation value
#' @param  ymin The upper truncation value
#' @param  N (Default 1000) The number of points use for identifying slope changes
#'
#' @return A list consisting of ylo / yhi (specifying the time periods), indPeak, ypeak, fpeak, and pattern (see Description)
#'
#' @export
bd_summarize_trunc_gauss_mix_sample <- function(th, ymin, ymax, N = 1000) {
  # (1) Calculate the derivative of the density
  K <- length(th) / 3 # Number of mixtures
  y <- seq(ymin, ymax, len = N)
  fprime <- bd_calc_gauss_mix_pdf(th, y, ymin, ymax, type = "derivative")

  # (2) Identify locations in y where the derivative changes sign. This happens
  #     if fprime[n] * fprime[n+1] is less than zero. Then, numerically
  #     estimate the exact y-value of the crossing.
  ind <- which(fprime[1:(length(fprime) - 1)] * fprime[2:length(fprime)] < 0)
  M <- length(ind) # Number of cross-overs

  # Vectors for y / f values of crossings
  ycross <- rep(NA, M)
  fcross <- rep(NA, M)

  if (M > 0) {
    # Objective function to maximize
    rootFun <- function(y) {
      return(bd_calc_gauss_mix_pdf(th, y, ymin, ymax, type = "derivative"))
    }

    # Iterate over crossings
    for (m in 1:M) {
      root <- uniroot(rootFun, lower = y[ind[m]], upper = y[ind[m] + 1])
      ycross[m] <- root$root
      fcross[m] <- bd_calc_gauss_mix_pdf(th, ycross[m], ymin, ymax)
    }
  }

  # (3-4) Create the vector of critical points, calculate densities, and
  #       identify peak
  ycrit <- c(ymin, ycross, ymax)
  fcrit <- c(bd_calc_gauss_mix_pdf(th, ymin, ymin, ymax), fcross, bd_calc_gauss_mix_pdf(th, ymin, ymin, ymax))
  indPeak <- which.max(fcrit)
  ypeak <- ycrit[indPeak]
  fpeak <- fcrit[indPeak]

  # (5) Create ylo, yhi, and slope
  numPer <- length(ycrit) - 1 # Number of periods
  ylo <- ycrit[1:numPer]
  yhi <- ycrit[2:(numPer + 1)]
  df <- diff(fcrit)
  slope <- rep("pos", numPer)
  slope[df < 0] <- "neg"

  # (6) Create the pattern (then return the result)
  pattern <- c(slope, as.character(indPeak))

  return(list(periods = data.frame(ylo = ylo, yhi = yhi, slope = slope), indPeak = indPeak, ypeak = ypeak, fpeak = fpeak, pattern = pattern))
}
