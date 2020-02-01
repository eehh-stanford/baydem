#' @title For each sample, calculate the time it takes for the density to decrease by half from the peak
#'
#' @details
#' For each sample, calculate the time it takes for the density to decrease by
#' half from the peak. Optionally, a different proportion can be used than the
#' default propChange = 0.5. For example, with propChange = 0.1 the time it
#' takes for the density to decrease by 10% is used. If the relative density is
#' not reached, the half life for the sample is set to NA. If there is no
#' interior peak in the range peakRange, which is taumin to taumax by default,
#' the half life is set to NA.
#'
#' @param soln The result of a call to bd_do_inference
#' @param propChange (Default 0.5) The relative decrease in density to use for the duration calculation
#' @param anal (Optional) The result of a call to bd_analyze_soln. If not provided, it is calculated
#' @param peakRange (default: `c(taumin, taumax)`) peakRange can be given so that the peak density used is on the range peakRange
#'
#' @return A vector of "half-lives" (proportional change set by propChange)
#'
#' @export
bd_calc_half_life_from_peak <- 
  function(soln, 
           propChange = 0.5, 
           anal = NA, 
           peakRange = NA) {
  TH <- bd_extract_param(soln$fit)
  N <- nrow(TH)
  taumin <- soln$prob$hp$taumin
  taumax <- soln$prob$hp$taumax
  dtau <- soln$prob$hp$dtau

  if (all(is.na(anal))) {
    anal <- bd_analyze_soln(soln)
  }
  summList <- anal$summList

  if (all(is.na(peakRange))) {
    peakRange <- c(taumin, taumax)
  }


  halfLife <- rep(NA, N)
  for (n in 1:N) {
    th <- TH[n, ]
    # Identify the peak, ensuring it is on peakRange

    # critical points
    tcrit <- c(summList[[n]]$periods$tlo, summList[[n]]$periods$thi[length(summList[[n]]$periods$thi)])
    tcrit <- tcrit[peakRange[1] <= tcrit & tcrit <= peakRange[2]]
    fcrit <- bd_calc_gauss_mix_pdf(th, tcrit, taumin, taumax)
    indPeak <- which.max(fcrit)
    tpeak <- tcrit[indPeak]
    fpeak <- fcrit[indPeak]
    isIn <- taumin < tpeak && tpeak < taumax
    if (isIn) {
      # Function for root finder
      rootFun <- function(t) {
        return(fpeak * propChange - bd_calc_gauss_mix_pdf(th, t, taumin, taumax, type = "density"))
      }

      # Find root. Catch any errors in case the half life does not exist on the
      # interval tpeak to taumax
      result <- tryCatch(
        {
          root <- stats::uniroot(rootFun, 
                                 lower = tpeak, 
                                 upper = peakRange[2])
          halfLife[n] <- min(root$root - tpeak)
        },
        error = function(e) {
          NA
        }
      )
    }
  }
  return(halfLife)
}
