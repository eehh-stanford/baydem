#' @title For each sample, calculate the time it takes for the density to decrease by half from the peak
#'
#' @details
#' For each sample, calculate the time it takes for the density to decrease by
#' half from the peak. Optionally, a different proportion can be used than the
#' default propChange = 0.5. For example, with propChange = 0.1 the time it
#' takes for the density to decrease by 10% is used. If the relative density is
#' not reached, the half life for the sample is set to NA. If there is no
#' interior peak in the range peakRange, which is ymin to ymax by default, the
#' half life is set to NA.
#'
#' @param soln The result of a call to bd_do_inference
#' @param anal (Optional) The result of a call to bd_analyze_soln. If not provided, it is calculated
#' @param propChange (Default 0.5) The relative decrease in density to use for the duration calculation
#' @param (default ymin to ymax) peakRange can be given so that the peak density used is on the range peakRange
#'
#' @return A vector of half-lives
#'
#' @export
bd_calc_half_life_from_peak <- function(soln, propChange = 0.5, anal = NA, peakRange = NA) {
  TH <- bd_extract_param(soln$fit)
  N <- nrow(TH)
  ymin <- soln$prob$hp$ymin
  ymax <- soln$prob$hp$ymax
  dy <- soln$prob$hp$dy

  if (all(is.na(anal))) {
    anal <- bd_analyze_soln(soln)
  }
  summList <- anal$summList

  if (all(is.na(peakRange))) {
    peakRange <- c(ymin, ymax)
  }


  halfLife <- rep(NA, N)
  for (n in 1:N) {
    th <- TH[n, ]
    # Identify the peak, ensuring it is on peakRange

    # critical points
    ycrit <- c(summList[[n]]$periods$ylo, summList[[n]]$periods$yhi[length(summList[[n]]$periods$yhi)])
    ycrit <- ycrit[peakRange[1] <= ycrit & ycrit <= peakRange[2]]
    fcrit <- bd_calc_gauss_mix_pdf(th, ycrit, ymin, ymax)
    indPeak <- which.max(fcrit)
    ypeak <- ycrit[indPeak]
    fpeak <- fcrit[indPeak]
    isIn <- ymin < ypeak && ypeak < ymax
    if (isIn) {
      # Function for root finder
      rootFun <- function(y) {
        return(fpeak * propChange - bd_calc_gauss_mix_pdf(th, y, ymin, ymax, type = "density"))
      }

      # Find root. Catch any errors in case the half life does not exist on the
      # interval ypeak to ymax
      result <- tryCatch({
        root <- uniroot(rootFun, lower = ypeak, upper = peakRange[2])
        halfLife[n] <- min(root$root - ypeak)
      })
    }
  }
  return(halfLife)
}
