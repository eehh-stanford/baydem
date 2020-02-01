#' @title Calculate the density of a possibly truncated Gaussian mixture for a matrix of samples
#'
#' @description TH is a matrix of samples of a possibly truncated Gaussian
#'              mixture with dimensions S x P, where S is the number of samples
#'              and P is the number of parameters. Repeatedly call
#'              \code{bd_calc_gauss_mix_pdf} to calculate the density function
#'              for each sample at the points in the vector tau, which has length
#'              G. The output density matrix, fMat, has dimensions S x G.
#'              \code{bd_calc_gauss_mix_pdf_mat} supports the same optional
#'              as \code{bd_calc_gauss_mix_pdf}: taumin / taumax to specify the
#'              boundaries of truncation and type to set whether the density,
#'              cumuluative distribution, derivative, or rate is calculated.
#'
#' @param TH A matrix of samples with dimensions S x P
#' @param tau A vector of points for the density calculation with length G
#' @param taumin (optional) The mininum time-value for truncation
#' @param taumax (optional) The maximum time-value for truncation
#' @param type (optional) The type of calculation: density (default), cumulative, derivative, or rate
#'
#' @return The output matrix with dimensions S x G
#'
#' @export

bd_calc_gauss_mix_pdf_mat <- function(TH, tau, taumin = NA, taumax = NA, type = "density") {
  S <- dim(TH)[1] # number of samples
  G <- length(tau) # number of time-values (usually grid points)
  fMat <- matrix(NA, S, G)

  # Iterate over samples to fill fMat
  for (s in 1:S) {
    fMat[s, ] <- bd_calc_gauss_mix_pdf(TH[s, ], tau, taumin = taumin, taumax = taumax, type = type)
  }
  return(fMat)
}
