#' @title Convert a Gaussian mixture parameterization from a list to a vector
#'
#' @description For the input, parameterized Gaussian mixture convert from a
#'              list (returned by stan) to a vector representation (which is
#'              more efficient for calculations and the form expected by, e.g.,
#'              bayDem_calcGaussMixPdf.R).
#'
#' @param samp A list that parameterizes the Gaussian mixture. See
#'             bayDem_samplePrior.R.
#' @return  th A vector that parameterizes the Gaussian mixture. See
#'             bayDem_calcGaussMixPdf.R
#' @export
bd_convert_gauss_mix_param_list_to_vect <- function(samp) {
  K <- length(samp$sig)
  th <- c(samp$pi, samp$mu, samp$sig)
  return(th)
}
