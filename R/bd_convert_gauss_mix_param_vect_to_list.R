#' @title Convert a Gaussian mixture parameterization from a vector to a list
#'
#' @description For the input, parameterized Gaussian mixture convert from a
#'              vector representation (which is more efficient for calculations
#'              and the form expected by, e.g., bayDem_calcGaussMixPdf.R) list
#'              (returned by stan).
#'
#' @param  th A vector that parameterizes the Gaussian mixture. See
#'             bayDem_calcGaussMixPdf.R
#' @return samp A list that parameterizes the Gaussian mixture. See  See
#' @export

bd_convert_gauss_mix_param_vect_to_list <- function(th) {
  K <- length(th) / 3 # Number of  mixtures
  samp <- list(pi = th[1:K], mu = th[(1 + K):(2 * K)], sig = th[(1 + 2 * K):(3 * K)])
  return(samp)
}
