#' @title Calculate the quantiles for an input matrix X
#'
#' @description The input matrix X has dimensions S x G, where S is the number
#'              of samples and G the number of grid points at which X was
#'              evaluated. Calculate quantiles for each grid point, g = 1,2,..G.
#'
#' @param X The matrix for which quantiles are calculated, with dimensions S x G
#' @param probs The probability values at which to calculate the quantiles (default: `c(0.025, 0.5, 0.975)`)
#'
#' @return The quantiles, a matrix with dimension length(probs) x G
#'
#' @export
bd_calc_quantiles <- function(X, probs = c(.025, .5, .975)) {
  numQuant <- length(probs) # Number of quantiles
  G <- dim(X)[2] # Number of grid points

  Q <- matrix(NA, numQuant, G) # Initialize Q with dimensions numQuant x G
  # Iterate over grid points to calculate quantiles
  for (g in 1:G) {
    Q[, g] <- stats::quantile(X[, g], probs = probs)
  }
  return(Q)
}
