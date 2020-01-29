#' @title Caculate the density of a possibly truncated Gaussian mixture
#'
#' @description \code{bd_calc_gauss_mix_pdf} calculates the probability density
#' for a possibly truncated Gaussian mixture. The density, f', and rate, f'/f,
#' can also be calculated by specifying the input \code{type}, which is
#' 'density' by default, but can also be 'cumulative', 'derivative', or 'rate'.
#' The parameter vector th has the ordering [wk, muk, sigk], where
#' wk, muk, and sigk are the weight, mean, and standard deviation of the k-th
#' mixture. There are K total mixtures. The truncation boundaries ymin and ymax
#' are optional.
#'
#' @param y Vector of locations at which to calculate the density
#' @param th A vector-like object that parameterizes the distribution
#' @param type (default density) An optional input specifying whether to
#'        calculate the density, cumulative distribution, derivative of the
#'        density, or rate
#'
#' @return The output vector with length G
#'
#' @export
bd_calc_gauss_mix_pdf <- function(th, y, ymin = NA, ymax = NA, type = "density") {
  # First, determine whether ymin and ymax are input
  doNorm <- !is.na(ymin)

  # Determine number of mixtures (and also do error checking on length(th))
  if (length(th) %% 3 == 0) {
    # Correct length
    K <- length(th) / 3
  } else {
    stop("Unsupported length of parameter vector th")
  }

  # Check that type is valid
  if (!type %in% c("density", "cumulative", "derivative", "rate")) {
    stop("type must be density, cumulative, derivative, or rate")
  }

  # If type is rate, no normalization is needed.
  # Hence, if type is rate then ymin and ymax are ignored if they are input.
  if (type == "rate") {
    doNorm <- F
  }

  # Replicate the inputs for efficient calculations
  G <- length(y)
  yRep <- rep(y, K)
  wRep <- as.vector(t(matrix(th[1:K], G, nrow = K)))
  muRep <- as.vector(t(matrix(th[(K + 1):(2 * K)], G, nrow = K)))
  sigRep <- as.vector(t(matrix(th[(2 * K + 1):(3 * K)], G, nrow = K)))

  # Do the calculation
  if (type == "density") {
    output <- rowSums(matrix(stats::dnorm(yRep, muRep, sigRep) * wRep, ncol = K))
  } else if (type == "cumulative") {
    output <- rowSums(matrix(stats::pnorm(yRep, muRep, sigRep) * wRep, ncol = K))
  } else if (type == "derivative") {
    output <- rowSums(matrix(miscTools::ddnorm(yRep, muRep, sigRep) * wRep, ncol = K))
  } else if (type == "rate") {
    output1 <- rowSums(matrix(miscTools::ddnorm(yRep, muRep, sigRep) * wRep, ncol = K)) # the numerator
    output2 <- rowSums(matrix(stats::dnorm(yRep, muRep, sigRep) * wRep, ncol = K)) # the denominator
    output <- output1 / output2
  } else { # This should not be reached. Throw an error just in case
    stop("type must be density, cumulative, derivative, or rate")
  }
  if (doNorm) {
    normFact <- diff(bd_calc_gauss_mix_pdf(th, c(ymin, ymax), type = "cumulative"))
    output <- output / normFact
  }
  return(output)
}
