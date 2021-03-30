#' @title Calculate the density of a possibly truncated Gaussian mixture for a
#' matrix of samples
#'
#' @description TH is a matrix of samples of a possibly truncated Gaussian
#'              mixture with dimensions S x P, where S is the number of samples
#'              and P is the number of parameters. Repeatedly call
#'              \code{calc_gauss_mix_pdf} to calculate the density function
#'              for each sample at the points in the vector tau, which has length
#'              G. The output density matrix, f_mat, has dimensions S x G.
#'              \code{calc_gauss_mix_pdf_mat} supports the same optional
#'              as \code{calc_gauss_mix_pdf}: tau_min / tau_max to specify the
#'              boundaries of truncation and type to set whether the density,
#'              cumuluative distribution, derivative, or rate is calculated.
#'
#' @param TH A matrix of samples with dimensions S x P
#' @param tau A vector of points for the density calculation with length G
#' @param tau_min (optional) The mininum time-value for truncation
#' @param tau_max (optional) The maximum time-value for truncation
#' @param type (optional) The type of calculation: density (default),
#'   cumulative, derivative, or rate
#'
#' @return The output matrix with dimensions S x G
#'
#' @export

calc_gauss_mix_pdf_mat <- function(TH, tau, tau_min = NA, tau_max = NA, type = "density") {
  S <- dim(TH)[1] # number of samples
  G <- length(tau) # number of time-values (usually grid points)
  f_mat <- matrix(NA, S, G)

  # Iterate over samples to fill fMat
  for (s in 1:S) {
    f_mat[s, ] <- calc_gauss_mix_pdf(TH[s, ],
                                     tau,
                                     tau_min = tau_min,
                                     tau_max = tau_max, type = type)
  }
  return(f_mat)
}

#' @title Caculate the density of a possibly truncated Gaussian mixture
#'
#' @description \code{calc_gauss_mix_pdf} calculates the probability density
#' for a possibly truncated Gaussian mixture. The density, f', and rate, f'/f,
#' can also be calculated by specifying the input \code{type}, which is
#' 'density' by default, but can also be 'cumulative', 'derivative', or 'rate'.
#' The parameter vector th has the ordering (pi_k, mu_k, s_k), where
#' pi_k, mu_k, and s_k are the weight, mean, and standard deviation of the k-th
#' mixture. There are K total mixtures. The truncation boundaries tau_min and
#' tau_max are optional.
#'
#' @param th A vector-like object that parameterizes the distribution
#' @param tau Vector of locations at which to calculate the density
#' @param tau_min Minimum truncation bound
#' @param tau_max Maximum truncation bound
#' @param type (default density) An optional input specifying whether to
#'        calculate the density, cumulative distribution, derivative of the
#'        density, or rate
#'
#' @return The output vector with length G
#'
#' @export
calc_gauss_mix_pdf <- function(th,
                               tau,
                               tau_min = NA,
                               tau_max = NA,
                               type = "density") {
  # First, determine whether taumin and taumax are input
  do_norm <- !is.na(tau_min)

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
  # Hence, if type is rate then taumin and taumax are ignored if they are input.
  if (type == "rate") {
    do_norm <- F
  }

  # Replicate the inputs for efficient calculations
  G <- length(tau)
  tau_rep <- rep(tau, K)
  pi_rep <- as.vector(t(matrix(th[1:K], G, nrow = K)))
  mu_rep <- as.vector(t(matrix(th[(K + 1):(2 * K)], G, nrow = K)))
  s_rep <- as.vector(t(matrix(th[(2 * K + 1):(3 * K)], G, nrow = K)))

  # Do the calculation
  if (type == "density") {
    output <- rowSums(matrix(
      stats::dnorm(tau_rep, mu_rep, s_rep) * pi_rep, ncol = K))
  } else if (type == "cumulative") {
    output <- rowSums(matrix(
      stats::pnorm(tau_rep, mu_rep, s_rep) * pi_rep, ncol = K))
  } else if (type == "derivative") {
    output <- rowSums(matrix(
      miscTools::ddnorm(tau_rep, mu_rep, s_rep) * pi_rep, ncol = K))
  } else if (type == "rate") {
    # the numerator
    output1 <- rowSums(matrix(
      miscTools::ddnorm(tau_rep, mu_rep, s_rep) * pi_rep, ncol = K))
    output2 <- rowSums(matrix(
      # the denominator
      stats::dnorm(tau_rep, mu_rep, s_rep) * pi_rep, ncol = K))
    output <- output1 / output2
  } else { # This should not be reached. Throw an error just in case
    stop("type must be density, cumulative, derivative, or rate")
  }
  if (do_norm) {
    norm_fact <- diff(calc_gauss_mix_pdf(th, c(tau_min, tau_max), type = "cumulative"))
    output <- output / norm_fact
  }
  return(output)
}