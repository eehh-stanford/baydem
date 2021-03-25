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

#' @title Caculate the density of a possibly truncated Gaussian mixture
#'
#' @description \code{bd_calc_gauss_mix_pdf} calculates the probability density
#' for a possibly truncated Gaussian mixture. The density, f', and rate, f'/f,
#' can also be calculated by specifying the input \code{type}, which is
#' 'density' by default, but can also be 'cumulative', 'derivative', or 'rate'.
#' The parameter vector th has the ordering (pik, muk, sigk), where
#' pik, muk, and sigk are the weight, mean, and standard deviation of the k-th
#' mixture. There are K total mixtures. The truncation boundaries taumin and
#' taumax are optional.
#'
#' @param th A vector-like object that parameterizes the distribution
#' @param tau Vector of locations at which to calculate the density
#' @param taumin Minimum truncation bound
#' @param taumax Maximum truncation bound
#' @param type (default density) An optional input specifying whether to
#'        calculate the density, cumulative distribution, derivative of the
#'        density, or rate
#'
#' @return The output vector with length G
#'
#' @export
bd_calc_gauss_mix_pdf <- function(th,
                                  tau,
                                  taumin = NA,
                                  taumax = NA,
                                  type = "density") {
  # First, determine whether taumin and taumax are input
  doNorm <- !is.na(taumin)

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
    doNorm <- F
  }

  # Replicate the inputs for efficient calculations
  G <- length(tau)
  tauRep <- rep(tau, K)
  piRep <- as.vector(t(matrix(th[1:K], G, nrow = K)))
  muRep <- as.vector(t(matrix(th[(K + 1):(2 * K)], G, nrow = K)))
  sigRep <- as.vector(t(matrix(th[(2 * K + 1):(3 * K)], G, nrow = K)))

  # Do the calculation
  if (type == "density") {
    output <- rowSums(matrix(stats::dnorm(tauRep, muRep, sigRep) * piRep, ncol = K))
  } else if (type == "cumulative") {
    output <- rowSums(matrix(stats::pnorm(tauRep, muRep, sigRep) * piRep, ncol = K))
  } else if (type == "derivative") {
    output <- rowSums(matrix(miscTools::ddnorm(tauRep, muRep, sigRep) * piRep, ncol = K))
  } else if (type == "rate") {
    output1 <- rowSums(matrix(miscTools::ddnorm(tauRep, muRep, sigRep) * piRep, ncol = K)) # the numerator
    output2 <- rowSums(matrix(stats::dnorm(tauRep, muRep, sigRep) * piRep, ncol = K)) # the denominator
    output <- output1 / output2
  } else { # This should not be reached. Throw an error just in case
    stop("type must be density, cumulative, derivative, or rate")
  }
  if (doNorm) {
    normFact <- diff(bd_calc_gauss_mix_pdf(th, c(taumin, taumax), type = "cumulative"))
    output <- output / normFact
  }
  return(output)
}

#' @title Create a possibly truncated Gaussian mixture
#'
#' @details `x` is a [`data.frame`] that specifies the Gaussian mixture.
#' Each row represents a gaussian with a mean of `mu`,  a standard deviation of `sig`,
#' and a probability for the mixing component of `pi`. Optionally, the samples are
#' drawn on the truncated interval `taumin` to `taumax`.
#'
#' @param x Parameterization [`data.frame`] for the Gaussian mixture
#' @param taumin (Optional) Lower bound for samples
#' @param taumax (Optional) Upper bound for samples
#'
#' @return An object of class [`distr::UnivarMixingDistribution`]
#' @importFrom magrittr `%$%` `%<>%`
#' @export
bd_create_gauss_mix <-
  function(x,
           taumin = -Inf,
           taumax = Inf) {
    suppressWarnings(
      x %>%
        dplyr::rowwise() %>%
        dplyr::mutate(norm = list(distr::Norm(mean = mu, sd = sig))) %$%
        distr::UnivarMixingDistribution(
          Dlist = norm,
          mixCoeff = pi
        ) %>%
        distr::Truncate(
          lower = taumin,
          upper = taumax
        )
    )
  }

#' @title Create a list of possibly truncated Gaussian mixtures for a matrix of samples
#'
#' @description `x` is a matrix of samples of a possibly truncated Gaussian
#'              mixture with dimensions S x P, where S is the number of samples
#'              and P is the number of parameters.
#'              \code{bd_create_gauss_mix_mat} supports the same optional
#'              as \code{bd_calc_gauss_mix}: taumin / taumax to specify the
#'              boundaries of truncation.
#'
#' @param x A matrix of samples with dimensions S x P
#' @param taumin (optional) The mininum time-value for truncation
#' @param taumax (optional) The maximum time-value for truncation
#' @param progress (optional) Whether to show a progress bar
#'
#' @return A list of objects of class [`distr::UnivarMixingDistribution`]
#'
#' @export

bd_create_gauss_mix_mat <- function(x,
                                    taumin = -Inf,
                                    taumax = Inf,
                                    progress = TRUE) {
  x %<>%
    tibble::as_tibble() %>%
    tidyr::unite("pi", dplyr::starts_with("pi")) %>%
    tidyr::unite("mu", dplyr::starts_with("mu")) %>%
    tidyr::unite("sig", dplyr::starts_with("sig")) %>%
    dplyr::mutate_all(
      function(x) {
        stringr::str_split(x, pattern = "_") %>%
          purrr::map(as.numeric)
      }
    ) %>%
    purrr::transpose() %>%
    purrr::map(tibble::as_tibble) %>%
    purrr::map(~ dplyr::mutate(.x, pi = pi / sum(pi)))

  pb <- dplyr::progress_estimated(length(x))

  x %>%
    purrr::map(function(x) {
      if (progress) {
        pb$tick()$print()
      }

      x %>%
        bd_create_gauss_mix(
          taumin = taumin,
          taumax = taumax
        )
    })
}
