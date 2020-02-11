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
