globalVariables(c("mu", "sig"))

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
