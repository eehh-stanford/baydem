#' @title Plot the summed probability density function (SPDF)
#'
#' @description Plot the summed probability density function (SPDF)
#'
#' @details The input, an, is the result of a call to bd_analyze_soln. It is a
#' list-like object of class bayDem_analysis with information on the quantiles
#' of the density function and growth rate. By default, a new plot is made, but
#' if add = TRUE the curve is added to the active plot.
#'
#' @param an A list-like object of class \code{baydem::bd_analysis} with information on the quantiles of the density function and growth rate
#' @param add (default: `FALSE`) Whether to make a new plot or add to the active plot
#' @param ... Additional parameters to pass to plot / lines
#'
#' @export
bd_plot_summed_density <- function(an, add = F, ...) {
  if (add) {
    graphics::lines(an$y, an$f_spdf, ...)
  } else {
    graphics::plot(an$y, an$f_spdf, type = "l", ...)
  }
}
