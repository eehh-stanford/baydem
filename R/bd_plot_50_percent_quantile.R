#' @title Plot the 50% quantile curve
#'
#' @description Plot the 50% quantile curve
#'
#' @details The input, an, is the result of a call to bd_analyze_soln. It is a
#' list-like object of class bayDem_analysis with information on the quantiles
#' of the density function and growth rate. By default, a new plot is made, but
#' if add = TRUE the curve is added to the active plot.
#'
#' @param an A list-like object of class \code{baydem::bd_analysis} with information on the quantiles of the density function and growth rate
#' @param add [default FALSE] Whether to make a new plot or add to the active plot
#' @param ... Additional parameters to pass to plot / lines
#'
#' @export
bd_plot_50_percent_quantile <- function(an,add=F,...) {

  # The index of the 50% quantile
  ind50 <- which(an$probs == 0.5) # The index in probs / Qdens of the 50% quantile
  if(add) {
    lines(an$y,an$Qdens[ind50,],...)
  } else {
    plot(an$y,an$Qdens[ind50,],type='l',...)
  }
}
