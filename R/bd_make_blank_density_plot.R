#' @title Prepare a blank plot for showing densities
#'
#' @description Prepare a blank plot for showing densities
#'
#' @details The input, an, is the result of a call to bd_analyze_soln. It is a
#' list-like object of class bayDem_analysis with information on the quantiles
#' of the density function and growth rate. Prepare a blank plot for adding
#' density information, e.g. by calling bd_plot_50_percent_quantile and
#' bd_add_shaded_quantiles.R. If the axis labels and limits are not specified,
#' sensible defaults are used.
#'
#' @param an A list-like object of class \code{baydem::bd_analysis} with information on the quantiles of the density function and growth rate
#' @param ... Additional parameters to pass to plot
#'
#' @export
bd_make_blank_density_plot <- function(an,xlim=NA,ylim=NA,xlab=NA,ylab=NA,...) {

  if(is.na(xlab)) {
    xlab <- 'Calendar Date [AD]'
  }

  if(is.na(ylab)) {
    ylab <- 'Density'
  }

  if(all(is.na(xlim))) {
    xlim <- range(an$tau)
  }

  if(all(is.na(ylim))) {
    # The following command works even if, e.g., f_sim is not in the list an
    # because an$f_sim results in NULL
    ylim <- c(0,max(c(max(an$Qdens),an$f_spdf,an$f_sim)))
  }

  plot(NULL, xlim=xlim, ylim=ylim, ylab=ylab, xlab=xlab,...)
}
