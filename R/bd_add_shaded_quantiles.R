#' @title Add a shaded region showing the min/max quantiles
#'
#' @description Add a shaded region showing the min/max quantiles (usually +/-2.5%) to a plot
#'
#' @details The input, an, is the result of a call to bd_analyze_soln. It is a
#' list-like object of class bayDem_analysis with information on the quantiles
#' of the density function and growth rate. Typically, the min/max quantiles
#' are 2.5% and 97.5%, but this can be changed via the input lev in the call to
#' bd_analyze_soln.
#'
#' @param an A list-like object of class \code{baydem::bd_analysis} with information on the quantiles of the density function and growth rate
#' @param add [default FALSE] Whether to make a new plot or add to the active plot
#' @param ... Additional parameters to pass to polygon
#'
#' @export
bd_add_shaded_quantiles <- function(an,col=adjustcolor("grey",alpha.f=0.5),...) {

  # The indidces of the minimum and maximum quantiles
  indMin <- which.min(an$probs)
  indMax <- which.max(an$probs)
  polygon(c(an$t,rev(an$t)),c(an$Qdens[indMin,],rev(an$Qdens[indMax,])),border=NA,xlab=NULL,col=col)
}
