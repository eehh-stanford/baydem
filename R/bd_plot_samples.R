# Description
# Add priors to a plot given the input hyperparameters and number of priors to
# add.
#
# Example calls(s)
#
# bd_plot_samples(hp,N,y)
#
# Input(s)
#   Name    Type           Description
#   samps   list           List of samples to plot
#   y       vector         Points at which to evaluate curves
#   rPlot   boolean        [optional] TRUE if the plot is the rate and FALSE if
#                          the plot is the probability density
#
# Output(s)
#   Name    Type           Description
#   NA      NA             NA
bd_plot_samples <- function(samps, y, rPlot = F, add = F, lineCol = "black", lineWid = 1) {
  N <- length(samps)
  for (n in 1:N) {
    th <- bd_convert_gauss_mix_param_list_to_vect(samps[[n]])
    f <- bd_calc_gauss_mix_pdf(th, y)
    if (rPlot) {
      r <- bd_calc_r_from_f(f, y[2] - y[1])
      if (add) {
        graphics::lines(y[1:(length(f) - 1)], r, col = lineCol, lwd = lineWid)
      } else {
        graphics::plot(y[1:(length(f) - 1)], r, type = "l", col = lineCol, lwd = lineWid)
      }
    } else {
      if (add) {
        graphics::lines(y, f, col = lineCol, lwd = lineWid)
      } else {
        graphics::plot(y, f, type = "l", col = lineCol, lwd = lineWid)
      }
    }
  }
}
