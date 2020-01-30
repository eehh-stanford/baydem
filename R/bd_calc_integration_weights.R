#' @title Calculate integration weights assuming the midpoint rule
#'
#' @description y is a vector of locations where a function to be integrated is
#'              evaluated (the function values, f, are not input). Let
#'              y_g be the locations of integration and f_g the corresponding
#'              function values for g=1, 2, ... G. Assuming trapezoidal
#'              integration at the midpoints between elements of y, the weghts
#'              to use for integration are dy_g = (y_(g+1) - y_(g-1)) / 2, where
#'              the conventions y_0 = y_1 and y_(G+1) = y_G are used.
#'
#' @param y A vector of locations where the function is sampled, possibly irregular
#'
#' @return A vector of integration weights the same length as y

#' @export
bd_calc_integration_weights <- function(y) {
  G <- length(y)
  weightVect <- rep(NA, length(y))
  indCent <- 2:(G - 1)
  weightVect[indCent] <- (y[indCent + 1] - y[indCent - 1]) / 2
  weightVect[1] <- (y[2] - y[1]) / 2
  weightVect[G] <- (y[G] - y[G - 1]) / 2
  return(weightVect)
}
