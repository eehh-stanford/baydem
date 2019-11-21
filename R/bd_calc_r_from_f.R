# Description
# Calculate the growth rate, r, given the probability density f.
#
# Example calls(s)
#
#   r <- bd_calc_r_from_f(f,dy=1)
#
# Input(s)
#   Name    Type           Description
#   f       vector         [G x 1] The probability density function
#   dy      scalar         [default=1] The spacing of the vector f
#
# Output(s)
#   Name    Type           Description
#   r       vector         [(G-1) x 1] The growth rate

bd_calc_r_from_f <- function(f, dy = 1) {
  G <- length(f)
  r <- log(f[2:G] / f[1:(G - 1)]) / dy
  return(r)
}
