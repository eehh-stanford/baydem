% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bayesian_inference_functions.R
\name{summarize_bayesian_inference}
\alias{summarize_bayesian_inference}
\title{Calculate some key summary measures using the result of a call to
\code{sample_theta}}
\usage{
summarize_bayesian_inference(
  bayesian_soln,
  rc_meas,
  density_model,
  calib_df,
  dtau = 5,
  th_sim = NA,
  lev = 0.025,
  rate_prop = NA,
  do_sample_summaries = T
)
}
\arguments{
\item{bayesian_soln}{The solution, a list-like object of class
bd_bayesian_soln (see sample_theta).}

\item{rc_meas}{The radiocarbon measurements (see import_rc_data).}

\item{density_model}{The density model (see set_density_model).}

\item{calib_df}{The calibration data frame (see load_calib_curve).}

\item{dtau}{The spacing of the sampling grid (default: 5).}

\item{th_sim}{The known parameters used to create simulation data (default:
NA, not provided).}

\item{lev}{The level to use for the quantile bands (default: 0.025).}

\item{rate_prop}{The cumulative density needed to define rate growth bands
(default: NA, not used).}

\item{do_sample_summaries}{Whether to calculate some summary information for
each sampled curve (Default: TRUE).}
}
\value{
A list with information on the quantiles of the density function and
growth rate (and sample summaries)
}
\description{
\code{sample_theta} calls Stan to do Bayesian inference by
generating a sample of parameters from the posterior of theta (or \code{th}).
\code{sample_theta} analyzes the result of that inference. Notably,
it calculates the quantiles of the density function and the growth rate.
}
\details{
\code{bayesian_soln} is the result of a call to
\code{sample_theta}. It contains posterior samples for the density
model. The primary thing \code{summarize_bayesian_inference} does is
calculate quantiles of both the parameterized density and growth rate. For
example, for a calendar date tau_g each sample yields a density and growth
rate. The quantile is the value of the density or growth rate such that a
given proportion of samples are smaller than that value. The probabilities
used to calculate these quantiles are \code{probs = c(lev, 0.5, 1-lev)}, where
\code{lev} is the level (0.025 by default, so that 95\% of the observations lie
between the first and last quantile bands).

In addition, \code{summarize_bayesian_inference} identifies calendar dates
for which the growth rate quantiles defined by \code{lev} and \code{1 - lev} do not
contain zero. This indicates significant positive or negative growth for the
density curve. The output vector \code{growth_state} codes calendar dates by
growth state as 'negative', 'zero', and 'positive'. For the Gaussian mixture
parameterization of the density, the rate is not typically meaningful near
the calendar date boundaries where it increases linearly as the calendar date
goes to positive or negative infinity. The parameter \code{rate_prop} provides
control on how calendar dates are classified by growth rate near these
boundaries. In particular, the calendar dates with a cumulative density (50\%
quantile) below \code{rate_prop} (for the lower boundary) or above \code{1 - rate_prop}
(for the upper boundary) are classified as 'missing' in \code{growth_state}. By
default, \code{rate_prop} is NA and no calendar dates are classified as missing.

By default, a summary is done for each sample by calling summarize_sample.
This is not done if do_summary is FALSE.
}
\seealso{
\itemize{
\item \code{\link[=import_rc_data]{import_rc_data()}} for the format of \code{rc_meas}
\item \code{\link[=set_density_model]{set_density_model()}} for the format of \code{density_model}
\item \code{\link[=load_calib_curve]{load_calib_curve()}} for the format of \code{calib_df}
}
}
