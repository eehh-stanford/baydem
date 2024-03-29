% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bayesian_inference_functions.R
\name{sample_theta}
\alias{sample_theta}
\title{Sample from the posterior of a density model for a set of radiocarbon
measurements}
\usage{
sample_theta(
  rc_meas,
  density_model,
  hp,
  calib_df,
  th0 = NA,
  init_seed = NA,
  stan_seed = NA,
  calibration_curve = "intcal20",
  control = list()
)
}
\arguments{
\item{rc_meas}{The radiocarbon measurements (see import_rc_data).}

\item{density_model}{The density model (see set_density_model).}

\item{hp}{Hyperparameters for the priors and to specify the spacing of the
Riemann sum that approximates the integral for the likelihood.}

\item{calib_df}{The calibration data frame (see load_calib_curve).}

\item{th0}{An optional parameter vector to initialize the Stan chains. If not
provided, it is set by calling init_trunc_gauss_mix.}

\item{init_seed}{An optional random number seed for determining the starting
parameter vector using a maximum likelihood fit. If not provided, it is
drawn. It should not be provided if th0 is provided.}

\item{stan_seed}{An optional random number seed for the call to Stan. If not
provided, it is drawn.}

\item{calibration_curve}{The calibration curve to use for masking (only used
control$mask is TRUE). The default is "intcal20". Other options are
"shcal20" and "marine20". For further options see Bchron::BchronCalibrate.}
}
\value{
\code{bayesian_soln}, a list-like object of class bd_bayesian_soln with the
following fields:
\itemize{
\item{\code{fit}}
{The result of the call to stan}
\item{\code{final_th0}}
{The final \code{th0} value; i.e., never NA.}
\item{\code{final_init_seed}}
{The final init_seed value; i.e., never NA unless \code{th0} is
provided.}
\item{\code{final_stan_seed}}
{The final \code{stan_seed} value; i.e., never NA.}
\item{\code{final_control}}
{The final control parameters used; i.e., if a parameter is not
provided.}
\item{\code{optional_inputs}}
{A record of the actual input values for the optional inputs, which are
\code{th0}, \code{init_seed}, \code{stan_seed}, and \code{control}.}
}
}
\description{
This is the core function that implements the Bayesian inference. Currently,
the only supported density model is a truncated Gaussian mixture. If a
starting parameter vector (\code{th0}) is not provided, it is set by calling
init_trunc_gauss_mix; the same vector is used for all sampling chains. Named
elements of the variable control must consist of one of the following four
options (defaults in parentheses):

\itemize{
\item{\code{num_chains}}

The calibration_curve to use for masking is separately input to maintain

consistency with previous versions of baydem.
{Number of chains (4)}
\item{\code{samps_per_chain}}
{Number of samples per chain (2000)}
\item{\code{warmup}}
{Number of warmup samples (\code{samps_per_chain/2})}
\item{\code{stan_control}}
{Additional control parameters to pass to stan (\code{list()})}
\item{\code{mask}}
{Whether to mask the likelihood sum based on individual calibration
(FALSE)}
}

The calibration_curve to use for masking is separately input to maintain
consistency with previous versions of baydem.
}
\seealso{
\itemize{
\item \code{\link[=import_rc_data]{import_rc_data()}} for the format of \code{rc_meas}
\item \code{\link[=set_density_model]{set_density_model()}} for the format of \code{density_model}
\item \code{\link[=load_calib_curve]{load_calib_curve()}} for the format of \code{calib_df}
}
}
