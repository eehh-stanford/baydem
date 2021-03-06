% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulation_functions.R
\name{simulate_rc_data}
\alias{simulate_rc_data}
\title{Create simulated radiocarbon data}
\usage{
simulate_rc_data(sim_spec)
}
\arguments{
\item{sim_spec}{A simulation specification (see Description)}
}
\value{
A list consisting of the input sim_spec and data (see Description)
}
\description{
The input is a multi-tiered list specifying all aspects of the simulation. In
particular, sim_spec (for simulation specification) consists of three
required fields and one optional field: model_spec (model specification),
N (number of samples), calib_curve (calibration curve), and seed (an optional
random number seed). model_spec, in turn, is a list than consists of three
required fields: density_type (the type of parametric model to use for the
target density), th (the parameter vector for the target density), is_AD (a
boolean variable indicating whether the target density is for years AD or
BP), and error_spec (a specification for how to model the measurement
errors). Currently, two model types are supported for the target density,
gauss_mix = a Gaussian mixture and trunc_gauss_mix = a truncated Gaussian
mixture. Currently, one specification of the error_spec (a list) is
supported: unif_fm = a uniform draw for the error of the fraction modern
value. For unif_fm, error_spec must contain the fields min and max (the
minimum and maximum values for the uniform draw).

The output is a list with the fields dates and rc_meas. Dates is a length N
vector of dates (on non-simulated data, cannot be directly known). rc_meas
is a list consisting of four length N vectors: phi_m, the fracture modern,
sig_m, the uncertainty of the fraction modern, trc_m, trc_m, the measurement
in uncalibrated radiocarbon years, and sig_trc_m, the error for trc_m. The
following hierarchical summary may be easier to digest:

sim_spec            A full specification of the simulation
model_spec          A specification for the model used to generate the
data
density_type    The model type to use for the target density (e.g., a
a truncated Gaussian mixture)
th              The parameter vector for the target density
is_AD           A boolean indicating whether the parameter vector
assumes AD or BP (AD is 1950 - BP, which allows
negative dates)
error_spec      A specification for the measurement errors
N               The number of random samples to make
calib_curve     The calibration curve (currently, must be a named intcal
curve, but support for arbitrary curves could be added
in the future)
seed            An optional random number seed to use to ensure
reproducibility

data                The simulated data
dates               A vector of original dates for the samples
rc_meas             A list with the radiocarbon measurements
phi_m               A vector of simulated fraction moderns
sig_m               A vector of simulated standard deviations for
phi_m
trc_m               A vector of simulated uncalibrated radiocarbon
years (BP)
sig_trc_m           A vector of simulated uncertainties for trc_m
}
