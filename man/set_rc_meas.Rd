% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_io_functions.R
\name{set_rc_meas}
\alias{set_rc_meas}
\title{Set the radiocarbon measurements (rc_meas) for an analysis}
\usage{
set_rc_meas(data_dir, analysis_name, rc_meas)
}
\arguments{
\item{data_dir}{The directory in which to store analysis data}

\item{analysis_name}{A unique name for a given analysis in data_dir}

\item{rc_meas}{The radiocarbon measurements to use for this analysis (the
format of rc_meas is as output by import_rc_data)}
}
\description{
This is one of a set of helper functions for undertaking a typical analysis
of radiocarbon dates. As the analysis proceeds, results are stored in a save
file called analysis_name.rds in the folder data_dir. Where results are
non-deterministic, random number seeds are set and stored to ensure that,
even if processing is interrupted, results are fully reproducible.

There are, in fact, two distinct pipelines for (a) simulations and (b)
analyses that use existing radiocarbon measurements saved in a .csv file.
Since this is the first helper function used in the standard pipeline for
using an input .csv file, that pipeline is described here. See set_sim for
the pipeline when doing a simulation.

The standard pipeline for analyses that use an input data stored in a .csv
file (see the vignette standard_pipeline for additional details):

(0) import_rc_data
Import data from a .csv file using import_rc_data.

(1) set_rc_meas
Use the data from the previous step to set the radiocarbon measurements
to be used for this analysis

(2) calc_tau_range
Determine the range of calendar dates that span the radiocarbon
measurements (or specify the range directly)

(3) set_density_model
Set the parametric model to be used for Bayesian inference. Currently,
only a truncated Gaussian mixture is supported, which requires setting
the calendar date range for truncation (likely using the result from step
(3)) and choosing the number of mixture components, K. If K is a vector,
Bayesian inference will be done for each element of the vector during
step (5), do_bayesian_inference.

(4) set_calib_curve
Set the calibration curve to be used for subsequent steps.

(5) do_bayesian_inference
Call sample_theta to do the Bayesian inference for each model to check.

(6) do_bayesian_summary
Calculate summary statistics for the best model sampled by calling
do_bayesian_inference (the best model is the one with the best WAIC).

(7) plot_best_solution
Create a plot of the best solution.
}
\seealso{
\code{\link[=import_rc_data]{import_rc_data()}} for the format of \code{rc_meas}'
}
