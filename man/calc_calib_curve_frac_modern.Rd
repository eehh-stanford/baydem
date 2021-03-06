% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/radiocarbon_functions.R
\name{calc_calib_curve_frac_modern}
\alias{calc_calib_curve_frac_modern}
\title{Calculation the calibration curve fraction modern value for the input
calendar dates.}
\usage{
calc_calib_curve_frac_modern(calib_df, tau = NA, is_BP = FALSE)
}
\arguments{
\item{calib_df}{The calibration data frame, with columns year_BP,
uncal_year_BP, and uncal_year_BP_error.}

\item{tau}{The calendar dates tau (if not input, calib_df$uncal_year_BP is
used).}

\item{is_BP}{Whether the input dates are before present (BP), as opposed to
AD (default: \code{FALSE}).}
}
\value{
The vector of calibration curve fraction modern values
}
\description{
The input calibration data frame has three columns: year_BP, uncal_year_BP,
and uncal_year_BP_error. For each input calendar date, tau, use the
calibration curve information to estimate the fraction modern of the
calibration curve. If tau is not specified, calib_df$uncal_year_BP is used as
the dates for the calculation. By default, dates are assumed to be AD (more
precisely: 1950 - years BP), but this can be changed using the optional input
is_BP, which is FALSE by default.
}
\seealso{
\code{\link[=load_calib_curve]{load_calib_curve()}} for the format of \code{calib_df}'
}
