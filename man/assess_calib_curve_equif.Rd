% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/identifiability_functions.R
\name{assess_calib_curve_equif}
\alias{assess_calib_curve_equif}
\title{Assess the equifinality or time spans in the input radiocarbon
calibration curve}
\usage{
assess_calib_curve_equif(calib_df, equi_list = NA)
}
\arguments{
\item{calib_df}{The calibration data frame, with columns year_BP,
uncal_year_BP, and uncal_year_BP_error}

\item{equi_list}{(Optional) The result of a call to
calc_calib_curve_equif_dates}
}
\value{
A list of with named variables can_invert and inv_span_list (see
details)
}
\description{
Assess the equifinality or time spans in the input radiocarbon
calibration curve
}
\details{
The input calibration data frame has three columns: year_BP, uncal_year_BP,
and uncal_year_BP_error. A time span is equifinal if, for each date in the
span, there is at least one other date in the radiocarbon calibration curve
with the same fraction modern value. Conversely, a time span is not equifinal
if this is not true (all dates in the time span correspond to a unique
fraction modern value). Thus, the calibration curve is divided into
alternating equifinal and non-equifinal spans. This function identifies these
regions for the input calibration data frame, calib_df.

Although calib_df uses year BP, all calculations and returned data use AD.

A list is returned with two named variables: can_invert and inv_span_list
(inv = invert, which is conceptually identical to non-equifinal).

can_invert is a boolean that indicates whether each entry (row) in the
calibration data frame, calib_df, is inside or outside an invertible region.

inv_span_list is a list summarizing information about all the invertible
(non-equifinal) time spans. It contains the following named variables:

ind       -- Indices (rows) of calib_df following within the invertible time
span
tau_left  -- Calendar date (AD) of the left (earlier) boundary of the
invertible time span
phi_left  -- Fraction modern value of the left (earlier) boundary of the
invertible time span
tau_right -- Calendar date (AD) of the right (later) boundary of the
invertible time span
phi_right -- Fraction modern value of the right (later) boundary of the
invertible time span
ii_prev   -- The index (row) of the calib_df earlier than the invertible
region with closest fraction modern value
ii_next   -- The index (row) of the calib_df later than the invertible
region with closest fraction modern value
}
