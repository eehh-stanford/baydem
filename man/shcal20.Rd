% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bd_data.R
\docType{data}
\name{shcal20}
\alias{shcal20}
\title{The SHCal20 Southern Hemisphere atmospheric radiocarbon calibration curve}
\format{
An object of class \code{tbl_df} (inherits from \code{tbl}, \code{data.frame}) with 9500 rows and 5 columns.
}
\source{
http://intcal.org/curves/shcal20.14c
}
\usage{
shcal20
}
\description{
A \code{tibble::tibble} containing the SHCal20 Southern Hemisphere
atmospheric radiocarbon calibration curve.
}
\details{
\code{shcal20} is a \link[tibble:tbl_df-class]{tibble::tbl_df} object with the following fields:
\itemize{
\item \strong{CAL BP} --- The calibrated year before present
\item \strong{14C age} --- The un-calibrated 14C age of the sample, in years before present
\item \strong{Error} --- The measurement error of the 14C age, in years before present
\item \strong{Delta 14C} --- The normalized delta 14C (Δ14C) measurement for the sample, in per mil (‰)
\item \strong{Sigma} --- The standard deviation (σ) of the Δ14C measurement, in per mil (‰)
}
}
\references{
Hogg, A. G., Heaton, T. J., Hua, Q., Palmer, J. G., Turney, C. S. M.,
Southon, J., Bayliss, A., Blackwell, P. G., Boswijk, G., Bronk Ramsey, C.,
Pearson, C., Petchey, F., Reimer, P., Reimer, R. and Wacker, L. (2020)
“SHCal20 SOUTHERN HEMISPHERE CALIBRATION, 0–55,000 YEARS CAL BP,”
Radiocarbon. Cambridge University Press, pp. 1–20.
DOI: \href{https://doi.org/10.1017/RDC.2020.59}{10.1017/RDC.2020.59}
}
\keyword{datasets}
