#' The IntCal13 calibration curve
#'
#' A [tibble][tibble::tibble()] containing the IntCal13 Northern Hemisphere
#' atmospheric radiocarbon calibration curve.
#'
#' `intcal13` is a [tibble::tbl_df-class] object with the following fields:
#'
#' * **CAL BP** --- The calibrated year before present
#' * **14C age** --- The un-calibrated 14C age of the sample, in years before present
#' * **Error** --- The measurement error of the 14C age, in years before present
#' * **Delta 14C** --- The normalized delta 14C (Δ14C) measurement for the sample, in per mil (‰)
#' * **Sigma** --- The standard deviation (σ) of the Δ14C measurement, in per mil (‰)
#'
#' @references
#' Reimer PJ, Bard E, Bayliss A, Beck JW, Blackwell PG, Bronk Ramsey C, Buck CE,
#' Cheng H, Edwards RL, Friedrich M, Grootes PM, Guilderson TP, Haflidason H,
#' Hajdas I, Hatté C, Heaton TJ, Hogg AG, Hughen KA, Kaiser KF, Kromer B,
#' Manning SW, Niu M, Reimer RW, Richards DA, Scott EM, Southon JR, Turney CSM,
#' van der Plicht J. IntCal13 and MARINE13 radiocarbon age calibration curves 0-50000 years calBP.
#' *Radiocarbon* 55(4). DOI: [10.2458/azu_js_rc.55.16947](https://doi.org/10.2458/azu_js_rc.55.16947)
#'
#' @source http://www.radiocarbon.org/IntCal13/
"intcal13"


#' The raw data for constructing the IntCal13 calibration curve
#'
#' A \code{base::list} of \code{tibble::tibble} objects containing the data for calculating the IntCal13 Northern Hemisphere
#' atmospheric radiocarbon calibration curve.
#'
#' These data are as available from the IntCal webservice
#' hosted at Queen's University, Belfast. Details on the structure of the data can be found at
#' [http://intcal.qub.ac.uk/intcal13/query/](http://intcal.qub.ac.uk/intcal13/query/).
#'
#' @references
#' Reimer PJ, Bard E, Bayliss A, Beck JW, Blackwell PG, Bronk Ramsey C, Buck CE,
#' Cheng H, Edwards RL, Friedrich M, Grootes PM, Guilderson TP, Haflidason H,
#' Hajdas I, Hatté C, Heaton TJ, Hogg AG, Hughen KA, Kaiser KF, Kromer B,
#' Manning SW, Niu M, Reimer RW, Richards DA, Scott EM, Southon JR, Turney CSM,
#' van der Plicht J. IntCal13 and MARINE13 radiocarbon age calibration curves 0-50000 years calBP.
#' *Radiocarbon* 55(4). DOI: [10.2458/azu_js_rc.55.16947](https://doi.org/10.2458/azu_js_rc.55.16947)
#'
#' @source http://intcal.qub.ac.uk/intcal13/query/query.php
"intcal13_raw"


#' The Marine13 calibration curve
#'
#' A [tibble][tibble::tibble()] containing the Marine13 marine radiocarbon calibration curve.
#'
#' `marine13` is a [tibble::tbl_df-class] object with the following fields:
#'
#' * **CAL BP** --- The calibrated year before present
#' * **14C age** --- The un-calibrated 14C age of the sample, in years before present
#' * **Error** --- The measurement error of the 14C age, in years before present
#' * **Delta 14C** --- The normalized delta 14C (Δ14C) measurement for the sample, in per mil (‰)
#' * **Sigma** --- The standard deviation (σ) of the Δ14C measurement, in per mil (‰)
#'
#' @references
#' Reimer PJ, Bard E, Bayliss A, Beck JW, Blackwell PG, Bronk Ramsey C, Buck CE,
#' Cheng H, Edwards RL, Friedrich M, Grootes PM, Guilderson TP, Haflidason H,
#' Hajdas I, Hatté C, Heaton TJ, Hogg AG, Hughen KA, Kaiser KF, Kromer B,
#' Manning SW, Niu M, Reimer RW, Richards DA, Scott EM, Southon JR, Turney CSM,
#' van der Plicht J. IntCal13 and MARINE13 radiocarbon age calibration curves 0-50000 years calBP.
#' *Radiocarbon* 55(4). DOI: [10.2458/azu_js_rc.55.16947](https://doi.org/10.2458/azu_js_rc.55.16947)
#'
#' @source http://www.radiocarbon.org/IntCal13/
"marine13"


#' The raw data for constructing the Marine13 calibration curve
#'
#' A \code{base::list} of \code{tibble::tibble} objects containing the data for calculating the
#' Marine13 marine radiocarbon calibration curve.
#'
#' These data are as available from the IntCal webservice
#' hosted at Queen's University, Belfast. Details on the structure of the data can be found at
#' [http://intcal.qub.ac.uk/marine/query/](http://intcal.qub.ac.uk/marine/query/).
#'
#' @references
#' Reimer PJ, Bard E, Bayliss A, Beck JW, Blackwell PG, Bronk Ramsey C, Buck CE,
#' Cheng H, Edwards RL, Friedrich M, Grootes PM, Guilderson TP, Haflidason H,
#' Hajdas I, Hatté C, Heaton TJ, Hogg AG, Hughen KA, Kaiser KF, Kromer B,
#' Manning SW, Niu M, Reimer RW, Richards DA, Scott EM, Southon JR, Turney CSM,
#' van der Plicht J. IntCal13 and MARINE13 radiocarbon age calibration curves 0-50000 years calBP.
#' *Radiocarbon* 55(4). DOI: [10.2458/azu_js_rc.55.16947](https://doi.org/10.2458/azu_js_rc.55.16947)
#'
#' @source http://intcal.qub.ac.uk/marine/query/query.php
"marine13_raw"


#' The SHCal13 Southern Hemisphere atmospheric radiocarbon calibration curve
#'
#' A \code{tibble::tibble} containing the SHCal13 Southern Hemisphere atmospheric radiocarbon calibration curve.
#'
#' `shcal13` is a [tibble::tbl_df-class] object with the following fields:
#'
#' * **CAL BP** --- The calibrated year before present
#' * **14C age** --- The un-calibrated 14C age of the sample, in years before present
#' * **Error** --- The measurement error of the 14C age, in years before present
#' * **Delta 14C** --- The normalized delta 14C (Δ14C) measurement for the sample, in per mil (‰)
#' * **Sigma** --- The standard deviation (σ) of the Δ14C measurement, in per mil (‰)
#'
#' @references
#' Alan G Hogg, Quan Hua, Paul G Blackwell, Caitlin E Buck, Thomas P Guilderson,
#' Timothy J  Heaton, Mu Niu, Jonathan G Palmer, Paula J Reimer, Ron W Reimer,
#' Christian S M Turney, Susan R H Zimmerman
#' *Radiocarbon* 55(4). DOI: [10.2458/azu_js_rc.55.16783](https://doi.org/10.2458/azu_js_rc.55.16783)
#'
#' @source http://www.radiocarbon.org/IntCal13/
"shcal13"


#' The raw data for constructing the SHCal13 Southern Hemisphere atmospheric calibration curve
#'
#' A \code{base::list} of \code{tibble::tibble} objects containing the data for calculating the
#' SHCal13 Southern Hemisphere atmospheric radiocarbon calibration curve.
#'
#' These data are as available from the IntCal webservice
#' hosted at Queen's University, Belfast. Details on the structure of the data can be found at
#' [http://intcal.qub.ac.uk/shcal13/query/](http://intcal.qub.ac.uk/shcal13/query/).
#'
#' @references
#' Alan G Hogg, Quan Hua, Paul G Blackwell, Caitlin E Buck, Thomas P Guilderson,
#' Timothy J  Heaton, Mu Niu, Jonathan G Palmer, Paula J Reimer, Ron W Reimer,
#' Christian S M Turney, Susan R H Zimmerman
#' *Radiocarbon* 55(4). DOI: [10.2458/azu_js_rc.55.16783](https://doi.org/10.2458/azu_js_rc.55.16783)
#'
#' @source http://intcal.qub.ac.uk/shcal13/query/query.php
"shcal13_raw"
