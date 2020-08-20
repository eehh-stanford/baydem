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




#' The IntCal20 calibration curve
#'
#' A [tibble][tibble::tibble()] containing the IntCal20 Northern Hemisphere
#' atmospheric radiocarbon calibration curve.
#'
#' `intcal20` is a [tibble::tbl_df-class] object with the following fields:
#'
#' * **CAL BP** --- The calibrated year before present
#' * **14C age** --- The un-calibrated 14C age of the sample, in years before present
#' * **Error** --- The measurement error of the 14C age, in years before present
#' * **Delta 14C** --- The normalized delta 14C (Δ14C) measurement for the sample, in per mil (‰)
#' * **Sigma** --- The standard deviation (σ) of the Δ14C measurement, in per mil (‰)
#'
#' @references
#' Reimer, P. J., Austin, W. E. N., Bard, E., Bayliss, A., Blackwell, P. G.,
#' Bronk Ramsey, C., Butzin, M., Cheng, H., Edwards, R. L., Friedrich, M.,
#' Grootes, P. M., Guilderson, T. P., Hajdas, I., Heaton, T. J., Hogg, A. G.,
#' Hughen, K. A., Kromer, B., Manning, S. W., Muscheler, R., Palmer, J. G.,
#' Pearson, C., van der Plicht, J., Reimer, R. W., Richards, D. A.,
#' Scott, E. M., Southon, J. R., Turney, C. S. M., Wacker, L., Adolphi, F.,
#' Büntgen, U., Capano, M., Fahrni, S. M., Fogtmann-Schulz, A., Friedrich, R.,
#' Köhler, P., Kudsk, S., Miyake, F., Olsen, J., Reinig, F., Sakamoto, M.,
#' Sookdeo, A. and Talamo, S. (2020) “THE INTCAL20 NORTHERN HEMISPHERE
#' RADIOCARBON AGE CALIBRATION CURVE (0–55 CAL kBP),” Radiocarbon.
#' Cambridge University Press, pp. 1–33.
#' DOI: [10.1017/RDC.2020.41](https://doi.org/10.1017/RDC.2020.41)
#'
#' @source http://intcal.org/curves/intcal20.14c
"intcal20"


#' The raw data for constructing the IntCal20 calibration curve
#'
#' A \code{base::list} of \code{tibble::tibble} objects containing the data for
#' calculating the IntCal20 Northern Hemisphere
#' atmospheric radiocarbon calibration curve.
#'
#' These data are as available from the IntCal webservice
#' hosted at Queen's University, Belfast.
#' Details on the structure of the data can be found at
#' [http://intcal.org/JS/JSintcal20/query/](http://intcal.org/JS/JSintcal20/query/).
#'
#' @references
#' Reimer, P. J., Austin, W. E. N., Bard, E., Bayliss, A., Blackwell, P. G.,
#' Bronk Ramsey, C., Butzin, M., Cheng, H., Edwards, R. L., Friedrich, M.,
#' Grootes, P. M., Guilderson, T. P., Hajdas, I., Heaton, T. J., Hogg, A. G.,
#' Hughen, K. A., Kromer, B., Manning, S. W., Muscheler, R., Palmer, J. G.,
#' Pearson, C., van der Plicht, J., Reimer, R. W., Richards, D. A.,
#' Scott, E. M., Southon, J. R., Turney, C. S. M., Wacker, L., Adolphi, F.,
#' Büntgen, U., Capano, M., Fahrni, S. M., Fogtmann-Schulz, A., Friedrich, R.,
#' Köhler, P., Kudsk, S., Miyake, F., Olsen, J., Reinig, F., Sakamoto, M.,
#' Sookdeo, A. and Talamo, S. (2020) “THE INTCAL20 NORTHERN HEMISPHERE
#' RADIOCARBON AGE CALIBRATION CURVE (0–55 CAL kBP),” Radiocarbon.
#' Cambridge University Press, pp. 1–33.
#' DOI: [10.1017/RDC.2020.41](https://doi.org/10.1017/RDC.2020.41)
#'
#' @source http://intcal.org/JS/JSintcal20/query/query.php
"intcal20_raw"


#' The Marine20 calibration curve
#'
#' A [tibble][tibble::tibble()] containing the Marine20 marine radiocarbon calibration curve.
#'
#' `marine20` is a [tibble::tbl_df-class] object with the following fields:
#'
#' * **CAL BP** --- The calibrated year before present
#' * **14C age** --- The un-calibrated 14C age of the sample, in years before present
#' * **Error** --- The measurement error of the 14C age, in years before present
#' * **Delta 14C** --- The normalized delta 14C (Δ14C) measurement for the sample, in per mil (‰)
#' * **Sigma** --- The standard deviation (σ) of the Δ14C measurement, in per mil (‰)
#'
#' @references
#' Heaton, T. J., Köhler, P., Butzin, M., Bard, E., Reimer, R. W.,
#' Austin, W. E. N., Bronk Ramsey, C., Grootes, P. M., Hughen, K. A.,
#' Kromer, B., Reimer, P. J., Adkins, J., Burke, A., Cook, M. S.,
#' Olsen, J. and Skinner, L. C. (2020) “MARINE20—THE MARINE RADIOCARBON AGE
#' CALIBRATION CURVE (0–55,000 CAL BP),” Radiocarbon. Cambridge University
#' Press, pp. 1–42. DOI: [10.1017/RDC.2020.68](https://doi.org/10.1017/RDC.2020.68)
#'
#' @source http://intcal.org/curves/marine20.14c
"marine20"


#' The SHCal20 Southern Hemisphere atmospheric radiocarbon calibration curve
#'
#' A \code{tibble::tibble} containing the SHCal20 Southern Hemisphere
#' atmospheric radiocarbon calibration curve.
#'
#' `shcal20` is a [tibble::tbl_df-class] object with the following fields:
#'
#' * **CAL BP** --- The calibrated year before present
#' * **14C age** --- The un-calibrated 14C age of the sample, in years before present
#' * **Error** --- The measurement error of the 14C age, in years before present
#' * **Delta 14C** --- The normalized delta 14C (Δ14C) measurement for the sample, in per mil (‰)
#' * **Sigma** --- The standard deviation (σ) of the Δ14C measurement, in per mil (‰)
#'
#' @references
#' Hogg, A. G., Heaton, T. J., Hua, Q., Palmer, J. G., Turney, C. S. M.,
#' Southon, J., Bayliss, A., Blackwell, P. G., Boswijk, G., Bronk Ramsey, C.,
#' Pearson, C., Petchey, F., Reimer, P., Reimer, R. and Wacker, L. (2020)
#' “SHCal20 SOUTHERN HEMISPHERE CALIBRATION, 0–55,000 YEARS CAL BP,”
#' Radiocarbon. Cambridge University Press, pp. 1–20.
#' DOI: [10.1017/RDC.2020.59](https://doi.org/10.1017/RDC.2020.59)
#'
#' @source http://intcal.org/curves/shcal20.14c
"shcal20"


#' The raw data for constructing the SHCal20 Southern Hemisphere atmospheric calibration curve
#'
#' A \code{base::list} of \code{tibble::tibble} objects containing the data for calculating the
#' SHCal20 Southern Hemisphere atmospheric radiocarbon calibration curve.
#'
#' These data are as available from the IntCal webservice
#' hosted at Queen's University, Belfast. Details on the structure of the data can be found at
#' [http://intcal.org/JS/JSshcal20/query/](http://intcal.org/JS/JSshcal20/query/).
#'
#' @references
#' Hogg, A. G., Heaton, T. J., Hua, Q., Palmer, J. G., Turney, C. S. M.,
#' Southon, J., Bayliss, A., Blackwell, P. G., Boswijk, G., Bronk Ramsey, C.,
#' Pearson, C., Petchey, F., Reimer, P., Reimer, R. and Wacker, L. (2020)
#' “SHCal20 SOUTHERN HEMISPHERE CALIBRATION, 0–55,000 YEARS CAL BP,”
#' Radiocarbon. Cambridge University Press, pp. 1–20.
#' DOI: [10.1017/RDC.2020.59](https://doi.org/10.1017/RDC.2020.59)
#'
#' @source http://intcal.org/JS/JSshcal20/query/query.php
"shcal20_raw"
