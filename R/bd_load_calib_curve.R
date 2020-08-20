if (getRversion() >= "2.15.1") {
  utils::globalVariables(c(
    "CAL BP",
    "14C age",
    "Error"
  ))
}

#' @title Load Calibration Curve
#'
#' @description Parse and return the radiocarbon calibration curve stored in data
#'
#' @param calibCurve Name of calibration curve
#'
#' @export
#'
#' @importFrom magrittr %>%
#'
#' @return The calibration dataframe, with columns yearBP, uncalYearBP, and uncalYearBPError
#'
#' @author Michael Holton Price <MichaelHoltonPrice@gmail.com>

bd_load_calib_curve <- function(calibCurve) {
  if (!(calibCurve %in% c("intcal20", "marine20", "shcal20", "intcal13", "marine13", "shcal13"))) {
    stop(paste("Unknown calibration curve name:", calibCurve))
  }

  # calibDf <- read.csv(calibFile,comment.char='#',header=F)
  # calibDf <- calibDf[,1:3]
  # colnames(calibDf) <- c('yearBP','uncalYearBP','uncalYearBPError')

  if (calibCurve == "intcal20") {
    calibCurve <- baydem::intcal20
  } else if (calibCurve == "marine20") {
    calibCurve <- baydem::marine20
  } else if (calibCurve == "shcal20") {
    calibCurve <- baydem::shcal20
  } else if (calibCurve == "intcal13") {
    calibCurve <- baydem::intcal13
  } else if (calibCurve == "marine13") {
    calibCurve <- baydem::marine13
  } else if (calibCurve == "shcal13") {
    calibCurve <- baydem::shcal13
  }

  calibDf <- calibCurve %>%
    dplyr::select(
      `CAL BP`,
      `14C age`,
      Error
    ) %>%
    dplyr::rename(
      yearBP = `CAL BP`,
      uncalYearBP = `14C age`,
      uncalYearBPError = `Error`
    )

  return(calibDf)
}
