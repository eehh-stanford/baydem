#' @title Calculate equifinal dates for each point in the calibration curve
#'
#' @details
#' The input calibration data frame has three columns: yearBP, uncalYearBP, and
#' uncalYearBPError. For each date in the column yearBP, determine all other
#' dates with the same fraction modern. These dates are equifinal since, in the
#' absence of additional information, there is no way to determine which year
#' a sample came from. Typically the actual equifinal date lies between two
#' observations in yearBP, so linear interpolation is used to estimate the
#' decimal year. A list of length nrow(calibDf) is returned with, for each
#' point in the calibration curve, the following information:
#'
#' indBase -- The row index in calibDf of the point
#' yBase   -- The calendar date (AD) of the point
#' indEqui -- The row index/indices in calibDf of the equifinal points
#' yEqui   -- The calendar date(s) (AD) of the equifinal points
#'
#' @param calibDf The calibration data frame, with columns yearBP, uncalYearBP, and uncalYearBPError
#'
#' @return A list of equifinality information for each point in the calibration curve (see details)
#'
#' @export
bd_calc_calib_curve_equif_dates <- function(calibDf) {
    # For all the calendar dates in the calibration dataframe, identify points
    # (other calendar dates) with the same fraction modern.
    y_curve <- 1950 - calibDf$yearBP
    phi_curve <- bd_calc_calib_curve_frac_modern(calibDf)
    outputList <- list()
    for(ii in 1:(length(phi_curve)-1)) {
        y <- y_curve[ii]
	phi <- phi_curve[ii]
	# Look for adjacent points that bracket y
	ind1 <- which(phi_curve[1:(length(phi_curve)-1)] >= phi_curve[ii] & phi_curve[2:length(phi_curve)] <= phi_curve[ii])
	ind2 <- which(phi_curve[1:(length(phi_curve)-1)] <= phi_curve[ii] & phi_curve[2:length(phi_curve)] >= phi_curve[ii])
	ind <- setdiff(c(ind1,ind2),ii)
	if(length(ind) > 1) {
            indEqui <- ind
            yEqui <- rep(NA,length(indEqui))
            for(jj in 1:length(indEqui)) {
	        ii_lo <- indEqui[jj]
	        ii_hi <- indEqui[jj]+1
                yEqui[jj] <- y_curve[ii_lo] + (y_curve[ii_hi] - y_curve[ii_lo])*(phi_curve[ii]-phi_curve[ii_lo])/(phi_curve[ii_hi]-phi_curve[ii_lo])
	    }
            newEntry = list(indBase=ii,yBase=y_curve[ii],indEqui=indEqui,yEqui=unique(yEqui))
            outputList[[length(outputList)+1]] <- newEntry
	}
    }
    return(outputList)
}
