#' @title Assess the equifinality or time spans in the input radiocarbon
#' calibration curve
#'
#' @details
#' The input calibration data frame has three columns: year_BP, uncal_year_BP,
#' and uncal_year_BP_error. A time span is equifinal if, for each date in the
#' span, there is at least one other date in the radiocarbon calibration curve
#' with the same fraction modern value. Conversely, a time span is not equifinal
#' if this is not true (all dates in the time span correspond to a unique
#' fraction modern value). Thus, the calibration curve is divided into
#' alternating equifinal and non-equifinal spans. This function identifies these
#' regions for the input calibration data frame, calib_df.
#'
#' Although calib_df uses year BP, all calculations and returned data use AD.
#'
#' A list is returned with two named variables: can_invert and inv_span_list
#' (inv = invert, which is conceptually identical to non-equifinal).
#'
#' can_invert is a boolean that indicates whether each entry (row) in the
#' calibration data frame, calib_df, is inside or outside an invertible region.
#'
#' inv_span_list is a list summarizing information about all the invertible
#' (non-equifinal) time spans. It contains the following named variables:
#'
#' ind       -- Indices (rows) of calib_df following within the invertible time
#'              span
#' tau_left  -- Calendar date (AD) of the left (earlier) boundary of the
#'              invertible time span
#' phi_left  -- Fraction modern value of the left (earlier) boundary of the
#'              invertible time span
#' tau_right -- Calendar date (AD) of the right (later) boundary of the
#'              invertible time span
#' phi_right -- Fraction modern value of the right (later) boundary of the
#'              invertible time span
#' ii_prev   -- The index (row) of the calib_df earlier than the invertible
#'              region with closest fraction modern value
#' ii_next   -- The index (row) of the calib_df later than the invertible
#'              region with closest fraction modern value
#'
#' @param calib_df The calibration data frame, with columns year_BP,
#'   uncal_year_BP, and uncal_year_BP_error
#' @param equi_list (Optional) The result of a call to
#'   calc_calib_curve_equif_dates
#'
#' @return A list of with named variables can_invert and inv_span_list (see
#'   details)
#'
#' @export
#'
assess_calib_curve_equif <- function(calib_df, equi_list = NA) {
  if (all(is.na(equi_list))) {
    equi_list <- calc_calib_curve_equif_dates(calib_df)
  }

  can_invert <- rep(T, length(calib_df$year_BP))
  tau_curve <- 1950 - calib_df$year_BP
  phi_curve <- calc_calib_curve_frac_modern(calib_df)
  for (equi_entry in equi_list) {
    can_invert[equi_entry$ind_base] <- F
  }
  clust_list <- evd::clusters(can_invert, .5)
  inv_span_list <- list()
  for (cc in 1:length(clust_list)) {
    clust <- clust_list[[cc]]
    lo <- as.numeric(names(clust)[1])
    hi <- as.numeric(names(clust)[length(clust)])
    # Find left boundary
    if (lo == 1) {
      tau_left <- tau_curve[1]
      phi_left <- phi_curve[1]
      ii_prev <- 1
    } else {
      ii_prev <- which.min(abs(phi_curve[1:(lo - 1)] - phi_curve[lo]))
      tau_left <- phi2tau(tau_curve, phi_curve, phi_curve[ii_prev], lo - 1, lo)
      phi_left <- phi_curve[ii_prev]
    }

    # Find right boundary
    if (hi == length(tau_curve)) {
      tau_right <- tau_curve[length(tau_curve)]
      phi_right <- phi_curve[length(tau_curve)]
      ii_next <- length(tau_curve)
    } else {
      ii_next <- hi + which.min(abs(phi_curve[(hi + 1):length(phi_curve)] - phi_curve[hi]))
      tau_right <- phi2tau(tau_curve, phi_curve, phi_curve[ii_next], hi, hi + 1)
      phi_right <- phi_curve[ii_next]
    }
    inv_span <- list(ind = lo:hi, tau_left = tau_left, phi_left = phi_left,
                     tau_right = tau_right, phi_right = phi_right,
                     ii_prev = ii_prev, ii_next = ii_next)
    inv_span_list[[cc]] <- inv_span
  }
  return(list(inv_span_list = inv_span_list, can_invert = can_invert))
}

#' @title Calculate equifinal dates for each point in the calibration curve
#'
#' @details
#' The input calibration data frame has three columns: year_BP, uncal_year_BP,
#' and uncal_year_BP_error. For each date in the column year_BP, determine all
#' other dates with the same fraction modern. These dates are equifinal since,
#' in the absence of additional information, there is no way to determine which
#' year a sample came from. Typically the actual equifinal date lies between two
#' observations in year_BP, so linear interpolation is used to estimate the
#' decimal year. A list of length nrow(calib_df) is returned with, for each
#' point in the calibration curve, the following information:
#'
#' ind_base -- The row index in calib_df of the point
#' tau_base -- The calendar date (AD) of the point
#' ind_equi -- The row index/indices in calib-df of the equifinal points
#' tau_equi -- The calendar date(s) (AD) of the equifinal points
#'
#' @param calib_df The calibration data frame, with columns year_BP,
#'   uncal_year_BP, and uncal_year_BP_error
#'
#' @return A list of equifinality information for each point in the calibration
#'   curve (see details)
#'
#' @export
calc_calib_curve_equif_dates <- function(calib_df) {
  # For all the calendar dates in the calibration dataframe, identify points
  # (other calendar dates) with the same fraction modern.
  tau_curve <- 1950 - calib_df$year_BP
  phi_curve <- calc_calib_curve_frac_modern(calib_df)
  output_list <- list()
  for (ii in 1:(length(phi_curve) - 1)) {
    tau <- tau_curve[ii]
    phi <- phi_curve[ii]
    # Look for adjacent points that bracket tau
    ind1 <- which(phi_curve[1:(length(phi_curve) - 1)] >= phi_curve[ii] &
                    phi_curve[2:length(phi_curve)] <= phi_curve[ii])
    ind2 <- which(phi_curve[1:(length(phi_curve) - 1)] <= phi_curve[ii] &
                    phi_curve[2:length(phi_curve)] >= phi_curve[ii])
    ind <- setdiff(c(ind1, ind2), ii)
    if (length(ind) > 1) {
      ind_equi <- ind
      tau_equi <- rep(NA, length(ind_equi))
      for (jj in 1:length(ind_equi)) {
        ii_lo <- ind_equi[jj]
        ii_hi <- ind_equi[jj] + 1
        tau_equi[jj] <- tau_curve[ii_lo] +
          (tau_curve[ii_hi] - tau_curve[ii_lo]) *
            (phi_curve[ii] - phi_curve[ii_lo]) /
             (phi_curve[ii_hi] - phi_curve[ii_lo])
      }
      new_entry <- list(ind_base = ii,
                        tau_base = tau_curve[ii],
                        ind_equi = ind_equi,
                        tau_equi = unique(tau_equi))
      output_list[[length(output_list) + 1]] <- new_entry
    }
  }
  return(output_list)
}
