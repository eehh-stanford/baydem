#' @title
#' Add a shaded region showing the min/max quantiles
#'
#' @description
#' Add a shaded region showing the min/max quantiles (usually +/-2.5%) to a
#' plot. The input, bayesian_summary, is the result of a call to
#' \code{summarize_bayesian_inference}, which is a list-like object of class
#' bd_bayesian_summary with information on the quantiles of the density function
#' and growth rate. Typically, the min/max quantiles are 2.5% and 97.5%, but
#' this can be changed via the input lev in the call to
#' \code{summarize_bayesian_inference}.
#'
#' @param bayesian_summary A list-like object of class
#'   \code{bd_bayesian_summary} with information on the quantiles of the density
#'   function and growth rate (see \code{summarize_bayesian_inference}).
#' @param col The color of the shaded region (default:
#'   \code{grDevices::adjustcolor("grey",alpha.f = 0.5)}.
#' @param ... Additional parameters to pass to \code{graphics::polygon}.
#'
#' @seealso [summarize_bayesian_inference()]
#'
#' @export
add_shaded_quantiles <-
  function(bayesian_summary,
           col = grDevices::adjustcolor("grey",
             alpha.f = 0.5
           ),
           ...) {

    # The indidces of the minimum and maximum quantiles
    ind_min <- which.min(bayesian_summary$probs)
    ind_max <- which.max(bayesian_summary$probs)
    graphics::polygon(c(bayesian_summary$tau, rev(bayesian_summary$tau)),
      c(
        bayesian_summary$Qdens[ind_min, ],
        rev(bayesian_summary$Qdens[ind_max, ])
      ),
      border = NA,
      xlab = NULL,
      col = col
    )
  }

#' @title
#' Prepare a blank plot for showing densities
#'
#' @description
#'
#' The input, bayesian_summary, is the result of a call to
#' \code{summarize_bayesian_inference}. Prepare a blank plot for adding density
#' information, e.g. by calling \code{plot_50_percent_quantile} and
#' \code{add_shaded_quantiles}. If the axis labels and limits are not specified,
#' sensible defaults are used.
#'
#' @param bayesian_summary A list-like object of class
#'   \code{bd_bayesian_summary} with information on the quantiles of the density
#'   function and growth rate (see \code{summarize_bayesian_inference}).
#' @param xlim Numeric vector of length 2, giving the x coordinate range.
#' @param ylim Numeric vector of length 2, giving the y coordinate range.
#' @param xlab x-axis label.
#' @param ylab y-axis label.
#' @param ... Additional parameters to pass to plot.
#'
#' @seealso [summarize_bayesian_inference()]
#'
#' @export
make_blank_density_plot <-
  function(bayesian_summary,
           xlim = NA,
           ylim = NA,
           xlab = "Calendar Date [AD]",
           ylab = "Density",
           ...) {
    if (all(is.na(xlim))) {
      xlim <- range(bayesian_summary$tau)
    }

    if (all(is.na(ylim))) {
      # The following command works even if, e.g., f_sim is not in the list an
      # because an$f_sim results in NULL
      ylim <- c(0,
                max(c(max(bayesian_summary$Qdens),
                      bayesian_summary$f_spdf, bayesian_summary$f_sim)))
    }

    graphics::plot(NULL,
      xlim = xlim,
      ylim = ylim,
      ylab = ylab,
      xlab = xlab,
      ...
    )
  }

#' @title
#' Plot the 50 percent quantile curve
#'
#' @description
#' The input, bayesian_summary, is the result of a call to
#' \code{summarize_bayesian_inference}. By default, a new plot is made, but
#' if \code{add = TRUE} the curve is added to the active plot.
#'
#' @param bayesian_summary A list-like object of class
#'   \code{bd_bayesian_summary} with information on the quantiles of the density
#'   function and growth rate (see \code{summarize_bayesian_inference}).
#' @param add Whether to make a new plot or add to the active plot (default:
#'   FALSE).
#' @param ... Additional parameters to pass to plot / lines
#'
#' @seealso [summarize_bayesian_inference()]
#'
#' @export
plot_50_percent_quantile <- function(bayesian_summary, add = F, ...) {

  # The index in probs / Qdens of the 50 percent quantile
  ind50 <- which(bayesian_summary$probs == 0.5)
  if (add) {
    graphics::lines(bayesian_summary$tau,
                    bayesian_summary$Qdens[ind50, ], ...)
  } else {
    graphics::plot(bayesian_summary$tau,
                   bayesian_summary$Qdens[ind50, ], type = "l", ...)
  }
}

#' @title
#' Plot the known simulation density
#'
#' @description
#' The input, bayesian_summary, is the result of a call to
#' \code{summarize_bayesian_inference}. By default, a new plot is made, but
#' if \code{add = TRUE} the curve is added to the active plot. By default, a
#' new plot is made, but if add = TRUE the curve is added to the active plot.
#'
#' @param bayesian_summary A list-like object of class
#'   \code{bd_bayesian_summary} with information on the quantiles of the density
#'   function and growth rate (see \code{summarize_bayesian_inference}).
#' @param add Whether to make a new plot or add to the active plot (default:
#'   FALSE).
#' @param ... Additional parameters to pass to plot / lines
#'
#' @seealso [summarize_bayesian_inference()]
#'
#' @export
plot_known_sim_density <- function(bayesian_summary, add = F, ...) {
  if (add) {
    graphics::lines(bayesian_summary$tau,
                    bayesian_summary$f_sim, ...)
  } else {
    graphics::plot(bayesian_summary$tau,
                   bayesian_summary$f_sim, type = "l", ...)
  }
}

#' @title Plot the summed probability density function (SPDF)
#'
#' @description
#' The input, bayesian_summary, is the result of a call to
#' \code{summarize_bayesian_inference}. By default, a new plot is made, but
#' if \code{add = TRUE} the curve is added to the active plot. By default, a
#' new plot is made, but if add = TRUE the curve is added to the active plot.
#'
#' @param bayesian_summary A list-like object of class
#'   \code{bd_bayesian_summary} with information on the quantiles of the density
#'   function and growth rate (see \code{summarize_bayesian_inference}).
#' @param add Whether to make a new plot or add to the active plot (default:
#'   FALSE).
#' @param ... Additional parameters to pass to plot / lines
#'
#' @seealso [summarize_bayesian_inference()]
#'
#' @export
plot_summed_density <- function(bayesian_summary, add = F, ...) {
  if (add) {
    graphics::lines(bayesian_summary$tau,
                    bayesian_summary$f_spdf, ...)
  } else {
    graphics::plot(bayesian_summary$tau,
                   bayesian_summary$f_spdf, type = "l", ...)
  }
}

#' @title
#' Visualize the calibration curve with equifinal and non-equifinal time spans
#'
#' @description
#' Vizualize the input calibration curve, calib_df, on the interval tau_min to
#' tau_max. This involves plotting the curve itself and shading invertible
#' (non-equifinal) and non-invertible (equifinal) regions.
#'
#' @param tau_min The minimum calendar date for plotting (AD)
#' @param tau_max The maximum calendar date for plotting (AD)
#' @param calib_df The calibration data frame, with columns year_BP,
#'   uncal_year_BP, and uncal_year_BP_error
#' @param invert_col The color for shading invertible regions
#'   (default: `gray90`)
#' @param point_col The color for calibration curve points (default: `black`)
#' @param point_pch The symbol for calibration curve points (default: `19`)
#' @param ... Additional inputs to plots
#'
#' @seealso [load_calib_curve()] for the format of \code{calib_df}
#'
#' @export
vis_calib_curve <-
  function(tau_min,
           tau_max,
           calib_df,
           invert_col = "gray90",
           point_col = "black",
           point_pch = 19,
           ...) {
    tau_curve <- 1950 - calib_df$year_BP
    phi_curve <- exp(-calib_df$uncal_year_BP / 8033)
    ind <- (tau_curve >= tau_min) & (tau_curve <= tau_max)
    tau_vect <- tau_curve[ind]
    phi_vect <- phi_curve[ind]
    phi_min <- min(phi_vect)
    phi_max <- max(phi_vect)

    # Create an empty plot
    graphics::plot(1, type = "n", xlim = c(tau_min, tau_max),
                   ylim = c(phi_min, phi_max), ...)

    equi_info <- assess_calib_curve_equif(calib_df)
    inv_span_list <- equi_info$inv_span_list
    for (ii in 1:length(inv_span_list)) {
      inv_span <- inv_span_list[[ii]]
      if (dplyr::between(inv_span$tau_left, tau_min, tau_max) ||
        dplyr::between(inv_span$tau_right, tau_min, tau_max)) {
        graphics::rect(inv_span$tau_left, phi_min, inv_span$tau_right,
                       phi_max, border = NA, col = invert_col)
      }
    }

    # Draw calibration curve points
    graphics::points(tau_vect, phi_vect, col = point_col, pch = point_pch)
  }
