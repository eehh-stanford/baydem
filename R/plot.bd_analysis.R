#' @title Plot the result of a call to bayDem_analyzeSoln
#'
#' @description Plot the result of a call to bayDem_analyzeSoln. Three types of plots can be created: density, rate, and shaded density. The default is density.
#'
#' @details an is the result of a call to bayDem_analyzeSoln.
#' It is a list-like object of class bayDem_analysis with information
#' on the quantiles of the density function and growth rate.
#' The optional variables plotType, showSPDF and showSim control the type of
#' plot that is made and characteristics of that plot.
#'
#' The default plotType is 'density'.
#' This plots the density the 50% density quantile (solid) and
#' the bounding quantiles given the level variable, lev,
#' using in bayDem_analyzeSoln (dashed).
#' Optionally, the summed probability density and known target distribution
#' from a simulation can be added.
#'
#' plotType 'rate' yields a plot of the rate (deriveative of density divided by density)
#' instead of density. This plots the 50% rate quantile (solid) and the
#' bounding quantiles given the level variable, lev, used in bayDem_analyzeSoln (dashed).
#' Optionally, the rate for the known target distribution from a simulation can be added.
#' Shaded rectangles are added to the plot using an$growthState with the coloring positive / green,
#' zero / grey, negative / red, and missing / white.
#'
#' Optionally, the known target distribution from a simulation can be added.
#'
#' @param an A list-like object of class \code{baydem::bd_analysis} with information on the quantiles of the density function and growth rate
#' @param plotType (default: 'density') The type of plot to make.
#' Options are: 'density' and 'rate'.
#' @param showSPDF (default: FALSE) Whether to plot the summed probability.
#' Ignored unless plotType is 'density'.
#' @param showSim (default: FALSE) Whether to add the target density or rate to the plot.
#' @param ... Additional parameters passed on to other functions.
#'
#' @export
plot.bd_analysis <-
  function(
           an,
           plotType = "density",
           showSPDF = F,
           showSim = F,
           ...) {

    # Is the true, target density available because this was a simulation?
    haveSim <- "f_sim" %in% names(an)
    if (showSim && !haveSim) {
      stop("showSim is TRUE but simulation target distribution is not available")
    }

    # Plot the density (rate=F; default) or the rate (rate=T)?
    if (tolower(plotType) == "density") {
      # Plot the density

      # The maximum value for the y-axis
      plotMaxF <- max(an$Qdens)

      if (showSPDF) {
        plotMaxF <- max(plotMaxF, an$f_spdf)
      }

      if (showSim) {
        plotMaxF <- max(plotMaxF, an$f_sim)
      }

      # Make an empty plot
      graphics::plot(1,
        xlim = range(an$tau),
        ylim = c(0, plotMaxF),
        xlab = "Calendar Date [AD]",
        ylab = "Density"
      )
      # Plot the summed probability density
      if (showSPDF) {
        graphics::lines(an$tau, an$f_spdf, col = "black", lwd = 3, ...)
      }

      # Plot the posterior quantiles
      indMin <- which.min(an$probs)
      ind50 <- which(an$probs == 0.5) # The index in probs / Qdens of the 50% quantile
      indMax <- which.max(an$probs)
      polygon(c(an$tau, rev(an$tau)), c(an$Qdens[indMin, ], rev(an$Qdens[indMax, ])), col = adjustcolor("grey", alpha.f = 0.5), border = NA, xlab = NULL)
      lines(an$tau, an$Qdens[ind50, ], lwd = 3)

      # If necessary, plot the true target density
      if (showSim) {
        graphics::lines(an$tau, an$f_sim, col = "blue", lwd = 3, lty = 1)
      }
    } else if (tolower(plotType) == "rate") {
      # Plot the rate
      # If showSPDF is true it is ignored for this plot type

      # The maximum value for the y-axis
      plotMinR <- min(an$Qrate)
      plotMaxR <- max(an$Qrate)

      if (showSim) {
        plotMinR <- min(plotMinR, an$rate_sim)
        plotMaxR <- max(plotMaxR, an$rate_sim)
      }

      # Make an empty plot
      graphics::plot(1,
        xlim = range(an$tau),
        ylim = c(plotMinR, plotMaxR),
        xlab = "Calendar Date [AD]",
        ylab = "Growth Rate"
      )

      # Add growth band rectangles
      # The rle command gives a run length encoding,
      # consisting of lengths and values of repeated elements of growthState
      bands <- rle(an$growthState)
      firstInd <- 1
      for (bb in 1:length(bands$values)) {
        if (bands$values[bb] == "missing") {
          bandCol <- "white"
        } else if (bands$values[bb] == "negative") {
          bandCol <- "red"
        } else if (bands$values[bb] == "zero") {
          bandCol <- "grey"
        } else if (bands$values[bb] == "positive") {
          bandCol <- "green"
        } else {
          stop("Unrecognized growthState")
        }
        lastInd <- firstInd + bands$lengths[bb] - 1
        tau0 <- an$tau[firstInd]
        tau1 <- an$tau[lastInd] + an$dtau
        r0 <- plotMinR
        r1 <- plotMaxR
        graphics::polygon(c(tau0, tau1, tau1, tau0), c(r0, r0, r1, r1),
          col = grDevices::adjustcolor(bandCol, alpha.f = .5),
          border = NA
        )
        firstInd <- lastInd + 1
      }

      # Plot the posterior quantiles
      ind50 <- which(an$probs == 0.5) # The index in probs / Qdens of the 50% quantile
      for (ii in 1:nrow(an$Qrate)) {
        lineCol <- "black"
        if (ii == ind50) {
          lineType <- 1 # Line type is 1 (solid) for the 50% quantile
        } else {
          lineType <- 2 # Line type is 2 (dashed) for other quantiles
        }
        graphics::lines(an$tau[an$rateInd], an$Qrate[ii, ], col = lineCol, lwd = 3, lty = lineType)
      }

      # If necessary, plot the true target density
      if (showSim) {
        graphics::lines(an$tau[an$rateInd], an$rate_sim[an$rateInd], col = "blue", lwd = 3, lty = 1)
      }
    } else {
      stop("Unrecognized plot type")
    }
  }
