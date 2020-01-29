#' @title For each sample, calculate the time it takes for the density to decrease by half from the peak
#'
#' @details
#' Calculate the relative density for two dates or, more generally, for two
#'
#' @param soln The result of a call to bd_do_inference
#' @param anal The result of a call to bd_analyze_soln
#' @param propChange (Default 0.5) The proportional
#' @param spec2 The specification for the second density (see details)
#' @param ind (Optional) Indices at which to do the calculation. By default,
#'            all the samples in anal are used.
#' @param anal (Optional) The result of a call to bd_analyze_soln. This is only
#'             needed if either spec1 or spec2 is 'peak'
#'
#' @return A vector of relative densities (f_spec1 / f_spec2)
#'
#' @export
bd_calc_relative_density <- 
  function(soln, 
           spec1, 
           spec2, 
           ind = NA, 
           anal = NA) {
  TH <- bd_extract_param(soln$fit)
  N <- nrow(TH)
  if (all(is.na(ind))) {
    ind <- 1:N
  }


  # Interpret and do error checking on inputs by calling helper function below
  spec1 <- unpack_spec(spec1, soln, T)
  spec2 <- unpack_spec(spec2, soln, F)

  if (spec1$type == "peak" || spec2$type == "peak") {
    if (all(is.na(anal))) {
      anal <- bd_analyze_soln(soln) # If ind is not NA, this may involve unused computation
    }
    summList <- anal$summList[ind]
  }

  # Calculate the density for spec1
  if (spec1$type == "point") {
    f1 <- calc_point_density(TH[ind, ], soln, spec1$value)
  } else if (spec1$type == "range") {
    f1 <- calc_range_density(TH[ind, ], soln, spec1$lower, spec1$upper)
  } else if (spec1$type == "peak") {
    f1 <- calc_peak_density(summList)
  } else {
    # This should not happen, but throw an error regardless
    stop("Unsupported spec type")
  }

  # Calculate the density for spec2
  if (spec2$type == "point") {
    f2 <- calc_point_density(TH[ind, ], soln, spec2$value)
  } else if (spec2$type == "range") {
    f2 <- calc_range_density(TH[ind, ], soln, spec2$lower, spec2$upper)
  } else if (spec2$type == "peak") {
    f2 <- calc_peak_density(summList)
  } else {
    # This should not happen, but throw an error regardless
    stop("Unsupported spec type")
  }

  return(f1 / f2)
}

# A helper function to unpack and do error checking on inputs spec1 / spec2
unpack_spec <- function(spec, soln, isOne) {
  # For more informative error messages, use the input isOne to set the string
  # s to spec1 or spec2
  if (isOne) {
    s <- "spec1"
  } else {
    s <- "spec2"
  }

  # Handle the supported cases, throwing an error if necessary
  if (is.numeric(spec)) {
    if (length(spec) == 1) { # Numeric / length 1
      point <- spec
      if (point < soln$prob$hp$ymin || soln$prob$hp$ymax < point) {
        stop(paste(s, "is a single date, but not in the range ymin to ymax"))
      }
      return(list(type = "point", value = point))
    } else if (length(spec) == 2) { # Numeric / length 2
      lower <- spec[1]
      if (lower < soln$prob$hp$ymin || soln$prob$hp$ymax < lower) {
        stop(paste(s, "is a date range, but lower value is not in the range ymin to ymax"))
      }
      upper <- spec[2]
      if (upper < soln$prob$hp$ymin || soln$prob$hp$ymax < upper) {
        stop(paste(s, "is a date range, but upper value is not in the range ymin to ymax"))
      }
      if (lower > upper) {
        stop(paste(s, "is a date range, but lower value is greater than upper value"))
      }
      return(list(type = "range", lower = lower, upper = upper))
    } else { # Numeirc / not length 1 or 2
      stop(paste(s, "is numeric, but is neither a single date nor a date range"))
    }
  } else if (is.character(spec)) { # Character
    if (spec == "peak") {
      return(list(type = "peak"))
    } else {
      stop(paste(s, "is a character, but not equal to peak"))
    }
  } else { # Neither character nor numeric
    stop(paste(s, "is neither numeric nor a character"))
  }
}

# A helper function to calculate point densities
calc_point_density <- function(TH, soln, y) {
  return(as.numeric(bd_calc_gauss_mix_pdf_mat(TH, y, ymin = soln$prob$hp$ymin, ymax = soln$prob$hp$ymax)))
}

# A helper function to calculate the mean density over a range
calc_range_density <- function(TH, soln, ylo, yhi) {
  flo <- as.numeric(bd_calc_gauss_mix_pdf_mat(TH, ylo, ymin = soln$prob$hp$ymin, ymax = soln$prob$hp$ymax, type = "cumulative"))
  fhi <- as.numeric(bd_calc_gauss_mix_pdf_mat(TH, yhi, ymin = soln$prob$hp$ymin, ymax = soln$prob$hp$ymax, type = "cumulative"))
  return((fhi - flo) / (yhi - ylo))
}

# A helper function to calculate the peak density
calc_peak_density <- function(summList) {
  return(unlist(lapply(summList, function(x) {
    x$fpeak
  })))
}
