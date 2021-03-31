#' Import radiocarbon data from a .csv file
#'
#' Use read.csv to load the file specified by the input file_name (file_name can
#' also be the full path to the file). The inputs phi_m_col, sig_m_col,
#' trc_m_col, and sig_trc_m_col specify which columns contain the radiocarbon
#' data. All are optional, but only three patterns of inputs are valid:
#'
#' (1) If none of the column specifiers are given, it is assumed that the first
#'     column in the .csv file contains trc_m and the second column contains
#'     sig_trc_m.
#' (2) If trc_m_col and sig_trc_m_col are given, phi_m_col and sig_m_col
#'     must not be given.
#' (3) If phi_m_col and sig_m_col are given, trc_m_col and sig_trc_m_col must
#'     not be given.
#'
#' @param sig_m_col The column containing the error for the fraction modern
#'   values (optional, but see details)
#' @param trc_m_col The column containing the radiocarbon years values
#'   (optional, but see details).
#' @param file_name The file name of the input .csv file (or full path to the
#'   file if not located in the current directory)
#' @param phi_m_col The column containing the fraction modern values (optional,
#'   but see details).
#' @param sig_m_col The column containing the error for the radiocarbon year
#'   values (optional, but see details)
#' @param ... Additional arguments to pass to read.csv (e.g., skip=2 would skip
#'   the first two lines when reading the .csv file)
#'
#' @return A list with four named entries, each a vector equal to the number of
#'  samples: phi_m, the fraction modern values; sig_m, the error for phi_m;
#'   trc_m, a vector of uncalibrated radiocarbon years; and sig_trc_m, the error
#'   for trc_m.
#'
#' @export
import_rc_data <- function(file_name,
                           trc_m_col=NA,
                           sig_trc_m_col=NA,
                           phi_m_col=NA,
                           sig_m_col=NA,
                           ...) {

  pattern <- !is.na(c(trc_m_col,sig_trc_m_col,phi_m_col,sig_m_col))

  if (all(pattern == c(F,F,F,F))) {
    # None of the columns are specified. Assume that the first column is trc_m
    # and the second column is sig_trc_m
    have_rc_years <- T
    trc_m_col     <- 1
    sig_trc_m_col <- 2
  } else if (all(pattern ==  c(T,T,F,F))) {
    # Radiocarbon years are input
    have_rc_years <- T
  } else if (all(pattern ==  c(F,F,T,T))) {
    # Fraction moderns are are input
    have_rc_years <- F
  } else {
    stop("Unsupported input pattern for specifying data columns")
  }

  data_frame <- read.csv(file_name,...)
  if(have_rc_years) {
    trc_m     <- data_frame[,trc_m_col]
    sig_trc_m <- data_frame[,sig_trc_m_col]
    phi_m     <- exp(-trc_m/8033)
    sig_m     <- phi_m * sig_trc_m / 8033
  } else {
    phi_m     <- data_frame[,phi_m_col]
    sig_m     <- data_frame[,sig_m_col]
    trc_m     <- -8033 * log(phi_m)
    sig_trc_m <- 8033 * sig_m / phi_m
  }

  return (list(phi_m=phi_m,sig_m=sig_m,trc_m=trc_m,sig_trc_m=sig_trc_m))
}

#' Set the radiocarbon measurements (rc_meas) for an analysis
#'
#' This is one of a set of helper functions for undertaking a typical analysis
#' of radiocarbon dates. As the analysis proceeds, results are stored in a save
#' file called analysis_name.rds in the folder data_dir. Where results are
#' non-determistic, random number seeds are set and stored to ensure that, even
#' if processing is interupted, results are fully reproducible.
#'
#' Since this is the first helper function that is used in the sequence of
#' analysis steps, a full description of the set of helper functions is
#' described here and the documentation in the other functions points here.
#'
#' The variables data_dir and analysis_name must be specified for each helper
#' function. data_dir is a directory that can store the data for multiple
#' analyses. The analysis_name is a unique specified of a given analysis within
#' data_dir. If a new analysis is initiated by calling set_rc_meas yet a data
#' file already exists in data_dir for that analysis name, an error is thrown.
#'
#' The typical set of steps for an analysis is:
#'
#' (1) set_rc_meas  Set the radiocarbon measurements
#' (2) set_model    Set the parametric model to be used for fitting and Bayesian
#'                  inference
#'
#' @param data_dir The directory in which to store analysis data
#' @param analysis_name A unique name for a given analysis in data_dir
#' @param rc_meas The radiocarbon measurements to use for this analysis (the
#'   format of rc_meas is as output by import_rc_data)
#'
#' @export
set_rc_meas <- function(data_dir,analysis_name,rc_meas) {
  data_file <- file.path(data_dir,paste0(analysis_name,".rds"))

  if (file.exists(data_file)) {
    stop("A save file for analysis_name already exists in data_dir")
  }

  # TODO: Consider doing error checking on rc_meas

  analysis <- list(rc_meas=rc_meas)
  saveRDS(analysis,data_file)
}

#' Set the density model (density_model) for an analysis
#'
#' This is one of a set of helper functions for undertaking a typical analysis
#' of radiocarbon dates. For details on the overall framework for these helper
#' function, see set_rc_meas.
#'
#' set_density_model sets the parameteric model to use for the density function
#' assumed to give rise to the set of radiocarbon dates (rc_meas). The input
#' variable density_model is a list in which the field "type" specifies the
#' type of parametric model and additional named fields are type-specific.
#' Currently, only one type is supported: "trunc_gauss_mix", which stands for
#' truncated Gaussian mixture. For type = "trunc_gauss_mix", the following
#' named fields must also be specified in density_model:
#'
#' tau_min The minimum calendar year for truncation (in AD or, more precisely,
#'         1950 - years BP)
#' tau_max The maximum calendar year for truncation (in AD or, more precisely,
#'         1950 - years BP)
#' K       The number of mixture components. K must be an integer greater than
#'         or equal to 2, and can further be a vector of such integers rather
#'         than a single integer (a vector should be specified if model
#'         selection on the number of mixtures will subsequently be done)
#'         and can be a vector
#'
#' @param data_dir The directory in which to store analysis data
#' @param analysis_name A unique name for a given analysis in data_dir
#' @param density_model A list object specifying the density model
#'
#' @export

set_density_model <- function(data_dir,analysis_name,density_model) {
  data_file <- file.path(data_dir,paste0(analysis_name,".rds"))

  if (!file.exists(data_file)) {
    stop("A save file for analysis_name does not exist in data_dir")
  }

  analysis <- readRDS(data_file)

  if ("density_model" %in% names(analysis)) {
    stop("A density model has already been defined for this analysis")
  }

  if (!("type" %in% names(density_model))) {
    stop("type must be a named field in the list density_model")
  }

  if (density_model$type == "trunc_gauss_mix") {
    # Do some error checking on the inputs
    if (!("tau_min" %in% names(density_model))) {
      stop("tau_min must be a field in density_model for trunc_gauss_mix")
    }

    if (!("tau_max" %in% names(density_model))) {
      stop("tau_max must be a field in density_model for trunc_gauss_mix")
    }

    if (!("K" %in% names(density_model))) {
      stop("K must be a field in density_model for trunc_gauss_mix")
    }

    if (density_model$tau_max <= density_model$tau_min) {
      stop("tau_min must be less than tau_max")
    }

    analysis$density_model <- density_model
    saveRDS(analysis,data_file)
  } else {
    stop("Unsupported type for density_model")
  }
}