#' @title
#' Set the radiocarbon measurements (rc_meas) for an analysis
#'
#' @description
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
#' (1) import_rc_data
#'     Import the radiocarbon measurements from file (or simulate them using
#'     simulate_rc_data)
#'
#' (2) set_rc_meas
#'     Set the radiocarbon measurements (or use set_sim)
#'
#' (3) calc_tau_range
#'     Determine the range of calendar dates that span the radiocarbon
#'     measurements
#'
#' (4) set_density_model
#'     Set the parametric model to be used for Bayesian inference. Currently,
#'     only a truncated Gaussian mixture is supported, which requires setting
#'     the calendar date range for truncation (likely using the result from step
#'     (3)) and choosing the number of mixture components, K. If K is a vector,
#'     Bayesian inference will be done for each element of the vector during
#'     step (6), do_bayesian_inference.
#'
#' (5) set_calib_curve
#'     Set the calibration curve to use for fitting and analysis
#'
#' (6) do_bayesian_inference
#'     Call sample_theta to do the Bayesian inference
#' @param data_dir The directory in which to store analysis data
#' @param analysis_name A unique name for a given analysis in data_dir
#' @param rc_meas The radiocarbon measurements to use for this analysis (the
#'   format of rc_meas is as output by import_rc_data)
#'
#' @seealso [import_rc_data()] for the format of \code{rc_meas}'
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

#' @title
#' Import radiocarbon data from a .csv file
#'
#' @description
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
#' @param file_name The file name of the input .csv file (or full path to the
#'   file if not located in the current directory)
#' @param trc_m_col The column containing the radiocarbon years values
#'   (optional, but see details).
#' @param sig_m_col The column containing the error for the fraction modern
#'   values (optional, but see details)
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

#' @title
#' Use the results of a simulation to set the radiocarbon measurements (rc_meas)
#' for an analysis (and store the settings used to do the simulation)
#'
#' @description
#' This is one of a set of helper functions for undertaking a typical analysis
#' of radiocarbon dates. As the analysis proceeds, results are stored in a save
#' file called analysis_name.rds in the folder data_dir. Where results are
#' non-determistic, random number seeds are set and stored to ensure that, even
#' if processing is interupted, results are fully reproducible.
#'
#' @param data_dir The directory in which to store analysis data
#' @param analysis_name A unique name for a given analysis in data_dir
#' @param sim The simulation object (see simulate_rc_data)
#'
#' @seealso [set_sim()] for an overview of the "standard pipeline"
#'
#' @export
set_sim <- function(data_dir,analysis_name,sim) {
  data_file <- file.path(data_dir,paste0(analysis_name,".rds"))

  if (file.exists(data_file)) {
    stop("A save file for analysis_name already exists in data_dir")
  }

  calib_df <- load_calib_curve(sim$sim_spec$calib_curve)

  # rc_meas are stored twice for convenience since they are also contained in
  # the input simulation object
  analysis <- list(rc_meas=sim$data$rc_meas,calib_df=calib_df,sim=sim)
  saveRDS(analysis,data_file)
}

#' @title
#' Calculate a calendar date range that spans the input radiocarbon measurements
#'
#' @description
#' To determine the calendar dates, first Bchron::BchronCalibrate is called to
#' calibrate the individual dates, which yields a date range for each. The
#' maximum range across dates is identified (BchronCalibrate adopts a spacing of
#' one year), then (if necessary) the range is slighty extended to multiples of
#' dtau. If intcal20 is being used, dtau should likely be 5, since that is the
#' spacing at which the calibration curve is specified. If dtau is not an
#' integer, it is ignored and a warning is thrown. If dtau is 1 (the default),
#' no rounding is done since that is already the spacing returned by Bchron.
#'
#' @param rc_meas The radiocarbon measurements (see import_rc_data for the
#'   expected format)
#' @param calibration_curve The name of the calibration curve to use (default:
#'   intcal20)
#' @param dtau An integer to round to in extending the calendar range on either
#'   end (default: 1, which has no effect)
#'
#' @return A list containing the minimum and maximum calendar dates, tau_min and
#'   tau_max
#'
#' @seealso
#' * [set_sim()] for an overview of the "standard pipeline"
#' * [import_rc_data()] for the format of \code{rc_meas}'
#'
#' @export
calc_tau_range <- function(rc_meas,calibration_curve="intcal20",dtau=1) {
  N <- length(rc_meas$trc_m)
  calibrations <- Bchron::BchronCalibrate(ages=round(rc_meas$trc_m),
                                          ageSds=round(rc_meas$sig_trc_m),
                                          calCurves=rep(calibration_curve,N))

  # The following range across calibrations is in years BP
  total_range <- range(as.vector(unlist(lapply(
    calibrations,function(calib){range(calib$ageGrid)}))))

  tau_min <- 1950-total_range[2]
  tau_max <- 1950-total_range[1]

  if (dtau %% 1 == 0) {
    # The spacing is already 1, so only do the rounding if dtau is not 1
    if(dtau != 1) {
      tau_min <- tau_min - ( tau_min %% dtau)
      tau_max <- tau_max + (-tau_max %% dtau)
    }
  } else {
    warning("dtau is being ignored because it is not an integer")
  }

  return(list(tau_min=tau_min,tau_max=tau_max))
}

#' @title
#' Set the density model (density_model) for an analysis
#'
#' @description
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
#' @seealso [set_sim()] for an overview of the "standard pipeline"
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

#' @title
#' Set the radiocarbon calibration curve for an analysis
#'
#' @description
#' This is one of a set of helper functions for undertaking a typical analysis
#' of radiocarbon dates. For details on the overall framework for these helper
#' function, see set_rc_meas.
#'
#' @param data_dir The directory in which to store analysis data
#' @param analysis_name A unique name for a given analysis in data_dir
#' @param calibration_curve Either the name of a calibration curve or a data
#'   frame specifying the calibration curve (see load_calib_curve format and
#'   choices).
#'
#' @seealso
#' * [set_sim()] for an overview of the "standard pipeline"
#' * [load_calib_curve()] for the format of \code{calib_df}
#'
#' @export
set_calib_curve <- function(data_dir,analysis_name,calibration_curve) {
  data_file <- file.path(data_dir,paste0(analysis_name,".rds"))

  if (!file.exists(data_file)) {
    stop("A save file for analysis_name does not exist in data_dir")
  }

  analysis <- readRDS(data_file)

  if ("calib_df" %in% names(analysis)) {
    stop("A calibration curve has already been defined for this analysis")
  }

  if (is(calibration_curve,"data.frame")) {
    analysis$calib_df <- calibration_curve
  } else {
    # Otherwise, calibration_curve should be the name of a calibration curve
    # that load_calib_curve recognizes.
    analysis$calib_df <- load_calib_curve(calibration_curve)
  }
  saveRDS(analysis,data_file)
}

#' @title
#' Do Bayesian inference given a set of radiocarbon measurements (rc_meas) and
#' a density model (density_model).
#'
#' @description
#' This is one of a set of helper functions for undertaking a typical analysis
#' of radiocarbon dates. For details on the overall framework for these helper
#' function, see set_rc_meas.
#'
#' #' To ensure reproducibility, the input_seed can be provided. The input seed
#' should either be (1) NA [the default], (2) a single integer, or (3) a matrix
#' with dimensions `num_models` by 2, where
#' `num_models = length(analysis$density_model$K)`. If no seed is provided, one
#' is drawn and treated as if it (a single integer) was input. If a single
#' integer is provided, it is used to generate a matrix of integers with
#' dimensions `num_models` by 2. The first column of the matrix provides the
#' seeds for initializing the parameter vector and the second column provides
#' the seeds for stan.
#'
#' The best model is identified based on the widely applicable information
#' criterion (WAIC) and the corresponding index (in density_model$K) and value
#' are stored in, respectively, m_K_best and K_best (currently, only a truncated
#' Gaussian mixture is supported for the density model specification).
#'
#' @param data_dir The directory in which to store analysis data.
#' @param analysis_name A unique name for a given analysis in data_dir.
#' @param hp Hyperparameters for the priors and to specify the spacing of the
#'   Riemann sum that approximates the integral for the likelihood (see
#'   sample_theta).
#' @param input_seed An optional seed that can be used to make results
#'   reproducible. The input_seed must be either (1) NA / not provided (the
#'   default), (2) a single integer, or (3) a matrix with dimensions
#'   `num_models` by 2 (see description for further details).
#' @param Control arguments to pass to the Bayesian inference function (see
#'   sample_theta).
#'
#' @seealso [set_sim()] for an overview of the "standard pipeline"
#'
#' @export
do_bayesian_inference <- function(data_dir,
                                   analysis_name,
                                   hp,
                                   input_seed=NA,
                                   control=list()) {

  data_file <- file.path(data_dir,paste0(analysis_name,".rds"))
  if (!file.exists(data_file)) {
    stop("A save file for analysis_name does not exist in data_dir")
  }

  analysis <- readRDS(data_file)

  if (!("rc_meas" %in% names(analysis))) {
    stop("Radiocarbon measurements have not specified for this analysis")
  }

  if (!("density_model" %in% names(analysis))) {
    stop("A density model has not been specified for this analysis")
  }

  if (!("calib_df" %in% names(analysis))) {
    stop("A calibration curve has not been specified for this analysis")
  }

  if (analysis$density_model$type != "trunc_gauss_mix") {
    stop("Unsupported type for density_model")
  }

  if ("bayesian_solutions" %in% names(analysis)) {
    stop("Bayesian inference has already been done for this analysis")
  }

  # If additional model types are added, the behavior of input_seed may need to
  # be changed (though it may not need to be changed, too).
  num_models <- length(analysis$density_model$K)
  if (is.vector(input_seed)) {
    if (is.na(input_seed)) {
      base_seed <- sample.int(1000000,1)
    } else {
      base_seed <- input_seed
    }
    set.seed(base_seed)
    seed_mat <- matrix(sample.int(1000000,2*num_models),ncol=2)
  } else if(is.matrix(input_seed)) {
      if(dim(input_seed) != c(num_models,2))
      stop(paste0("If input_seed is a matrix, it must have dimensions ",
                  "num_models x 2 (see function details)"))
    base_seed <- NA
    seed_mat <- input_seed
  } else {
    stop(paste0("input_seed must be NA, a single integer, or a matrix (see ",
                "function details)"))
  }

  # Save the random number seed information
  analysis$input_seed <- input_seed
  analysis$base_seed  <- base_seed
  analysis$seed_mat   <- seed_mat

  # analysis$density_model may be a vector. Loop over values of K to do
  # inference
  analysis$bayesian_solutions <- list()
  for (m_K in 1:num_models) {
    K <- analysis$density_model$K[m_K]
    modified_density_model <- analysis$density_model
    modified_density_model$K <- K
    analysis$bayesian_solutions[[m_K]] <-
      sample_theta(analysis$rc_meas,
                   modified_density_model,
                   hp,
                   analysis$calib_df,
                   init_seed=seed_mat[m_K,1],
                   stan_seed=seed_mat[m_K,2],
                   control=control)
  }
  analysis$hp <- hp
  analysis$K_best <- get_best_K(analysis$bayesian_solutions)
  analysis$m_K_best <- which(analysis$density_model$K == analysis$K_best)
  saveRDS(analysis,data_file)
}

#' @title
#' For a mixture model, get the best value of K for an input list of Bayesian
#' solutions created (based on the widely applicable information criterion,
#' WAIC).
#'
#'
#' @description
#' Get the best number of mixture components (K) from the Bayesian inference
#' based on the widely applicable information criterion (WAIC).
#'
#' @param bayesian_solutions The list of Bayesian "solutions" (see
#' do_bayesian_inference).
#'
#' @seealso [do_bayesian_inference()]
#'
#' @returns The best value of K
#'
#' @export
get_best_K <- function(bayesian_solutions) {
  waic_vect <- rep(NA,length(bayesian_solutions))
  for (m_K in 1:length(bayesian_solutions)) {
    log_lik_mat <- rstan::extract(bayesian_solutions[[m_K]]$fit,"logh")[[1]]
    waic_analysis <- loo::waic(log_lik_mat)
    waic_vect[m_K] <- waic_analysis$estimates["waic","Estimate"]
  }
  m_K_best <- which.min(waic_vect)
  TH <- extract_param(bayesian_solutions[[m_K_best]]$fit)
  return(ncol(TH)/3)
}

#' @title
#' Call summarize_bayesian_inference for the best model created by
#' do_bayesian_inference (the previous step in the standard pipeline).
#'
#' Plot the best model among those for which Bayesian inference was done with
#' the standard pipeline. The plot is a standard plot shows the 50% and +/- 2.5%
#' quantiles. If an output file name for is provided (plot_file_name), the plot
#' is written to file. The format is detected from the extension (e.g., if
#' plot_file_name is "analysis.pdf", a PDF file is written. If the extension is
#' not one of (pdf, png, or jpg, an error is thrown).
#'
#' @param data_dir The directory in which to store analysis data.
#' @param analysis_name A unique name for a given analysis in data_dir.
#' @param plot_file_name An optional file name for saving the file
#'
#' @seealso [do_bayesian_inference()]
#'
#' @export
do_bayesian_summary <- function(data_dir,
                               analysis_name) {
  # TODO: consider writing a stand-alone helper function that does error
  #  checking for all the functions in the standard pipeline.

  data_file <- file.path(data_dir,paste0(analysis_name,".rds"))
  if (!file.exists(data_file)) {
    stop("A save file for analysis_name does not exist in data_dir")
  }

  analysis <- readRDS(data_file)

  if (!("rc_meas" %in% names(analysis))) {
    stop("Radiocarbon measurements have not specified for this analysis")
  }

  if (!("density_model" %in% names(analysis))) {
    stop("A density model has not been specified for this analysis")
  }

  if (!("calib_df" %in% names(analysis))) {
    stop("A calibration curve has not been specified for this analysis")
  }

  if (!("bayesian_solutions" %in% names(analysis))) {
    stop("Bayesian inference has not been done for this analysis")
  }

  if ("bayesian_summary" %in% names(analysis)) {
    stop("A summary has already been done for this analysis")
  }

  modified_density_model <- analysis$density_model
  modified_density_model$K <- analysis$K_best

  analysis$bayesian_summary <- summarize_bayesian_inference (
     analysis$bayesian_solutions[[analysis$m_K_best]],
     analysis$rc_meas,
     modified_density_model,
     analysis$calib_df,
     analysis$hp$dtau)

  saveRDS(analysis,data_file)
}

#' @title
#' Plot the best model created by do_bayesian_inference
#' Plot the best model among those for which Bayesian inference was done with
#' the standard pipeline. The plot is a standard plot shows the 50% and +/- 2.5%
#' quantiles. If an output file name for is provided (plot_file_name), the plot
#' is written to file. The format is detected from the extension (e.g., if
#' plot_file_name is "analysis.pdf", a PDF file is written. If the extension is
#' not one of (pdf, png, or jpg, an error is thrown).
#'
#' @param data_dir The directory in which to store analysis data.
#' @param analysis_name A unique name for a given analysis in data_dir.
#' @param plot_file_name An optional file name for saving the file
#'
#' @seealso [do_bayesian_inference()]
#'
#' @export
plot_best_solution <- function(data_dir,
                               analysis_name,
                               plot_file_name=NA) {

  data_file <- file.path(data_dir,paste0(analysis_name,".rds"))
  if (!file.exists(data_file)) {
    stop("A save file for analysis_name does not exist in data_dir")
  }

  analysis <- readRDS(data_file)

  if (!("rc_meas" %in% names(analysis))) {
    stop("Radiocarbon measurements have not specified for this analysis")
  }

  if (!("density_model" %in% names(analysis))) {
    stop("A density model has not been specified for this analysis")
  }

  if (!("calib_df" %in% names(analysis))) {
    stop("A calibration curve has not been specified for this analysis")
  }

  if ("bayesian_solutions" %in% names(analysis)) {
    stop("Bayesian inference has not been done for this analysis")
  }


#  # Make a blank plot
#  bd_make_blank_density_plot(sim_anal[[r]],
#    ylim = c(0, 0.01),
#    xlab = "",
#    ylab = "Density",
#    xaxt = "n",
#    yaxt = "n"
#  )
#
#  # Add the shaded quantiles
#  bd_add_shaded_quantiles(sim_anal[[r]],
#    col = "gray80"
#  )
#
#  # Add the summed probability density
#  bd_plot_summed_density(sim_anal[[r]],
#    lwd = 2,
#    add = T,
#    col = "black"
#  )
#
#  # Add the Bchron fit
#  lines(tau,sim_bc[[r]]$f_bc,
#    lwd = 2,
#    col = "black",
#    lty=3
#  )
#
#  # Add solid 50% quantile
#  bd_plot_50_percent_quantile(sim_anal[[r]],
#    lwd = 2,
#    add = T,
#    col = "red"
#  )
#
#  # Plot the known, target distribution
#  bd_plot_known_sim_density(sim_anal[[r]],
#    lwd = 2,
#    add = T,
#    col = "blue"
#  )



  K <- get_best_K(analysis$bayesian_solutions)

  return(ncol(TH)/3)
}
