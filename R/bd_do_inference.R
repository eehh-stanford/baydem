#' @title Do Demographic Bayesian Inference
#'
#' This is the core function that implements the Bayesian inference.
#'
#' The input is a problem statement object (a list),
#' prob, that consists of the input data (the vectors phi_m and sig_m)
#' and the hyperparameters (hp).
#' stan is called via the rstan package to sample from the posterior.
#' The output is the variable soln of class bayDem_soln,
#' which is a list with the fields prob (the input) and fit (the result of the stan fit).
#' prob can also have an optional field control that specifies the
#' following control parameters for the Bayesian inference (default in parentheses):
#'   numChains     -- (4)    Number of chains
#'   sampsPerChain -- (2000) Number of samples per chain
#'   initList      --        The initializations for each chain. The default is
#'                           to sample from prior using hyperparameters
#'   adaptDelta    -- (.99)  The value of the adapt_delta parameter, an input to
#'                           stan. The stan default is 0.8, but a higher value
#'                           is sensible for a mixture model with more than a
#'                           few parameters, hence use 0.99 by default
#'
#' @param prob List with the fields `phi_m` (vector of radiocarbon measurements as fraction modern),
#' `sig_m` (vector of measurement errors for phi_m), and `hp` (list of hyperparameters).
#' In addition, the field control is optional (see above).
#' @param calibDf A dataframe with radiocarbon calibration curve information
#'
#' @export
#'
#' @return soln, a list with three fields: prob (the input),
#' fit (the result of the call to stan),
#' and control (the control parameters used)
#'
bd_do_inference <- function(prob, calibDf,saveFile=NA) {
  # Unpack and/or define the control parameters
  if (exists("control", where = prob) == T) {
    haveNumChains <- exists("numChains", where = prob$control) == T
    haveSampsPerChain <- exists("sampsPerChain", where = prob$control) == T
    haveWarmup <- exists("warmup", where = prob$control) == T
    haveInitList <- exists("initList", where = prob$control) == T
    haveStanControl <- exists("stanControl", where = prob$control) == T
    #haveAdaptDelta <- exists("adaptDelta", where = prob$control) == T
  } else {
    haveNumChains <- F
    haveSampsPerChain <- F
    haveWarmup <- F
    haveInitList <- F
    haveStanControl <- F
    #haveAdaptDelta <- F
  }

  if (haveNumChains) {
    numChains <- prob$control$numChains
  } else {
    numChains <- 4
  }

  if (haveSampsPerChain) {
    sampsPerChain <- prob$control$sampsPerChain
  } else {
    sampsPerChain <- 2000
  }

  if (haveWarmup) {
    warmup <- prob$control$warmup
  } else {
    warmup <- floor(sampsPerChain/2)
  }


  if (haveInitList) {
    initList <- prob$control$initList
  } else {
    initList <- bd_sample_prior(prob$hp, numChains)
  }

  if (haveStanControl) {
    stanControl <- prob$control$stanControl
  } else {
    stanControl <- NA
  }

  controlFinal <- list(numChains = numChains,
                       sampsPerChain = sampsPerChain,
                       warmup = warmup,
                       initList = initList,
                       stanControl = stanControl)

  if (prob$hp$fitType == "gaussmix") {
    # Stan needs all the inputs and hyperparameters as variables in R's workspace
    ymin <- prob$hp$ymin
    ymax <- prob$hp$ymax
    mumin <- prob$hp$mumin
    mumax <- prob$hp$mumax
    ygrid <- seq(ymin, ymax, by = prob$hp$dy)
    M <- bd_calc_meas_matrix(ygrid, prob$phi_m, prob$sig_m, calibDf)
    Mt <- t(M)
    N <- dim(M)[1]
    G <- dim(M)[2]
    sigAlpha <- prob$hp$sigAlpha
    sigBeta <- prob$hp$sigBeta
    dirichParam <- prob$hp$dirichParam
    K <- prob$hp$K
    filePath <- system.file("stan/gaussmix.stan",
      package = "baydem"
    )
    options(mc.cores = parallel::detectCores())
    if (haveStanControl) {
      fit <- rstan::stan(filePath, chains = numChains, iter = sampsPerChain, warmup = warmup, init = initList, control = stanControl)
    } else {
      fit <- rstan::stan(filePath, chains = numChains, iter = sampsPerChain, warmup = warmup, init = initList)
    }
  } else {
    stop(paste("Unrecognized fit type:", prob$hp$fitType))
  }

  soln <- list(prob = prob, fit = fit, control = controlFinal)
  class(soln) <- "bd_soln"

  # If a save file was input, save the result to file. This is especially
  # useful for parallel batch runs
  if(!is.na(saveFile)) {
    saveRDS(soln,saveFile)
  }
  return(soln)
}
