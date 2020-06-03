#' @title Do Demographic Bayesian Inference
#'
#' @description This is the core function that implements the Bayesian inference.
#'
#' @details The input is a problem statement object (a list),
#' prob, that consists of the input data (the vectors phi_m and sig_m)
#' and the hyperparameters (hp).
#' stan is called via the rstan package to sample from the posterior.
#' The output is the variable soln of class bd_soln,
#' which is a list with the fields prob (the input) and fit (the result of the stan fit).
#' prob can also have an optional field control that specifies the
#' following control parameters for the Bayesian inference (default in parentheses):
#'   numChains     -- (4)      Number of chains
#'   sampsPerChain -- (2000)   Number of samples per chain
#'   warmup        -- (samp/2) Number of warmup samples (default is sampsPerChain/2)
#'   initSeed      --          An optional random number seed for initializing the chains
#'   stanSeed      --          An optional random number seed for the call to Stan
#'   initList      --          The initializations for each chain. The default is
#'                             to set this using a mixture fit to the summed
#'                             density
#'
#' @param prob List with the fields `phi_m` (vector of radiocarbon measurements as fraction modern),
#' `sig_m` (vector of measurement errors for phi_m), `hp` (list of hyperparameters), and
#' `calibDf`.
#' In addition, the field control is optional in prob (see above).
#'
#' @export
#'
#' @return soln, a list with three fields: prob (the input),
#' fit (the result of the call to stan),
#' and control (the control parameters used)
#'
bd_do_inference <- function(prob) {
  # Unpack and/or define the control parameters
  if (exists("control", where = prob) == T) {
    haveNumChains <- exists("numChains", where = prob$control) == T
    haveSampsPerChain <- exists("sampsPerChain", where = prob$control) == T
    haveWarmup <- exists("warmup", where = prob$control) == T
    haveInitSeed <- exists("initSeed", where = prob$control) == T
    haveStanSeed <- exists("stanSeed", where = prob$control) == T
    haveInitList <- exists("initList", where = prob$control) == T
    haveStanControl <- exists("stanControl", where = prob$control) == T
  } else {
    haveNumChains <- F
    haveSampsPerChain <- F
    haveWarmup <- F
    haveInitSeed <- F
    haveStanSeed <- F
    haveInitList <- F
    haveStanControl <- F
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
    warmup <- floor(sampsPerChain / 2)
  }

  if (haveInitSeed) {
    set.seed(prob$control$initSeed)
  }

  if (haveStanSeed) {
    stanSeed <- prob$control$stanSeed
  }

  if (haveStanControl) {
    stanControl <- prob$control$stanControl
  } else {
    stanControl <- NA
  }

  controlFinal <- list(
    numChains = numChains,
    sampsPerChain = sampsPerChain,
    warmup = warmup,
    stanControl = stanControl
  )

  if (prob$hp$fitType == "gaussmix") {
    # Stan needs all the inputs and hyperparameters as variables in R's workspace
    taumin <- prob$hp$taumin
    taumax <- prob$hp$taumax
    mumin <- prob$hp$mumin
    mumax <- prob$hp$mumax
    tau <- seq(taumin, taumax, by = prob$hp$dtau)
    M <- bd_calc_meas_matrix(tau, prob$phi_m, prob$sig_m, prob$calibDf)

    if (haveInitList) {
      initList <- prob$control$initList
    } else {
      # Set it using the summed density
      f_spdf <- colSums(M)
      # Sample 1000 times from the summed density to do a mixture fit
      xmix <- sample.int(length(f_spdf), 1000, replace = T, prob = f_spdf)
      gaussMix <- mixtools::normalmixEM(xmix, k = prob$hp$K, maxit = 20000)
      indSort <- order(gaussMix$mu)
      init0 <- list()
      init0$pi <- gaussMix$lambda[indSort]
      init0$mu <- taumin + (taumax - taumin) * (gaussMix$mu[indSort] - 1) / length(tau - 1)
      init0$sig <- gaussMix$sig[indSort] * prob$hp$dtau

      # Each chain needs an initialization for stan
      initList <- list()
      for (cc in 1:numChains) {
        initList[[cc]] <- init0
      }
    }
    controlFinal$initList <- initList

    Mt <- t(M)
    N <- dim(M)[1]
    G <- dim(M)[2]
    alpha_s <- prob$hp$alpha_s
    alpha_r <- prob$hp$alpha_r
    alpha_d <- prob$hp$alpha_d
    K <- prob$hp$K
    filePath <- system.file("stan/gaussmix.stan",
      package = "baydem"
    )
    options(mc.cores = parallel::detectCores())
    # There are four possible calls depending on whether haveStanControl is
    # TRUE and haveStanSeed is TRUE
    if(!haveStanSeed) {
      if (haveStanControl) {
        fit <- rstan::stan(filePath, chains = numChains, iter = sampsPerChain, warmup = warmup, init = initList, control = stanControl)
      } else {
        fit <- rstan::stan(filePath, chains = numChains, iter = sampsPerChain, warmup = warmup, init = initList)
      }
    } else { # do have Stan seed
      if (haveStanControl) {
        fit <- rstan::stan(filePath, chains = numChains, iter = sampsPerChain, warmup = warmup, seed=stanSeed, init = initList, control = stanControl)
      } else {
        fit <- rstan::stan(filePath, chains = numChains, iter = sampsPerChain, warmup = warmup, seed=stanSeed,init = initList)
      }
    }
  } else {
    stop(paste("Unrecognized fit type:", prob$hp$fitType))
  }

  soln <- list(prob = prob, fit = fit, control = controlFinal)
  class(soln) <- "bd_soln"

  # If a save file was input, save the result to file. This is especially
  # useful for parallel batch runs
  return(soln)
}
