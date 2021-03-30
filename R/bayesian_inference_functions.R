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
do_inference <- function(prob) {
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
    M <- calc_meas_matrix(tau, prob$phi_m, prob$sig_m, prob$calibDf)

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

#' @title Analyze the result of a call to \code{do_inference}
#'
#' @description \code{do_inference} calls Stan to do Bayesian inference by generating
#' a sample of parameters from the posterior of theta (or th).
#' This function analyzes the result of that inference.
#' In particular, it calculates the quantiles of the density function and growth rate.
#'
#' @details `soln` is the result of a call to \code{do_inference}.
#' It contains both the resulting samples and the parameters used in the inference,
#' such as the hyperparameters (see \code{do_inference} for further details).
#' The primary thing \code{analyze_soln} does is calculate quantiles of both
#' the parameterized density and growth rate. For example,
#' for a calendar date y each sample yields a density and growth rate.
#' The quantile is the value of the density or growth rate such that
#' a given proportion of samples are smaller than that value.
#' The probabilities used to calculate these quantiles are `probs = c(lev, 0.5, 1-lev)`,
#' where `lev` is the level (0.025 by default, so that 95% of the observations
#' lie between the first and last quantile bands).
#'
#' In addition, \code{analyze_soln} identifies calendar dates for which
#' the growth rate quantiles defined by `lev` and `1 - lev` do not contain zero.
#' This indicates significant positive or negative growth for the density curve.
#' The output vector `growthState` codes calendar dates by growth state as 'negative',
#' 'zero', and 'positive'. For the Gaussian mixture parameterization of the density,
#' the rate is not typically meaningful near the calendar date boundaries where it
#' increases linearly as the calendar date goes to positive or negative infinity.
#' The parameter `rateProp` provides control on how calendar dates are classified by
#' growth rate near these boundaries. In particular, the calendar dates with a cumulative
#' density (50% quantile) below `rateProp` (for the lower boundary) or above `1 - rateProp`
#' (for the upper boundary) are classified as 'missing' in `growthState`.
#' By default, `rateProp` is NA and no calendar dates are classified as missing.
#'
#' By default, a summary is done for each sample by calling
#' summarize_sample. This is not done of doSummary is FALSE
#'
#' @param soln The solution, a list-like object of class bd_soln (see \code{do_inference})
#' @param tau (optional) The calendar dates at which to evaluate densities.
#' If tau is not input, tau is built from the hyperparameters.
#' @param th_sim (optional) The known parameters used to create simulation data
#' @param lev (default: 0.025) The level to use for the quantile bands
#' @param rateProp (optional) The cumulative density needed to define rate growth bands
#' @param doSummary (default: `TRUE`) Whether to summarize each sample by calling summarize_sample
#'
#' @return A list with information on the quantiles of the density function and growth rate (and sample summaries)
#'
#' @export
analyze_soln <-
  function(soln,
           tau = NA,
           th_sim = NA,
           lev = 0.025,
           rateProp = NA,
           doSummary = T) {
    if (all(is.na(tau))) {
      tau <- seq(soln$prob$hp$taumin, soln$prob$hp$taumax, by = soln$prob$hp$dtau)
    }

    probs <- c(lev, 0.5, 1 - lev) # The probabilities to use for quantiles

    # Determine tau spacing, dtau, and ensure that tau is evenly spaced
    dtau <- unique(diff(tau))
    if (length(dtau) > 1) {
      stop("tau should by uniformily spaced")
    }

    # Extract the samples of theta in the variable TH. TH is matrix like object,
    # possibly of a specific class (e.g., gaussmix) with dimensions
    # numSamp x numParam, where numSamp is the number of samples and numParam is
    # the number of parameters (length of th).
    TH <- extract_param(soln$fit)

    numMix <- ncol(TH) / 3 # This assumes a gassian mixture fit. Future updates may generalize this
    numSamp <- nrow(TH)
    numGrid <- length(tau)

    # Calculate the pdf matrix, which is the density of the parametric model for
    # theta for each sample and each grid point. fMat has dimensions
    # N x G, where N is the number of samples in TH and G is the length of the vector tau.
    # Because calc_gauss_mix_pdf_mat is called with taumin and taumax, the density
    # is normalized to integrate to 1 on the interval taumin to taumax.
    fMat <- calc_gauss_mix_pdf_mat(TH, tau, taumin = soln$prob$hp$taumin, taumax = soln$prob$hp$taumax)

    # Calculate the rate for each sample and grid point (f' / f, where f is density)
    rateMat <- calc_gauss_mix_pdf_mat(TH, tau, taumin = soln$prob$hp$taumin, taumax = soln$prob$hp$taumax, type = "rate")

    # Calculate the quantiles of the normalized density matrix
    Qdens <- calc_quantiles(fMat, probs)

    # Normalized 50% densities (not normalized to integrate to 1)
    f50 <- Qdens[2, ] # The second row gives the 50% quantiles

    # Restrict to indices with enough probability mass (if necessary)
    if (!is.na(rateProp)) {
      rateInd <- which(cumsum(f50 * dtau) > rateProp & rev(cumsum(rev(f50) * dtau)) > rateProp)
    } else {
      rateInd <- 1:length(f50)
    }

    # Identify regions with growth rates that differ from zero per the input quantile level (lev)
    # growthState0 is -1 for significant negative growth, 1 for significant positive growth, and 0 otherwise
    Qrate <- calc_quantiles(rateMat[, rateInd], probs)
    growthState0 <- rep("zero", length(rateInd)) # growthState0 indices in rateInd
    growthState0[Qrate[2, ] > 0 & Qrate[1, ] > 0] <- "positive"
    growthState0[Qrate[2, ] < 0 & Qrate[3, ] < 0] <- "negative"
    growthState <- rep("missing", length(tau))
    growthState[rateInd] <- growthState0 # growthState for all indices


    # Calculate the measurement matrix
    M <- calc_meas_matrix(tau, soln$prob$phi_m, soln$prob$sig_m, soln$prob$calibDf)

    # Normalize by row
    M <- M / replicate(length(tau),rowSums(M)*dtau)
    # Calculate and normalize the summed probability density vector
    f_spdf <- colMeans(M)

    out <- list(
      tau = tau,
      f_spdf = f_spdf,
      Qdens = Qdens,
      Qrate = Qrate,
      probs = probs,
      rateProp = rateProp,
      rateInd = rateInd,
      growthState = growthState,
      dtau = dtau
    )
    class(out) <- "bd_analysis"

    if (doSummary) {
      summList <- list()
      for (n in 1:numSamp) {
        th <- TH[n, ]
        summList[[n]] <- summarize_trunc_gauss_mix_sample(th, soln$prob$hp$taumin, soln$prob$hp$taumax)
      }
      out$summList <- summList
    }

    haveSim <- !all(is.na(th_sim))
    if (haveSim) {
      f_sim <- calc_gauss_mix_pdf(th_sim, tau, taumin = soln$prob$hp$taumin, taumax = soln$prob$hp$taumax)
      rate_sim <- calc_gauss_mix_pdf(th_sim, tau, taumin = soln$prob$hp$taumin, taumax = soln$prob$hp$taumax, type = "rate")
      out$f_sim <- f_sim
      out$rate_sim <- rate_sim
    }
    return(out)
  }
#' @title Extract the Bayesian samples for a Gaussian mixture model generated by do_inference
#'
#' @description
#' The input fit is the result of a call to stan by \code{do_inference}, of
#' class stanfit. Return a matrix TH with dimensions S x (3*K), where S is the
#' number of samples (across all chains, and excluding warmup), and K is the
#' number of mixtures. The final column of as.matrix for class stanfit is the
#' log-posterior, which must be removed.
#'
#' @param fit The fit from stan, of class stanfit
#'
#' @return A matrix or list of samples
#' @export

extract_param <- function(fit) {
  if (class(fit) != "stanfit") {
    stop(paste("Expected fit to be class stanfit, but it is", class(fit)))
  }
  # as.matrix is defined for class stanfit and excludes warmup samples
  TH <- as.matrix(fit)
  # Remove the final column, which is the log-posterior, not a parameter
  TH <- TH[, -ncol(TH)]

  return(TH)
}

#' @title Identify growth periods and the peak value for a truncated Gaussian mixture
#'
#' @description
#' The input vector th parameterizes a Gaussian mixture, and taumin / taumax
#' give the limits of truncation. Summarize the sample by identifying growth /
#' decay periods and the peak value using the following procedure.
#'
#' (1) Calculate the derivative, f'(t), at the points t = seq(taumin,taumax,len=N),
#'     where N is 1000 by default.
#'
#' (2) Identify points where f'(t) changes sign, then numerically estimate the
#'     crossing point between the two t values where there was a sign change.
#'
#' (3) Create a vector of critical points, tcrit, which includes taumin / taumax
#'     as well as the crossing points found in the preceding step.
#'
#' (4) Calculate the density at the critical points to identify the peak value,
#'     fpeak, and corresponding calendar date, tpeak, as well as the index of
#'     of the peak in tcrit, indPeak.
#'
#' (5) For each time period (the length(tpeak)-1 durations defined by ypeak)
#'     determine the sign of the density function, f(t), and create a character
#'     vector, slope, that has the value 'pos' if f(t) is positive and 'neg' if
#'     f(t) is negative.
#'
#' (6) Finally, create a character vector, pattern, that appends the index of
#'     the peak in tcrit (converted to a character) to the character vector
#'     slope. This defines a unique pattern of the sample that takes into
#'     account periods of growth / decline and the relative location of the
#'     peak.
#'
#' @param  th The Gaussian mixture parameterization
#' @param  taumin The lower truncation value
#' @param  taumax The upper truncation value
#' @param  N (Default 1000) The number of points use for identifying slope changes
#'
#' @return A list consisting of tlo / thi (specifying the time periods), indPeak, tpeak, fpeak, and pattern (see Description)
#'
#' @export
summarize_trunc_gauss_mix_sample <- function(th, taumin, taumax, N = 1000) {
  # (1) Calculate the derivative of the density
  K <- length(th) / 3 # Number of mixtures
  t <- seq(taumin, taumax, len = N)
  fprime <- calc_gauss_mix_pdf(th, t, taumin, taumax, type = "derivative")

  # (2) Identify locations in t where the derivative changes sign. This happens
  #     if fprime[n] * fprime[n+1] is less than zero. Then, numerically
  #     estimate the exact t-value of the crossing.
  ind <- which(fprime[1:(length(fprime) - 1)] * fprime[2:length(fprime)] < 0)
  M <- length(ind) # Number of cross-overs

  # Vectors for t / f values of crossings
  tcross <- rep(NA, M)
  fcross <- rep(NA, M)

  if (M > 0) {
    # Objective function to maximize
    rootFun <- function(t) {
      return(calc_gauss_mix_pdf(th, t, taumin, taumax, type = "derivative"))
    }

    # Iterate over crossings
    for (m in 1:M) {
      root <- stats::uniroot(rootFun, lower = t[ind[m]], upper = t[ind[m] + 1])
      tcross[m] <- root$root
      fcross[m] <- calc_gauss_mix_pdf(th, tcross[m], taumin, taumax)
    }
  }

  # (3-4) Create the vector of critical points, calculate densities, and
  #       identify peak
  tcrit <- c(taumin, tcross, taumax)
  fcrit <- c(calc_gauss_mix_pdf(th, taumin, taumin, taumax), fcross, calc_gauss_mix_pdf(th, taumin, taumin, taumax))
  indPeak <- which.max(fcrit)
  tpeak <- tcrit[indPeak]
  fpeak <- fcrit[indPeak]

  # (5) Create tlo, thi, and slope
  numPer <- length(tcrit) - 1 # Number of periods
  tlo <- tcrit[1:numPer]
  thi <- tcrit[2:(numPer + 1)]
  df <- diff(fcrit)
  slope <- rep("pos", numPer)
  slope[df < 0] <- "neg"

  # (6) Create the pattern (then return the result)
  pattern <- c(slope, as.character(indPeak))

  return(list(periods = data.frame(tlo = tlo, thi = thi, slope = slope), indPeak = indPeak, tpeak = tpeak, fpeak = fpeak, pattern = pattern))
}

#' @title Calculate the quantiles for an input matrix X
#'
#' @description The input matrix X has dimensions S x G, where S is the number
#'              of samples and G the number of grid points at which X was
#'              evaluated. Calculate quantiles for each grid point, g = 1,2,..G.
#'
#' @param X The matrix for which quantiles are calculated, with dimension S x G
#' @param probs The probability values at which to calculate the quantiles (default: `c(0.025, 0.5, 0.975)`)
#'
#' @return The quantiles, a matrix with dimension length(probs) x G
#'
#' @export
calc_quantiles <- function(X, probs = c(.025, .5, .975)) {
  numQuant <- length(probs) # Number of quantiles
  G <- dim(X)[2] # Number of grid points

  Q <- matrix(NA, numQuant, G) # Initialize Q with dimensions numQuant x G
  # Iterate over grid points to calculate quantiles
  for (g in 1:G) {
    Q[, g] <- stats::quantile(X[, g], probs = probs)
  }
  return(Q)
}

#' @title For each sample, calculate the time it takes for the density to decrease by half from the peak
#'
#' @details
#' For each sample, calculate the time it takes for the density to decrease by
#' half from the peak. Optionally, a different proportion can be used than the
#' default propChange = 0.5. For example, with propChange = 0.1 the time it
#' takes for the density to decrease by 10% is used. If the relative density is
#' not reached, the half life for the sample is set to NA. If there is no
#' interior peak in the range peakRange, which is taumin to taumax by default,
#' the half life is set to NA.
#'
#' @param soln The result of a call to do_inference
#' @param propChange (Default 0.5) The relative decrease in density to use for the duration calculation
#' @param anal (Optional) The result of a call to analyze_soln. If not provided, it is calculated
#' @param peakRange (default: `c(taumin, taumax)`) peakRange can be given so that the peak density used is on the range peakRange
#'
#' @return A vector of "half-lives" (proportional change set by propChange)
#'
#' @export
calc_half_life_from_peak <-
  function(soln,
           propChange = 0.5,
           anal = NA,
           peakRange = NA) {
    TH <- extract_param(soln$fit)
    N <- nrow(TH)
    taumin <- soln$prob$hp$taumin
    taumax <- soln$prob$hp$taumax
    dtau <- soln$prob$hp$dtau

    if (all(is.na(anal))) {
      anal <- analyze_soln(soln)
    }
    summList <- anal$summList

    if (all(is.na(peakRange))) {
      peakRange <- c(taumin, taumax)
    }


    halfLife <- rep(NA, N)
    for (n in 1:N) {
      th <- TH[n, ]
      # Identify the peak, ensuring it is on peakRange

      # critical points
      tcrit <- c(summList[[n]]$periods$tlo, summList[[n]]$periods$thi[length(summList[[n]]$periods$thi)])
      tcrit <- tcrit[peakRange[1] <= tcrit & tcrit <= peakRange[2]]
      fcrit <- calc_gauss_mix_pdf(th, tcrit, taumin, taumax)
      indPeak <- which.max(fcrit)
      tpeak <- tcrit[indPeak]
      fpeak <- fcrit[indPeak]
      isIn <- taumin < tpeak && tpeak < taumax
      if (isIn) {
        # Function for root finder
        rootFun <- function(t) {
          return(fpeak * propChange - calc_gauss_mix_pdf(th, t, taumin, taumax, type = "density"))
        }

        # Find root. Catch any errors in case the half life does not exist on the
        # interval tpeak to taumax
        result <- tryCatch(
          {
            root <- stats::uniroot(rootFun,
              lower = tpeak,
              upper = peakRange[2]
            )
            halfLife[n] <- min(root$root - tpeak)
          },
          error = function(e) {
            NA
          }
        )
      }
    }
    return(halfLife)
  }

#' @title Calculate the relative density at two dates (or a range of dates / the peak)
#'
#' @description
#'
#' Calculate the relative density at two dates, and/or a range of dates and/or the peak value (see details).
#'
#' @details
#' Calculate the relative density for two dates or, more generally, for two
#' different specifications of the density aside from a simple date. The
#' additional specifications that are supported are the peak value and the mean
#' density on an interval. For a simple date, spec1/spec2 should be scalar
#' real numbers. For a date range, spec1/spec2 should be real vectors with a
#' length of 2. For the peak, spec1/spec2 should be the string 'peak'.
#'
#' By default, this calculation is done for all the Bayesian samples in soln,
#' which is the result of a call to do_inference. Optionally, a subset can be
#' specified via the input ind, which should be a vector of integer indices at
#' which to do the calculation. To save computation if either spec1 or spec2 is
#' 'peak', the result of a call to analyze_soln for which doSummary was T
#' can be input.
#'
#'
#' @param soln The result of a call to do_inference
#' @param anal The result of a call to analyze_soln
#' @param spec1 The specification for the first density (see details)
#' @param spec2 The specification for the second density (see details)
#' @param ind (Optional) Indices at which to do the calculation. By default,
#'            all the samples in anal are used.
#' @param anal (Optional) The result of a call to analyze_soln. This is only
#'             needed if either spec1 or spec2 is 'peak'
#'
#' @return A vector of relative densities (f_spec1 / f_spec2)
#'
#' @export
calc_relative_density <- function(soln, spec1, spec2, ind = NA, anal = NA) {
  TH <- extract_param(soln$fit)
  N <- nrow(TH)
  if (all(is.na(ind))) {
    ind <- 1:N
  }


  # Interpret and do error checking on inputs by calling helper function below
  spec1 <- unpack_spec(spec1, soln, T)
  spec2 <- unpack_spec(spec2, soln, F)

  if (spec1$type == "peak" || spec2$type == "peak") {
    if (all(is.na(anal))) {
      anal <- analyze_soln(soln) # If ind is not NA, this may involve unused computation
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
      if (point < soln$prob$hp$taumin || soln$prob$hp$taumax < point) {
        stop(paste(s, "is a single date, but not in the range taumin to taumax"))
      }
      return(list(type = "point", value = point))
    } else if (length(spec) == 2) { # Numeric / length 2
      lower <- spec[1]
      if (lower < soln$prob$hp$taumin || soln$prob$hp$taumax < lower) {
        stop(paste(s, "is a date range, but lower value is not in the range taumin to taumax"))
      }
      upper <- spec[2]
      if (upper < soln$prob$hp$taumin || soln$prob$hp$taumax < upper) {
        stop(paste(s, "is a date range, but upper value is not in the range taumin to taumax"))
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
calc_point_density <- function(TH, soln, t) {
  return(as.numeric(calc_gauss_mix_pdf_mat(TH, t, taumin = soln$prob$hp$taumin, taumax = soln$prob$hp$taumax)))
}


# A helper function to calculate the mean density over a range
calc_range_density <- function(TH, soln, tlo, thi) {
  flo <- as.numeric(calc_gauss_mix_pdf_mat(TH, tlo, taumin = soln$prob$hp$taumin, taumax = soln$prob$hp$taumax, type = "cumulative"))
  fhi <- as.numeric(calc_gauss_mix_pdf_mat(TH, thi, taumin = soln$prob$hp$taumin, taumax = soln$prob$hp$taumax, type = "cumulative"))
  return((fhi - flo) / (thi - tlo))
}


# A helper function to calculate the peak density
calc_peak_density <- function(summList) {
  return(unlist(lapply(summList, function(x) {
    x$fpeak
  })))
}
