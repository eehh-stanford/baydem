#' @title Do Demographic Bayesian Inference
#'
#' @description This is the core function that implements the Bayesian
#' inference.
#'
#' @details The input is a problem statement object (a list),
#' prob, that consists of the input data (the vectors phi_m and sig_m)
#' and the hyperparameters (hp).
#' stan is called via the rstan package to sample from the posterior.
#' The output is the variable soln of class bd_soln,
#' which is a list with the fields prob (the input) and fit (the result of the
#' stan fit).
#'
#' prob can also have an optional field control that specifies the
#' following control parameters for the Bayesian inference (defaults in
#' parentheses):
#'
#'   num_chains      -- Number of chains (4)
#'   samps_per_chain -- Number of samples per chain (2000)
#'   warmup          -- Number of warmup samples(default is samps_per_chain/2)
#'   init_seed        -- An optional random number seed for initializing the
#'                      chains
#'   stan_seed        -- An optional random number seed for the call to Stan
#'   init_list        -- The initializations for each chain. The default is to
#'                      set this using a mixture fit of the summed density
#'
#' @param prob List with the fields `phi_m` (vector of fraction modern values
#'   for the radiocarbon measurements, `sig_m` (measurement errors for phi_m),
#'   `hp` (list of hyperparameters), and `calibDf`.
#' @param control Control parameters (see details)
#'
#' @return soln, a list with three fields: prob (the input), fit (the result of
#'   the call to stan), and control (the control parameters used)
#'
#' @export
do_inference <- function(prob) {
  # Unpack and/or define the control parameters
  if (exists("control", where = prob) == T) {
    have_num_chains <- exists("num_chains", where = prob$control) == T
    have_samps_per_chain <- exists("samps_per_chain", where = prob$control) == T
    have_warmup <- exists("warmup", where = prob$control) == T
    have_init_seed <- exists("init_seed", where = prob$control) == T
    have_stan_seed <- exists("stan_seed", where = prob$control) == T
    have_init_list <- exists("init_list", where = prob$control) == T
    have_stan_control <- exists("stan_control", where = prob$control) == T
  } else {
    have_num_chains <- F
    have_samps_per_chain <- F
    have_warmup <- F
    have_init_seed <- F
    have_stan_seed <- F
    have_init_list <- F
    have_stan_control <- F
  }

  if (have_num_chains) {
    num_chains <- prob$control$num_chains
  } else {
    num_chains <- 4
  }

  if (have_samps_per_chain) {
    samps_per_chain <- prob$control$samps_per_chain
  } else {
    samps_per_chain <- 2000
  }

  if (have_warmup) {
    warmup <- prob$control$warmup
  } else {
    warmup <- floor(samps_per_chain / 2)
  }

  if (have_init_seed) {
    set.seed(prob$control$init_seed)
  }

  if (have_stan_seed) {
    stan_seed <- prob$control$stan_seed
  }

  if (have_stan_control) {
    stan_control <- prob$control$stan_control
  } else {
    stan_control <- NA
  }

  control_final <- list(
    num_chains = num_chains,
    samps_per_chain = samps_per_chain,
    warmup = warmup,
    stan_control = stan_control
  )

  if (prob$hp$fit_type == "gauss_mix") {
    # Stan needs all the inputs and hyperparameters as variables in R's
    # workspace
    tau_min <- prob$hp$tau_min
    tau_max <- prob$hp$tau_max
    mu_min <- prob$hp$mu_min
    mu_max <- prob$hp$mu_max
    tau <- seq(tau_min, tau_max, by = prob$hp$dtau)
    M <- calc_meas_matrix(tau, prob$phi_m, prob$sig_m, prob$calib_df)

    if (have_init_list) {
      init_list <- prob$control$init_list
    } else {
      # Set it using the summed density
      f_spdf <- colSums(M)
      # Sample 1000 times from the summed density to do a mixture fit
      x_mix <- sample.int(length(f_spdf), 1000, replace = T, prob = f_spdf)
      gauss_mix <- mixtools::normalmixEM(x_mix, k = prob$hp$K, maxit = 20000)
      ind_sort <- order(gauss_mix$mu)
      init0 <- list()
      init0$pi <- gauss_mix$lambda[ind_sort]
      init0$mu <- tau_min + (tau_max - tau_min) *
        (gauss_mix$mu[ind_sort] - 1) / length(tau - 1)
      init0$sig <- gauss_mix$sig[ind_sort] * prob$hp$dtau

      # Each chain needs an initialization for stan
      init_list <- list()
      for (cc in 1:num_chains) {
        init_list[[cc]] <- init0
      }
    }
    control_final$init_list <- init_list

    Mt <- t(M)
    N <- dim(M)[1]
    G <- dim(M)[2]
    alpha_s <- prob$hp$alpha_s
    alpha_r <- prob$hp$alpha_r
    alpha_d <- prob$hp$alpha_d
    K <- prob$hp$K
    file_path <- system.file("stan/gaussmix.stan",
      package = "baydem"
    )
    options(mc.cores = parallel::detectCores())
    # There are four possible calls depending on whether haveStanControl is
    # TRUE and haveStanSeed is TRUE
    if(!have_stan_seed) {
      if (have_stan_control) {
        fit <- rstan::stan(file_path,
                           chains = num_chains,
                           iter = samps_per_chain,
                           warmup = warmup,
                           init = init_list,
                           control = stan_control)
      } else {
        fit <- rstan::stan(file_path,
                           chains = num_chains,
                           iter = samps_per_chain,
                           warmup = warmup,
                           init = init_list)
      }
    } else { # do have Stan seed
      if (have_stan_control) {
        fit <- rstan::stan(file_path,
                           chains = num_chains,
                           iter = samps_per_chain,
                           warmup = warmup,
                           seed=stan_seed,
                           init = init_list,
                           control = stan_control)
      } else {
        fit <- rstan::stan(file_path,
                           chains = num_chains,
                           iter = samps_per_chain,
                           warmup = warmup,
                           seed=stan_seed,
                           init = init_list)
      }
    }
  } else {
    stop(paste("Unrecognized fit type:", prob$hp$fit_type))
  }

  soln <- list(prob = prob, fit = fit, control = control_final)
  class(soln) <- "bd_soln"

  # If a save file was input, save the result to file. This is especially
  # useful for parallel batch runs
  return(soln)
}

#' @title Analyze the result of a call to \code{do_inference}
#'
#' @description \code{do_inference} calls Stan to do Bayesian inference by
#' generating a sample of parameters from the posterior of theta (or th).
#' This function analyzes the result of that inference. In particular, it
#' calculates the quantiles of the density function and growth rate.
#'
#' @details `soln` is the result of a call to \code{do_inference}.
#' It contains both the resulting samples and the parameters used in the
#' inference, such as the hyperparameters (see \code{do_inference} for further
#' details). The primary thing \code{analyze_soln} does is calculate quantiles
#' of both the parameterized density and growth rate. For example, for a
#' calendar date y each sample yields a density and growth rate. The quantile is
#' the value of the density or growth rate such that a given proportion of
#' samples are smaller than that value. The probabilities used to calculate
#' these quantiles are `probs = c(lev, 0.5, 1-lev)`, where `lev` is the level
#' (0.025 by default, so that 95% of the observations lie between the first and
#' last quantile bands).
#'
#' In addition, \code{analyze_soln} identifies calendar dates for which
#' the growth rate quantiles defined by `lev` and `1 - lev` do not contain zero.
#' This indicates significant positive or negative growth for the density curve.
#' The output vector `growth_state` codes calendar dates by growth state as
#' 'negative', 'zero', and 'positive'. For the Gaussian mixture parameterization
#' of the density, the rate is not typically meaningful near the calendar date
#' boundaries where it increases linearly as the calendar date goes to positive
#' or negative infinity. The parameter `rate_prop` provides control on how
#' calendar dates are classified by growth rate near these boundaries. In
#' particular, the calendar dates with a cumulative density (50% quantile) below
#' `rate_prop` (for the lower boundary) or above `1 - rate_prop` (for the upper
#' boundary) are classified as 'missing' in `growth_state`. By default,
#' `rate_prop` is NA and no calendar dates are classified as missing.
#'
#' By default, a summary is done for each sample by calling summarize_sample.
#' This is not done if doSummary is FALSE.
#'
#' @param soln The solution, a list-like object of class bd_soln (see
#'   \code{do_inference})
#' @param tau (optional) The calendar dates at which to evaluate densities. If
#'   tau is not input, tau is built from the hyperparameters.
#' @param th_sim (optional) The known parameters used to create simulation data
#' @param lev (default: 0.025) The level to use for the quantile bands
#' @param rate_prop (optional) The cumulative density needed to define rate
#'   growth bands
#' @param do_summary (default: `TRUE`) Whether to summarize each sample by
#'   calling summarize_sample
#'
#' @return A list with information on the quantiles of the density function and
#'   growth rate (and sample summaries)
#'
#' @export

analyze_soln <- function(soln,
                         tau = NA,
                         th_sim = NA,
                         lev = 0.025,
                         rate_prop = NA,
                         do_summary = T) {

    if (all(is.na(tau))) {
      tau <- seq(soln$prob$hp$tau_min,
                 soln$prob$hp$tau_max,
                 by = soln$prob$hp$dtau)
    }

    probs <- c(lev, 0.5, 1 - lev) # The probabilities to use for quantiles

    # Determine tau spacing, dtau, and ensure that tau is evenly spaced
    dtau <- unique(diff(tau))
    if (length(dtau) > 1) {
      stop("tau should by uniformily spaced")
    }

    # Extract the samples of theta in the variable TH. TH is matrix like object,
    # possibly of a specific class (e.g., gauss_mix) with dimensions
    # num_samp x num_param, where num_samp is the number of samples and
    # num_param is the number of parameters (length of th).
    TH <- extract_param(soln$fit)

    num_mix <- ncol(TH) / 3 # This assumes a gassian mixture fit. Future updates
                            # may generalize this
    num_samp <- nrow(TH)
    num_grid <- length(tau)

    # Calculate the pdf matrix, which is the density of the parametric model for
    # theta for each sample and each grid point. f_mat has dimensions N x G,
    # where N is the number of samples in TH and G is the length of the vector
    # tau. Because calc_gauss_mix_pdf_mat is called with tau_min and tau_max,
    # the density is normalized to integrate to 1 on the interval tau_min to
    # tau_max.
    f_mat <- calc_gauss_mix_pdf_mat(TH,
                                    tau,
                                    tau_min = soln$prob$hp$tau_min,
                                    tau_max = soln$prob$hp$tau_max)

    # Calculate the rate for each sample and grid point (f' / f, where f is
    # density)
    rate_mat <- calc_gauss_mix_pdf_mat(TH,
                                       tau,
                                       tau_min = soln$prob$hp$tau_min,
                                       tau_max = soln$prob$hp$tau_max,
                                       type = "rate")

    # Calculate the quantiles of the normalized density matrix
    Qdens <- calc_quantiles(f_mat, probs)

    # Normalized 50% densities (not normalized to integrate to 1)
    f50 <- Qdens[2, ] # The second row gives the 50% quantiles

    # Restrict to indices with enough probability mass (if necessary)
    if (!is.na(rate_prop)) {
      rate_ind <- which(cumsum(f50 * dtau) > rate_prop &
                          rev(cumsum(rev(f50) * dtau)) > rate_prop)
    } else {
      rate_ind <- 1:length(f50)
    }

    # Identify regions with growth rates that differ from zero per the input
    # quantile level (lev) growth_state0 is -1 for significant negative growth,
    # 1 for significant positive growth, and 0 otherwise
    Qrate <- calc_quantiles(rate_mat[, rate_ind], probs)
    growth_state0 <- rep("zero", length(rate_ind)) # growthState0 indices in
                                                  # rate_ind
    growth_state0[Qrate[2, ] > 0 & Qrate[1, ] > 0] <- "positive"
    growth_state0[Qrate[2, ] < 0 & Qrate[3, ] < 0] <- "negative"
    growth_state <- rep("missing", length(tau))
    growth_state[rate_ind] <- growth_state0 # growthState for all indices


    # Calculate the measurement matrix
    M <- calc_meas_matrix(tau,
                          soln$prob$phi_m,
                          soln$prob$sig_m,
                          soln$prob$calib_df)

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
      rate_prop = rate_prop,
      rate_ind = rate_ind,
      growth_state = growth_state,
      dtau = dtau
    )
    class(out) <- "bd_analysis"

    if (do_summary) {
      summ_list <- list()
      for (n in 1:num_samp) {
        th <- TH[n, ]
        summ_list[[n]] <- summarize_trunc_gauss_mix_sample(th,
                                                           soln$prob$hp$tau_min,
                                                           soln$prob$hp$tau_max)
      }
      out$summ_list <- summ_list
    }

    have_sim <- !all(is.na(th_sim))
    if (have_sim) {
      f_sim <- calc_gauss_mix_pdf(th_sim,
                                  tau,
                                  tau_min = soln$prob$hp$tau_min,
                                  tau_max = soln$prob$hp$tau_max)
      rate_sim <- calc_gauss_mix_pdf(th_sim,
                                     tau,
                                     tau_min = soln$prob$hp$tau_min,
                                     tau_max = soln$prob$hp$tau_max,
                                     type = "rate")
      out$f_sim <- f_sim
      out$rate_sim <- rate_sim
    }
    return(out)
  }
#' @title Extract the Bayesian samples for a Gaussian mixture model generated by
#' do_inference
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

#' @title Identify growth periods and the peak value for a truncated Gaussian
#' mixture
#'
#' @description
#' The input vector th parameterizes a Gaussian mixture, and tau_min / tau_max
#' give the limits of truncation. Summarize the sample by identifying growth /
#' decay periods and the peak value using the following procedure.
#'
#' (1) Calculate the derivative, f'(t), at the points t =
#'     seq(tau_min,tau_max,len=N), where N is 1000 by default.
#'
#' (2) Identify points where f'(t) changes sign, then numerically estimate the
#'     crossing point between the two t values where there was a sign change.
#'
#' (3) Create a vector of critical points, t_crit, which includes
#'     tau_min / tau_max as well as the crossing points found in the preceding
#'     step.
#'
#' (4) Calculate the density at the critical points to identify the peak value,
#'     f_peak, and corresponding calendar date, t_peak, as well as the index of
#'     of the peak in t_crit, ind_peak.
#'
#' (5) For each time period (the length(t_peak)-1 durations defined by y_peak)
#'     determine the sign of the density function, f(t), and create a character
#'     vector, slope, that has the value 'pos' if f(t) is positive and 'neg' if
#'     f(t) is negative.
#'
#' (6) Finally, create a character vector, pattern, that appends the index of
#'     the peak in t_crit (converted to a character) to the character vector
#'     slope. This defines a unique pattern of the sample that takes into
#'     account periods of growth / decline and the relative location of the
#'     peak.
#'
#' @param  th The Gaussian mixture parameterization
#' @param  tau_min The lower truncation value
#' @param  tau_max The upper truncation value
#' @param  N (Default 1000) The number of points use for identifying slope
#'   changes
#'
#' @return A list consisting of t_lo / t_hi (specifying the time periods),
#'   ind_peak, t_peak, f_peak, and pattern (see Description)
#'
#' @export
summarize_trunc_gauss_mix_sample <- function(th,
                                             tau_min,
                                             tau_max,
                                             N = 1000) {
  # (1) Calculate the derivative of the density
  K <- length(th) / 3 # Number of mixtures
  t <- seq(tau_min, tau_max, len = N)
  f_prime <- calc_gauss_mix_pdf(th, t, tau_min, tau_max, type = "derivative")

  # (2) Identify locations in t where the derivative changes sign. This happens
  #     if f_prime[n] * f_prime[n+1] is less than zero. Then, numerically
  #     estimate the exact t-value of the crossing.
  ind <- which(f_prime[1:(length(f_prime) - 1)] *
                 f_prime[2:length(f_prime)] < 0)
  M <- length(ind) # Number of cross-overs

  # Vectors for t / f values of crossings
  t_cross <- rep(NA, M)
  f_cross <- rep(NA, M)

  if (M > 0) {
    # Objective function to maximize
    root_fun <- function(t) {
      return(calc_gauss_mix_pdf(th, t, tau_min, tau_max, type = "derivative"))
    }

    # Iterate over crossings
    for (m in 1:M) {
      root <- stats::uniroot(root_fun, lower = t[ind[m]], upper = t[ind[m] + 1])
      t_cross[m] <- root$root
      f_cross[m] <- calc_gauss_mix_pdf(th, t_cross[m], tau_min, tau_max)
    }
  }

  # (3-4) Create the vector of critical points, calculate densities, and
  #       identify peak
  t_crit <- c(tau_min, t_cross, tau_max)
  f_crit <- c(calc_gauss_mix_pdf(th, tau_min, tau_min, tau_max),
              f_cross,
              calc_gauss_mix_pdf(th, tau_max, tau_min, tau_max))
  ind_peak <- which.max(f_crit)
  t_peak <- t_crit[ind_peak]
  f_peak <- f_crit[ind_peak]

  # (5) Create tlo, thi, and slope
  num_per <- length(t_crit) - 1 # Number of periods
  t_lo <- t_crit[1:num_per]
  t_hi <- t_crit[2:(num_per + 1)]
  df <- diff(f_crit)
  slope <- rep("pos", num_per)
  slope[df < 0] <- "neg"

  # (6) Create the pattern (then return the result)
  pattern <- c(slope, as.character(ind_peak))

  return(list(periods = data.frame(t_lo = t_lo,
                                   t_hi = t_hi,
                                   slope = slope),
              ind_peak = ind_peak,
              t_peak = t_peak,
              f_peak = f_peak,
              pattern = pattern))
}

#' @title Calculate the quantiles for an input matrix X
#'
#' @description The input matrix X has dimensions S x G, where S is the number
#'              of samples and G the number of grid points at which X was
#'              evaluated. Calculate quantiles for each grid point, g = 1,2,..G.
#'
#' @param X The matrix for which quantiles are calculated, with dimension S x G
#' @param probs The probability values at which to calculate the quantiles
#'   (default: `c(0.025, 0.5, 0.975)`)
#'
#' @return The quantiles, a matrix with dimension length(probs) x G
#'
#' @export
calc_quantiles <- function(X, probs = c(.025, .5, .975)) {
  num_quant <- length(probs) # Number of quantiles
  G <- dim(X)[2] # Number of grid points

  Q <- matrix(NA, num_quant, G) # Initialize Q with dimensions numQuant x G
  # Iterate over grid points to calculate quantiles
  for (g in 1:G) {
    Q[, g] <- stats::quantile(X[, g], probs = probs)
  }
  return(Q)
}

#' @title For each sample, calculate the time it takes for the density to
#' decrease by half from the peak
#'
#' @details
#' For each sample, calculate the time it takes for the density to decrease by
#' half from the peak. Optionally, a different proportion can be used than the
#' default prop_change = 0.5. For example, with prop_change = 0.1 the time it
#' takes for the density to decrease by 10% is used. If the relative density is
#' not reached, the half life for the sample is set to NA. If there is no
#' interior peak in the range peak_range, which is tau_min to tau_max by
#' default, the half life is set to NA.
#'
#' @param soln The result of a call to do_inference
#' @param prop_change (Default 0.5) The relative decrease in density to use for
#'   the duration calculation
#' @param anal (Optional) The result of a call to analyze_soln. If not provided,
#'   it is calculated
#' @param peak_range (default: `c(tau_min, tau_max)`) peak_range can be given so
#'   that the peak density used is on the range peak_range
#'
#' @return A vector of "half-lives" (proportional change set by prop_change)
#'
#' @export
calc_half_life_from_peak <-
  function(soln,
           prop_change = 0.5,
           anal = NA,
           peak_range = NA) {
    TH <- extract_param(soln$fit)
    N <- nrow(TH)
    tau_min <- soln$prob$hp$tau_min
    tau_max <- soln$prob$hp$tau_max
    dtau <- soln$prob$hp$dtau

    if (all(is.na(anal))) {
      anal <- analyze_soln(soln)
    }
    summ_list <- anal$summ_list

    if (all(is.na(peak_range))) {
      peak_range <- c(tau_min, tau_max)
    }

    half_life <- rep(NA, N)
    for (n in 1:N) {
      th <- TH[n, ]
      # Identify the peak, ensuring it is in peak_range

      # critical points
      t_crit <-
        c(summ_list[[n]]$periods$t_lo,
          summ_list[[n]]$periods$t_hi[length(summ_list[[n]]$periods$t_hi)])
      t_crit <- t_crit[peak_range[1] <= t_crit & t_crit <= peak_range[2]]
      f_crit <- calc_gauss_mix_pdf(th, t_crit, tau_min, tau_max)
      ind_peak <- which.max(f_crit)
      t_peak <- t_crit[ind_peak]
      f_peak <- f_crit[ind_peak]
      is_in <- tau_min < t_peak && t_peak < tau_max
      if (is_in) {
        # Function for root finder
        root_fun <- function(t) {
          return(f_peak * prop_change - calc_gauss_mix_pdf(th,
                                                           t,
                                                           tau_min,
                                                           tau_max,
                                                           type = "density"))
        }

        # Find root. Catch any errors in case the half life does not exist on the
        # interval tpeak to taumax
        result <- tryCatch(
          {
            root <- stats::uniroot(root_fun,
                                   lower = t_peak,
                                   upper = peak_range[2]
            )
            half_life[n] <- min(root$root - t_peak)
          },
          error = function(e) {
            NA
          }
        )
      }
    }
    return(half_life)
  }

#' @title Calculate the relative density at two dates (or a range of dates / the
#' peak)
#'
#' @description
#'
#' Calculate the relative density at two dates, and/or a range of dates and/or
#' the peak value (see details).
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
#' 'peak', the result of a call to analyze_soln for which do_summary was TRUE
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
    summ_list <- anal$summ_list[ind]
  }

  # Calculate the density for spec1
  if (spec1$type == "point") {
    f1 <- calc_point_density(TH[ind, ], soln, spec1$value)
  } else if (spec1$type == "range") {
    f1 <- calc_range_density(TH[ind, ], soln, spec1$lower, spec1$upper)
  } else if (spec1$type == "peak") {
    f1 <- calc_peak_density(summ_list)
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
    f2 <- calc_peak_density(summ_list)
  } else {
    # This should not happen, but throw an error regardless
    stop("Unsupported spec type")
  }

  return(f1 / f2)
}

# A helper function to unpack and do error checking on inputs spec1 / spec2
unpack_spec <- function(spec, soln, is_one) {
  # For more informative error messages, use the input is_one to set the string
  # s to spec1 or spec2
  if (is_one) {
    s <- "spec1"
  } else {
    s <- "spec2"
  }

  # Handle the supported cases, throwing an error if necessary
  if (is.numeric(spec)) {
    if (length(spec) == 1) { # Numeric / length 1
      point <- spec
      if (point < soln$prob$hp$tau_min || soln$prob$hp$tau_max < point) {
        stop(
          paste(s, "is a single date, but not in the range tau_min to tau_max"))
      }
      return(list(type = "point", value = point))
    } else if (length(spec) == 2) { # Numeric / length 2
      lower <- spec[1]
      if (lower < soln$prob$hp$tau_min || soln$prob$hp$tau_max < lower) {
        stop(paste(s,"is a date range, but lower value is not ",
                   "in the range tau_min to tau_max"))
      }
      upper <- spec[2]
      if (upper < soln$prob$hp$tau_min || soln$prob$hp$tau_max < upper) {
        stop(paste(s, "is a date range, but upper value is not ",
                   "in the range taumin to taumax"))
      }
      if (lower > upper) {
        stop(paste(s, "is a date range, but lower value is ",
                   "greater than upper value"))
      }
      return(list(type = "range", lower = lower, upper = upper))
    } else { # Numeric / not length 1 or 2
      stop(
        paste(s, "is numeric, but is neither a single date nor a date range"))
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
  return(as.numeric(calc_gauss_mix_pdf_mat(TH,
                                           t,
                                           tau_min = soln$prob$hp$tau_min,
                                           tau_max = soln$prob$hp$tau_max)))
}


# A helper function to calculate the mean density over a range
calc_range_density <- function(TH, soln, t_lo, t_hi) {
  f_lo <- as.numeric(calc_gauss_mix_pdf_mat(TH,
                                            t_lo,
                                            tau_min = soln$prob$hp$tau_min,
                                            tau_max = soln$prob$hp$tau_max,
                                            type = "cumulative"))
  f_hi <- as.numeric(calc_gauss_mix_pdf_mat(TH,
                                            t_hi,
                                            tau_min = soln$prob$hp$tau_min,
                                            tau_max = soln$prob$hp$tau_max,
                                            type = "cumulative"))
  return((f_hi - f_lo) / (t_hi - t_lo))
}


# A helper function to calculate the peak density
calc_peak_density <- function(summ_list) {
  return(unlist(lapply(summ_list, function(summ) {
    summ$f_peak
  })))
}