#' @title Sample from the posterior of a density model for a set of radiocarbon
#' measurements
#'
#' @description This is the core function that implements the Bayesian
#' inference. Currently, the only supported density model is a truncated
#' Gaussian mixture.
#'
#' @details The input consists of the radiocarbon measurements (rc_meas), the
#' density model (density_model), hyperparameters that specify the prior
#' distribution as well as the spacing of the sampling grid (hp), and optional
#' control variables. stan is called via the rstan package to sample from the
#' posterior. The output is the variable bayesian_soln, which is a list
#' TODO: continue editing documentation from here once edits are complete
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
#' @return bayesian_soln, a list with three fields: prob (the input), fit (the result of
#'   the call to stan), and control (the control parameters used)
#'
#' @export
do_bayesian_inference <- function(rc_meas,
                                  density_model,
                                  hp,
                                  calib_df,
                                  th0=NA,
                                  init_seed=NA,
                                  stan_seed=NA,
                                  control=list()) {

  # Save the the optional inputs, which are stored in the return value
  optional_inputs <- list(th0=th0,
                          init_seed=init_seed,
                          stan_seed=stan_seed,
                          control=control)

  have_th0 <- !all(is.na(th0))
  have_init_seed <- !is.na(init_seed)

  # Raise an error if both th0 and init_seed are provided
  if(have_th0 && have_init_seed) {
    stop("init_seed should not be provided if th0 is provided")
  }

  # If necessary, draw the initialization seed
  if (!have_th0 && !have_init_seed) {
    init_seed <- sample.int(1000000,1)
  }

  # If necessary, draw the stan seed
  if (is.na(stan_seed)) {
    stan_seed <- sample.int(1000000,1)
  }

  # Unpack and/or define the control parameters
  have_num_chains      <- "num_chains"           %in% control
  have_samps_per_chain <- "have_samps_per_chain" %in% control
  have_warmup          <- "have_warmup"          %in% control
  have_stan_control    <- "have_stan_control"    %in% control

  if (have_num_chains) {
    num_chains <- control$num_chains
  } else {
    num_chains <- 4
  }

  if (have_samps_per_chain) {
    samps_per_chain <- control$samps_per_chain
  } else {
    samps_per_chain <- 2000
  }

  if (have_warmup) {
    warmup <- control$warmup
  } else {
    warmup <- floor(samps_per_chain / 2)
  }

  if (have_stan_control) {
    stan_control <- control$stan_control
  } else {
    stan_control <- NA
  }

  control_final <- list(
    num_chains = num_chains,
    samps_per_chain = samps_per_chain,
    warmup = warmup,
    stan_control = stan_control
  )

  if (density_model$type == "trunc_gauss_mix") {
    # Stan needs all the inputs and hyperparameters as variables in R's
    # workspace
    tau_min <- density_model$tau_min
    tau_max <- density_model$tau_max
    tau <- seq(tau_min, tau_max, by = hp$dtau)
    M <- calc_meas_matrix(tau, rc_meas$phi_m, rc_meas$sig_m, calib_df)

    if (all(is.na(th0))) {
      # Then an initialization vector has not been provided. Do a maximum
      # likelihood fit to get a starting parameter. Use init_seed for this to
      # allow reproducibility.
      dtau <- hp$dtau
      K <- density_model$K
      max_lik_fit <- fit_trunc_gauss_mix(K,
                                         rc_meas$phi_m,
                                         rc_meas$sig_m,
                                         tau_min,
                                         tau_max,
                                         dtau,
                                         calib_df,
                                         input_seed=init_seed)
      th0 <- max_lik_fit$th
    }

    # stan expects the initial parameter to be a named list

    init0 <- list()
    init0$pi <- th0[  0 + (1:K)]
    init0$mu <- th0[  K + (1:K)]
    init0$s  <- th0[2*K + (1:K)]

    init_list <- list()
    for (cc in 1:num_chains) {
      init_list[[cc]] <- init0
    }

    Mt <- t(M)
    N <- dim(M)[1]
    G <- dim(M)[2]
    alpha_s <- hp$alpha_s
    alpha_r <- hp$alpha_r
    alpha_d <- hp$alpha_d
    K <- density_model$K
    file_path <- system.file("stan/gaussmix.stan",
      package = "baydem"
    )
    options(mc.cores = parallel::detectCores())
    # There are four possible calls depending on whether have_stan_control is
    # TRUE and have_stan_seed is TRUE

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
  } else {
    stop(paste("Unrecognized fit type:", density_model$fit_type))
  }

  bayesian_soln <- list(fit=fit,
                        control=control_final,
                        init_seed=init_seed,
                        stan_seed=stan_seed,
                        optional_inputs=optional_inputs)
  class(bayesian_soln) <- "bd_bayesian_soln"

  # If a save file was input, save the result to file. This is especially
  # useful for parallel batch runs
  return(bayesian_soln)
}

#' @title Calculate some key summary measures using the result of a call to
#' \code{do_bayesian_inference}
#'
#' @description \code{do_bayesian_inference} calls Stan to do Bayesian inference
#' by generating a sample of parameters from the posterior of theta (or th).
#' do_bayesian_inference analyzes the result of that inference. Notably,
#' it calculates the quantiles of the density function and the growth rate.
#'
#' @details `bayesian_soln` is the result of a call to \code{do_bayesian_inference}.
#' It contains both the resulting samples and the parameters used in the
#' inference, such as the hyperparameters (see \code{do_bayseian_inference} for
#' further details). The primary thing \code{summarize_bayesian_inference} does
#' is calculate quantiles of both the parameterized density and growth rate. For
#' example, for a calendar date tau_g each sample yields a density and growth
#' rate. The quantile is the value of the density or growth rate such that a
#' given proportion of samples are smaller than that value. The probabilities
#' used to calculate these quantiles are `probs = c(lev, 0.5, 1-lev)`, where
#' `lev` is the level (0.025 by default, so that 95% of the observations lie
#' between the first and last quantile bands).
#'
#' In addition, \code{summarize_bayesian_inference} identifies calendar dates
#' for which the growth rate quantiles defined by `lev` and `1 - lev` do not
#' contain zero. This indicates significant positive or negative growth for the
#' density curve. The output vector `growth_state` codes calendar dates by
#' growth state as 'negative', 'zero', and 'positive'. For the Gaussian mixture
#' parameterization of the density, the rate is not typically meaningful near
#' the calendar date boundaries where it increases linearly as the calendar date
#' goes to positive or negative infinity. The parameter `rate_prop` provides
#' control on how calendar dates are classified by growth rate near these
#' boundaries. In particular, the calendar dates with a cumulative density (50%
#' quantile) below `rate_prop` (for the lower boundary) or above `1 - rate_prop`
#' (for the upper boundary) are classified as 'missing' in `growth_state`. By
#' default, `rate_prop` is NA and no calendar dates are classified as missing.
#'
#' By default, a summary is done for each sample by calling summarize_sample.
#' This is not done if doSummary is FALSE.
#'
#' @param bayesian_soln The solution, a list-like object of class bd_bayesian_soln (see
#'   \code{do_bayesian_inference})
#' @param tau (optional) The calendar dates at which to evaluate densities. If
#'   tau is not input, tau is built from the hyperparameters.
#' @param th_sim (optional) The known parameters used to create simulation data
#' @param lev (default: 0.025) The level to use for the quantile bands
#' @param rate_prop (optional) The cumulative density needed to define rate
#'   growth bands
#' @param do_sample_summaries (default: `TRUE`) Whether to calculate some
#'   summary summarize for each sampled curve.
#'
#' @return A list with information on the quantiles of the density function and
#'   growth rate (and sample summaries)
#'
#' @export

summarize_bayesian_inference <- function(bayesian_soln,
                                         rc_meas,
                                         density_model,
                                         calib_df,
                                         dtau=5,
                                         th_sim = NA,
                                         lev = 0.025,
                                         rate_prop = NA,
                                         do_sample_summaries = T) {

  if (density_model$type != "trunc_gauss_mix") {
    stop(paste0("Currently, only truncated Gaussian mixtures are supported ",
                "for the density model"))
  }
  tau_min <- density_model$tau_min
  tau_max <- density_model$tau_max
  tau <- seq(tau_min,tau_max,dtau)

  # The probabilities to use to calculate the quantiles
  probs <- c(lev, 0.5, 1 - lev)

  # Extract the samples of theta in the variable TH to create TH, a matrix with
  # dimensions N x Number of Parameters, where N is the number of Bayesian
  # samples
  TH <- extract_param(bayesian_soln$fit)

  num_samp <- nrow(TH)

  # Calculate the pdf matrix, which is the probability density for the input
  # density_model evaluated for each sample and at each grid point in the vector
  # tau. f_mat has dimensions N x G, where G is the number of grid points
  # (length of the vector tau).
  f_mat <- calc_gauss_mix_pdf_mat(TH,
                                  tau,
                                  tau_min=tau_min,
                                  tau_max=tau_max)

  # Calculate the rate for each sample and grid point. The rate is f' / f, where
  # f is the probability density and f' is the derivative of f.
  rate_mat <- calc_gauss_mix_pdf_mat(TH,
                                     tau,
                                     tau_min=tau_min,
                                     tau_max=tau_max,
                                     type="rate")

  # Calculate the quantiles of the probability density
  Qdens <- calc_quantiles(f_mat, probs)

  # Extract the 50% quantile of the probability density (which is not itself a
  # probability density; that is, it does not integrate to 1). By construction,
  # the 50% quantile is the second row of the matrix Qdens.
  f50 <- Qdens[2, ]

  # Restrict to indices with enough probability mass (if necessary)
  if (!is.na(rate_prop)) {
    rate_ind <- which(cumsum(f50 * dtau) > rate_prop &
                        rev(cumsum(rev(f50) * dtau)) > rate_prop)
  } else {
    rate_ind <- 1:length(f50)
  }

  # Identify regions with growth rates that differ from zero per the input
  # quantile level (lev).
  Qrate <- calc_quantiles(rate_mat[, rate_ind], probs)

  # Create growth_state0, a base vector summarizing the growth rate states.
  # growth_state0 has the value "negative" for significant negative growth,
  # "positive" for significant positive growth, and "zero" otherwise.
  # growth_state is identical to growth_state0, except that it has the value
  # "missing" for observations without enough probability mass per the rate_prop
  # condition.
  growth_state0 <- rep("zero", length(rate_ind))
  growth_state0[Qrate[2, ] > 0 & Qrate[1, ] > 0] <- "positive"
  growth_state0[Qrate[2, ] < 0 & Qrate[3, ] < 0] <- "negative"
  growth_state <- rep("missing", length(tau))
  growth_state[rate_ind] <- growth_state0


  # Calculate the measurement matrix
  M <- calc_meas_matrix(tau,
                        rc_meas$phi_m,
                        rc_meas$sig_m,
                        calib_df)

  # TODO: consider adding a separate function to calculate the SPD
  # Normalize the measurement matrix by row
  M <- M / replicate(length(tau),rowSums(M)*dtau)
  # Calculate summed probability density (SPD) vector
  f_spdf <- colMeans(M)

  # Create the output list
  out <- list(
    tau = tau,
    f_spdf = f_spdf,
    Qdens = Qdens,
    Qrate = Qrate,
    probs = probs,
    rate_prop = rate_prop,
    rate_ind = rate_ind,
    growth_state = growth_state
  )

  class(out) <- "bd_bayesian_summary"

  if (do_sample_summaries) {
    summ_list <- list()
    for (n in 1:num_samp) {
      th <- TH[n, ]
      summ_list[[n]] <- summarize_trunc_gauss_mix_sample(th,
                                                         tau_min,
                                                         tau_max)
    }
    out$summ_list <- summ_list
  }

  have_sim <- !all(is.na(th_sim))
  if (have_sim) {
    f_sim <- calc_gauss_mix_pdf(th_sim,
                                tau,
                                tau_min=tau_min,
                                tau_max=tau_max)
    rate_sim <- calc_gauss_mix_pdf(th_sim,
                                   tau,
                                   tau_min=tau_min,
                                   tau_max=tau_max,
                                   type = "rate")
    out$f_sim <- f_sim
    out$rate_sim <- rate_sim
  }
  return(out)
}

#' @title Extract the Bayesian samples for a Gaussian mixture model generated by
#' do_bayesian_inference
#'
#' @description
#' The input fit is the result of a call to stan by
#' \code{do_bayesian_inference}, of class stanfit. Return a matrix TH with
#' dimensions S x (3*K), where S is the number of samples (across all chains,
#' and excluding warmup), and K is the number of mixtures. The final column of
#' as.matrix for class stanfit is the log-posterior, which must be removed.
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
#' @param bayesian_soln The result of a call to do_bayesian_inference
#' @param prop_change (Default 0.5) The relative decrease in density to use for
#'   the duration calculation
#' @param bayesian_summary (Optional) The result of a call to summarize_bayesian_inference.
#'   If not provided, it is calculated
#' @param peak_range (default: `c(tau_min, tau_max)`) peak_range can be given so
#'   that the peak density used is on the range peak_range
#'
#' @return A vector of "half-lives" (proportional change set by prop_change)
#'
#' @export
calc_half_life_from_peak <- function(bayesian_soln,
                                     density_model,
                                     rc_meas=list(),
                                     calib_df=list(),
                                     prop_change=0.5,
                                     bayesian_summary=NA,
                                     peak_range = NA) {

  if (density_model$type != "trunc_gauss_mix") {
    stop(paste0("Currently, only truncated Gaussian mixtures are supported ",
                "for the density model"))
  }

  tau_min <- density_model$tau_min
  tau_max <- density_model$tau_max

  TH <- extract_param(bayesian_soln$fit)
  N <- nrow(TH)

  if (all(is.na(bayesian_summary))) {
    if (length(rc_meas) == 0) {
      stop("rc_meas must be provided if bayesian_summary is not provided")
    }
    if (length(calib_df) == 0) {
      stop("calib_df must be provided if bayesian_summary is not provided")
    }
    bayesian_summary <- summarize_bayesian_inference(bayesian_soln,
                                                     rc_meas,
                                                     density_model,
                                                     calib_df)
  }

  summ_list <- bayesian_summary$summ_list

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
#' By default, this calculation is done for all the Bayesian samples in bayesian_soln,
#' which is the result of a call to do_bayesian_inference. Optionally, a subset
#' can be specified via the input ind, which should be a vector of integer
#' indices at which to do the calculation. To save computation if either spec1
#' or spec2 is 'peak', the result of a call to summarize_bayesian_inference for
#' which do_summary was TRUE can be input.
#'
#'
#' @param bayesian_soln The result of a call to do_bayesian_inference
#' @param spec1 The specification for the first density (see details)
#' @param spec2 The specification for the second density (see details)
#' @param ind (Optional) Indices at which to do the calculation. By default,
#'            all the samples in bayesian_summary are used.
#' @param bayesian_summary (Optional) The result of a call to summarize_bayesian_inference.
#'   This is only needed if either spec1 or spec2 is 'peak'
#'
#' @return A vector of relative densities (f_spec1 / f_spec2)
#'
#' @export
calc_relative_density <- function(bayesian_soln,
                                  density_model,
                                  spec1,
                                  spec2,
                                  rc_meas=list(),
                                  calib_df=list(),
                                  ind=NA,
                                  bayesian_summary=NA) {

  TH <- extract_param(bayesian_soln$fit)
  N <- nrow(TH)
  if (all(is.na(ind))) {
    ind <- 1:N
  }

  # Interpret and do error checking on inputs by calling the helper function
  # unpack_spec
  spec1 <- unpack_spec(spec1,density_model,T)
  spec2 <- unpack_spec(spec2,density_model,F)

  if (spec1$type == "peak" || spec2$type == "peak") {
    if (all(is.na(bayesian_summary))) {
      if (length(rc_meas) == 0) {
        stop("rc_meas must be provided if bayesian_summary is not provided")
      }
      if (length(calib_df) == 0) {
        stop("calib_df must be provided if bayesian_summary is not provided")
      }
      # If ind is not NA, the following line may involve un-utilized computation
      bayesian_summary <- summarize_bayesian_inference(bayesian_soln,
                                                       rc_meas,
                                                       density_model,
                                                       calib_df)
    }
    summ_list <- bayesian_summary$summ_list[ind]
  }

  # Calculate the density for spec1
  if (spec1$type == "point") {
    f1 <- calc_point_density(TH[ind, ],
                             density_model,
                             spec1$value)
  } else if (spec1$type == "range") {
    f1 <- calc_range_density(TH[ind, ],
                             density_model,
                             spec1$lower,
                             spec1$upper)
  } else if (spec1$type == "peak") {
    f1 <- calc_peak_density(summ_list)
  } else {
    # This should not happen, but throw an error regardless
    stop("Unsupported spec type")
  }

  # Calculate the density for spec2
  if (spec2$type == "point") {
    f2 <- calc_point_density(TH[ind, ],
                             density_model,
                             spec2$value)
  } else if (spec2$type == "range") {
    f2 <- calc_range_density(TH[ind, ],
                             density_model,
                             spec2$lower,
                             spec2$upper)
  } else if (spec2$type == "peak") {
    f2 <- calc_peak_density(summ_list)
  } else {
    # This should not happen, but throw an error regardless
    stop("Unsupported spec type")
  }

  return(f1 / f2)
}

# A helper function to unpack and do error checking on inputs spec1 / spec2
unpack_spec <- function(spec,density_model,is_one) {
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
      if (point < density_model$tau_min || density_model$tau_max < point) {
        stop(
          paste(s, "is a single date, but not in the range tau_min to tau_max"))
      }
      return(list(type = "point", value = point))
    } else if (length(spec) == 2) { # Numeric / length 2
      lower <- spec[1]
      if (lower < density_model$tau_min || density_model$tau_max < lower) {
        stop(paste(s,"is a date range, but lower value is not ",
                   "in the range tau_min to tau_max"))
      }
      upper <- spec[2]
      if (upper < density_model$tau_min || density_model$tau_max < upper) {
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
calc_point_density <- function(TH,density_model,t) {
  return(as.numeric(calc_gauss_mix_pdf_mat(TH,
                                           t,
                                           tau_min = density_model$tau_min,
                                           tau_max = density_model$tau_max)))
}


# A helper function to calculate the mean density over a range
calc_range_density <- function(TH,density_model,t_lo,t_hi) {
  f_lo <- as.numeric(calc_gauss_mix_pdf_mat(TH,
                                            t_lo,
                                            tau_min = density_model$tau_min,
                                            tau_max = density_model$tau_max,
                                            type = "cumulative"))
  f_hi <- as.numeric(calc_gauss_mix_pdf_mat(TH,
                                            t_hi,
                                            tau_min = density_model$tau_min,
                                            tau_max = density_model$tau_max,
                                            type = "cumulative"))
  return((f_hi - f_lo) / (t_hi - t_lo))
}


# A helper function to calculate the peak density
calc_peak_density <- function(summ_list) {
  return(unlist(lapply(summ_list, function(summ) {
    summ$f_peak
  })))
}