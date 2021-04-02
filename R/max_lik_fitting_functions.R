#' Fit a truncated Gaussian mixture to radiocarbon data using the Hooke-Jeeves
#' algorithm with multiple restarts
#'
#' The radiocarbon data are specified by phi_m, the vector of fraction modern
#' values, and sig_m, the vector of associated standard deviations (both indexed
#' by i). The likelihood involves an integral over calendar dates that is
#' calculated using a Riemann integral at the points in the "sampling grid",
#' tau, which ranges from tau_min to tau_max and has a spacing of dtau (tau is
#' indexed by the variable g).
#'
#' Constraints on parameters: The parameter vector, theta (th), has the ordering
#' th = (pi_1,...,pi_K,mu_1,..mu_K,s_1,...s_K), where K is the number of
#' mixture components and pi_k / mu_k / s_k are, respectively, the weighting /
#' mean / standard deviation of the k-th mixture. Optimization is actually done
#' using a "reduced" vector from which pi_1 has been removed, th = (pi_2,...),
#' since the mixture weights must add to 1. Constraints are placed on each
#' variable in the parameter vector. All the weightings must lie between 0 and
#' 1, including pi_1 (even though it is not in th_reduced). The means must lie
#' between tau_min and tau_max. Since a numerical integral with a defined
#' spacing is used to calculate the likelihood, the standard deviations of the
#' mixture components cannot be too small, and are constrained to be greater
#' than 10*dtau; they are also constrained to be less than tau_max-tau_min. In
#' summary:
#'
#'       0 < pi_k < 1
#' tau_min < mu_k < tau_max
#' 10*dtau < s_k  < tau_max-tau_min.
#'
#' The bounded Hooke-Jeeves algorithm in the dfoptim R package (dfoptim::hjkb)
#' is used to do the optimization. The Hooke-Jeeves algorithm is perhaps the
#' best known varient of derivative free pattern-searches. Using multiple
#' restarts helps identify the global optimum (rather than a local optimum).
#'
#' To ensure reproducibility, the input_seed can be specified in one of two
#' ways. (1) input_seed is a single integer; if so, set the base_seed equal to
#' equal to input_seed (see below for how base_seed is used). (2) A vector of
#' integers of length 1 + num_restarts is specified; no base seed is needed if a
#' vector is provided. If no input_seed is provided (input_seed = NA), the base
#' seed is randomly generated and stored. (1) If a base_seed is used, the
#' command set.seed(base_seed) is run, then a vector of integers is drawn to
#' populate seed_vect. (2) If a vector is given for input_seed, seed_vect is set
#' equal to it. The first entry in seed_vect is used with set.seed() prior to
#' drawing the starting parameter vectors for the restarts. The subsequent
#' entries are used with set.seed() prior to calling the bounded Hooke-Jeeves
#' optimization algorithm for each restart. These additional seeds are necessary
#' to ensure reproducibility is achieved when multiple cores are used to do the
#' optimizations in a parallel for loop.
#'
#' @param K The number of mixture components to use for the fit
#' @param phi_m A vector of fraction moderns
#' @param sig_m A vector of standard deviations for phi_m
#' @param tau_min The minimum calendar date (AD) for the sampling grid
#' @param tau_max The maximum calendar date (AD) for the sampling grid
#' @param dtau The spacing of the sampling grid
#' @param calib_df Calibration curve (see load_calib_curve)
#' @param num_restarts Number of restarts to use (default: 100)
#' @param maxfeval Maximum number of function evaluations for each restart
#'   (default: (K-1)*10000)
#' @param num_cores The number of cores to use for parallelization (default: NA,
#'   do not parallelize)
#' @param input_seed Either a single integer or a vector of integers of length
#'   1 + num_restarts that is used to ensure reproducilibity (default: NA, a
#'   the base seed will be generated and saved).
#'
#' @return A list containing the best fit parameter vector (th), the
#'   corresponding value of the negative log-likelihood (neg_log_lik), tau (the
#'   sampling grid), f (the density of the best fit paramater vector evaluated
#'   at the points in the sampling grid), bic (the Bayesian information
#'   Criterion value), aic (the Akaike information criterion), hjkb_fit_list (
#'   a list containing the fits for each restart), input_seed (the value of
#'   input parameter input_seed), base_seed (the value of the base_seed, which
#'   is NA if input_seed is a vector), and seed_vect (the vector of integers
#'   of length 1 + num_restarts; see the function Description).
#'
#' @import doParallel
#' @import foreach
#' @export

fit_trunc_gauss_mix <- function(K,
                                phi_m,
                                sig_m,
                                tau_min,
                                tau_max,
                                dtau,
                                calib_df,
                                num_restarts=100,
                                maxfeval=(K-1)*10000,
                                num_cores=NA,
                                input_seed=NA) {

  # Create the sampling grid, tau
  tau <- seq(tau_min,tau_max,dtau)

  # Calculate the measurement matrix, M
  M <- calc_meas_matrix(tau, phi_m, sig_m, calib_df)

  # A total of 1 + num_restarts seeds are need to ensure reproducibility with
  # the parallel procesing: One seed is needed to randomly draw the restart
  # parameter vector and a seed is needed for the hjkb optimization for each
  # restart. Therefore, the input_seed should be one of three things:
  #
  # (1) NA              No input seed set. A base seed is randomly chosen and
  #                     flow proceeds as with option (2)
  # (2) Single integer  A single integer is provided as the seed. Use this to
  #                     set the random number seed generator and sample a vector
  #                     of 1 + num_restarts integers; flow then proceeds as with
  #                     option (3).
  # (3) Vector of       A vector of integers of length 1 + num_restarts is
  #     integers        provided. Use these directly.

  if (length(input_seed) == 1) {
    if (is.na(input_seed)) {
      base_seed <- sample.int(1000000,1)
    } else {
      base_seed <- input_seed
    }
    set.seed(base_seed)
    seed_vect <- sample.int(1000000,1+num_restarts)
  } else if(length(input_seed) == 1 + num_restarts) {
    base_seed <- NA
    seed_vect <- input_seed
  } else {
    stop(paste0("input_seed must be NA, a single integer, or a vector of ",
                "integers with length 1 + num_restarts"))
  }

  # Iterate over restarts. To generate the starting vectors, assume the
  # following:
  # (1) pi, the mixture proportions, are drawn from a Dirichlet distribution.
  # (2) mu, the means, are restricted to the interval tau_min to tau_max.
  # (3) s, the standard deviations, are restricted to the interval 10*dtau to
  #     tau_max-tau_min
  set.seed(seed_vect[1])
  TH0_reduced <- matrix(NA,3*K-1,num_restarts)
  for(m in 1:num_restarts) {
    th0 <- c(MCMCpack::rdirichlet(1,rep(1,K)),
             sort(runif(K,tau_min,tau_max)),
             runif(K,10*dtau,tau_max-tau_min))
    th0_reduced <- th0[-1]
    TH0_reduced[,m] <- th0_reduced
  }

  # Create the lower and upper bounds for optimization
  lower_bounds <- c(rep(0,K-1),rep(tau_min,K),rep(dtau*10,K))
  upper_bounds <- c(rep(1,K-1),rep(tau_max,K),rep(tau_max-tau_min,K))

  # Re-wrap the negative log-likelihood (the objective function) to enforce
  # the bounds on the parameter vector.
  neg_log_lik_rewrap <- function(th_reduced,M,tau,
                                 lower_bounds,upper_bounds) {
    if (any(th_reduced < lower_bounds)) {
      return(Inf)
    }

    if (any(th_reduced > upper_bounds)) {
      return(Inf)
    }

    return(calc_trunc_gauss_mix_neg_log_lik(th_reduced,M,tau))
  }

  # Create the control list variable to set the maximum number of function
  # evaluations
  control_list <- list(maxfeval=maxfeval)

  # Solve the problem using either a conventional for loop (if num_cores is NA)
  # or a parallel for loop. For the parallel for loops, it faster to preload
  # libraries and avoid using namespacing with :: (as is called for by some
  # style guides); hence, @import tags have been added to the function
  # documentation.
  if(is.na(num_cores)) {
    # Use a conventional for loop to do the fits

    # A list in which to store the fits for the restarts
    hjkb_fit_list <- list()
    for(m in 1:num_restarts) {
      set.seed(seed_vect[1 + m])
      hjkb_fit_list[[m]] <- dfoptim::hjkb(TH0_reduced[,m],
                                          neg_log_lik_rewrap,
                                          M=M,
                                          tau=tau,
                                          lower_bounds=lower_bounds,
                                          upper_bounds=upper_bounds,
                                          control=control_list)
    }
  } else {
    # Use a parallel for loop to do the fits
   doParallel::registerDoParallel(num_cores)
    hjkb_fit_list <-
      foreach(m=1:num_restarts,.packages=c('baydem')) %dopar% {
        set.seed(seed_vect[1 + m])
        dfoptim::hjkb(TH0_reduced[,m],
                      neg_log_lik_rewrap,
                      M=M,
                      tau=tau,
                      lower_bounds=lower_bounds,
                      upper_bounds=upper_bounds,
                      control=control_list)
      }
  }
  doParallel::stopImplicitCluster()

  # Identify the best solution across restarts and get the corresponding best
  # parameter vector and negative log-likelihood value.
  m_best <- which.min(unlist(lapply(
    hjkb_fit_list,function(hjkb_fit){hjkb_fit$value})))
  th_best <- hjkb_fit_list[[m_best]]$par
  th_best <- c(1-sum(th_best[1:(K-1)]),th_best)
  neg_log_lik_best <- hjkb_fit_list[[m_best]]$value


  # Calculate the probability density for the best fit parameter
  f  <- calc_gauss_mix_pdf(th_best, tau, tau_min, tau_max)

  # Calculate the Bayesian information criterion (BIC) and Akaike information
  # criterion (AIC).
  N <- length(phi_m)
  bic <- log(N)*(3*K-1) + 2*neg_log_lik_best
  aic <-      2*(3*K-1) + 2*neg_log_lik_best

  # Return a list with:
  # th          The best fit paramater vector
  # neg_log_lik The best value of the objective function
  # tau         The sampling grid
  # f           The probability density evaluted at locations in tau
  # bic         The Bayesian information criterion
  # aic         The Akaike information criterion
  return(list(th=th_best,neg_log_lik=neg_log_lik_best,tau=tau,f=f,
              bic=bic,aic=aic,hjkb_fit_list=hjkb_fit_list,
              input_seed=input_seed,base_seed=base_seed,seed_vect=seed_vect))
}

#' Calculate the negative log-likelihood of a truncuated Gaussian mixture fit
#'
#' The input parameter vector is the parameter vector without the first mixture
#' weight. The minimum and maximum calendar dates are assumed to be the minumum
#' and maximum values of the sampling grid, tau.
#'
#' @param th_reduced The paramater vector without the first mixture weight
#' @param M The measurement matrix
#' @param tau The sampling grid
#' @param s_min A minimum value for the mixture standard deviations
#'   (default: 0)
#'
#' @return The negative log-likelihood
#'
#' @export
calc_trunc_gauss_mix_neg_log_lik <- function(th_reduced,M,tau,s_min=0) {
  # If the parameter vector is invalid, return infinity
  if (!is_th_reduced_valid(th_reduced,s_min)) {
    return(Inf)
  }

  tau_min <- min(tau)
  tau_max <- max(tau)

  # Add an undersore to pi, the mixture proportions, since pi is 3.14... in R
  K <- (length(th_reduced) + 1)/3
  pi_ <- th_reduced[1:(K-1)]

  th <- c(1-sum(pi_),th_reduced)
  # Calculate v, the vector of densities
  v <- calc_gauss_mix_pdf(th, tau, tau_min, tau_max)
  h <- M %*% v
  return(-sum(log(h)))
}

#' A helper function to determine if the reduced form parameter vector is valid
#'
#' The paramater vector is invalid if:
#'
#' (a) Any mixture weight does not lie between 0 and 1
#' (b) mu is not ordered from small to large values
#' (c) A standard deviation is lower than s_min
#'
#' Neither bounds on the values of the means or the upper values of the standard
#' deviations are checked, but they are enforced in the optimization in
#' fit_trunc_gauss_mix.
#'
#' @param th_reduced The paramater vector without the first mixture weight
#' @param s_min A minimum value for the mixture standard deviations
#'   (default: 0)
#'
#' @return TRUE or FALSE
#'
#' @export
is_th_reduced_valid <- function(th_reduced,s_min=0) {
  # TODO: Consider adding a check on the length of th_reduced
  K <- (length(th_reduced) + 1)/3

  # Add an undersore to pi, the mixture proportions, since pi is 3.14149... in R
  pi_ <- th_reduced[1:(K-1)]
  pi_ <- c(sum(pi_),pi_) # mixture proportions must sum to 1

  # Mixture propotions must all lie between 0 and 1
  if (any(pi_ < 0)) {
    return(FALSE)
  }

  if (any(pi_ > 1)) {
    return(FALSE)
  }

  # The means must be ordered
  mu <- th_reduced[K:(2*K-1)]
  if (is.unsorted(mu)) {
    return(FALSE)
  }

  # The standard deviations must be strictly positive (or greater than s_min)
  s <- th_reduced[(2*K):(3*K-1)]
  if(any(s <= s_min)) {
    return(FALSE)
  }
  return(TRUE)
}
