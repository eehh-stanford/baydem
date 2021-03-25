#' Use parallel tempering to fit a truncated Gaussian mixture to radiocarbon
#' data
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
#' The parallel tempering algorithm in the enneal R package is used to do the
#' optimization. This ensures that the space of possible curves is well-explored
#' and possible local optima are avoided.
#'
#' @param K The number of mixture components to use for the fit
#' @param phi_m A vector of fraction moderns
#' @param sig_m A vector of standard deviations for phi_m
#' @param tau_min The minimum calendar date (AD) for the sampling grid
#' @param tau_max The maximum calendar date (AD) for the sampling grid
#' @param dtau The spacing of the sampling grid
#' @param calib_df Calibration curve (see bd_load_calib_curve)
#' @param samps_per_cyc Number of samples per cycle for the tempering
#'  (default:10)
#' @param num_cyc Number of cycles for the tempering (default:10)
#' @param num_cores The number of cores to use for parallelization (default: NA,
#'   do not parallelize)
#'
#' @return A list containing the best fit parameter vector (th), the
#'   corresponding value of the negative log-likelihood (neg_log_lik), tau (the
#'   sampling grid), and f (the density of the best fit paramater vector
#'   evaluated at the points in the sampling grid)
#'
#' @import doParallel
#' @import foreach
#' @export

temper_trunc_gauss_mix <- function(K,
                                   phi_m,
                                   sig_m,
                                   tau_min,
                                   tau_max,
                                   dtau,
                                   calib_df,
                                   samps_per_cyc=10,
                                   num_cyc=1000,
                                   num_cores=NA) {

  # Create the sampling grid, tau
  tau <- seq(tau_min,tau_max,dtau)

  # Calculate the measurement matrix, M
  M <- bd_calc_meas_matrix(tau, phi_m, sig_m, calib_df)

  # Iterate over restarts. To generate the starting vectors, assume the
  # following:
  # (1) pi, the mixture proportions, are drawn from a Dirichlet distribution.
  # (2) mu, the means, are restricted to the interval tau_min to tau_max.
  # (3) s, the standard deviations, are restricted to the interval 10*dtau to
  #     tau_max-tau_min

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

  # For the parallel tempering, use a temperature vector with temperatures
  # between sqrt(10) and 1/sqrt(10)
  temp_vect <- rev(10^seq(-0.5,0.5,len=10))

  # For the largest temperature, use the following for the standard deviation of
  # the proposal distribution:
  # pi_k  0.5
  # mu_k  (tau_max-tau_min)/2
  #  s_k  (tau_max-tau_min)/4
  tau_range <- tau_max - tau_min
  scale_vect0 <- c(rep(0.5,K-1),rep(tau_range/2,K),rep(tau_range/4,K))

  # For other temperatures, rescale the standard deviation used for the proposal
  # distribution such that the the smallest temperature has a scale 0.01 that
  # of the largest.
  scale_matrix <- matrix(NA,3*K-1,length(temp_vect))
  rescale_vect <- seq(1,.01,len=length(temp_vect))
  for(tt in 1:length(temp_vect)) {
    scale_matrix[,tt] <- scale_vect0 * rescale_vect[tt]
  }

  th0 <- c(MCMCpack::rdirichlet(1,rep(1,K)),
           sort(runif(K,tau_min,tau_max)),
           runif(K,10*dtau,tau_max-tau_min))
  th0_reduced <- th0[-1]

  if (is.na(num_cores)) {
    # If num_cores is NA, iterate over restarts using a regular for loop
    # A list of fits
    tempering <- enneal::par_temper(th0_reduced,
                                    neg_log_lik_rewrap,
                                    samps_per_cyc=samps_per_cyc,
                                    temp_vect=temp_vect,
                                    prop_scale=scale_matrix,
                                    num_cyc=num_cyc,
                                    M=M,
                                    tau=tau,
                                    lower_bounds=lower_bounds,
                                    upper_bounds=upper_bounds)
  } else {
    # For speed, do not use :: to call functions. This means that doParallel
    # (for registerDoParallel) and foreach (for %dopar%) must have been
    # pre-loaded; hence, @import statements have been added to the function
    # documentation.
    registerDoParallel(num_cores)
    # Must explicitly export functions so that dopar works on Windows
    #function_exports <- c("calc_trunc_gauss_mix_neg_log_lik",
    #                      "is_th_reduced_valid")
    tempering <- enneal::par_temper(th0_reduced,
                                    neg_log_lik_rewrap,
                                    samps_per_cyc=samps_per_cyc,
                                    temp_vect=temp_vect,
                                    prop_scale=scale_matrix,
                                    num_cyc=num_cyc,
                                    num_cores=num_cores,
                                    M=M,
                                    tau=tau,
                                    lower_bounds=lower_bounds,
                                    upper_bounds=upper_bounds)
  }

  ind_best <- which.min(unlist(lapply(
    tempering$chains,function(fit){fit$eta_best})))
  th_best  <- tempering$chains[[ind_best]]$theta_best
  th_best <- c(1-sum(th_best[1:(K-1)]),th_best)
  eta_best <- tempering$chains[[ind_best]]$eta_best

  # Calculate the probability density for the best fit parameter
  f  <- bd_calc_gauss_mix_pdf(th_best,tau,tau_min,tau_max)

  # Calculate the Bayesian information criterion (BIC) and Akaike information
  # criterion (AIC). General wisdom holds that BIC is preferred if the true
  # model is in the set of models tested and AIC is preferred if it may not be.
  N <- length(phi_m)
  bic <- log(N)*(3*K-1) + 2*eta_best
  aic <-      2*(3*K-1) + 2*eta_best

  # Return a list with:
  # th          The best fit paramater vector
  # neg_log_lik The best value of the objective function
  # tau         The sampling grid
  # f           The probability density evaluted at locations in tau
  # bic         The Bayesian information criterion
  # aic         The Akaike information criterion
  return(list(th=th_best,neg_log_lik=eta_best,tau=tau,f=f,bic=bic,aic=aic))
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
#' @param sig_min A minimum value for the mixture standard deviations
#'   (default: 0)
#'
#' @return The negative log-likelihood
#'
#' @export
calc_trunc_gauss_mix_neg_log_lik <- function(th_reduced,M,tau,sig_min=0) {
  # If the parameter vector is invalid, return infinity
  if (!is_th_reduced_valid(th_reduced,sig_min)) {
    return(Inf)
  }

  tau_min <- min(tau)
  tau_max <- max(tau)

  # Add an undersore to pi, the mixture proportions, since pi is 3.14... in R
  K <- (length(th_reduced) + 1)/3
  pi_ <- th_reduced[1:(K-1)]

  th <- c(1-sum(pi_),th_reduced)
  # Calculate v, the vector of densities
  v <- bd_calc_gauss_mix_pdf(th,tau,tau_min,tau_max)
  h <- M %*% v
  return(-sum(log(h)))
}

#' A helper function to determine if the reduced form parameter vector is valid
#'
#' The paramater vector is invalid if:
#'
#' (a) Any mixture weight does not lie between 0 and 1
#' (b) mu is not ordered from small to large values
#' (c) A standard deviation is lower than sig_min
#'
#' Neither bounds on the values of the means or the upper values of the standard
#' deviations are checked, but they are enforced in the optimization in
#' bd_fit_trunc_gauss_mix.
#'
#' @param th_reduced The paramater vector without the first mixture weight
#' @param sig_min A minimum value for the mixture standard deviations
#'   (default: 0)
#'
#' @return TRUE or FALSE
#'
#' @export
is_th_reduced_valid <- function(th_reduced,sig_min=0) {
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

  # The standard deviations must be strictly positive (or greater than sig_min)
  sig <- th_reduced[(2*K):(3*K-1)]
  if(any(sig <= sig_min)) {
    return(FALSE)
  }
  return(TRUE)
}