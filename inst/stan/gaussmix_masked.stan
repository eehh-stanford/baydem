functions {

  int bin_search(real x, int min_val, int max_val){
    // This assumes that min_val >= 0 is the minimum integer in range,
    //  max_val > min_val,
    // and that x has already been rounded.
    //  It should find the integer equivalent to x.
    int range = (max_val - min_val+1)/2; // We add 1 to make sure that truncation doesn't exclude a number
    int mid_pt = min_val + range;
    int out;
    while(range > 0) {
      if(x == mid_pt){
        out = mid_pt;
        range = 0;
      } else {
        // figure out if range == 1
        range =  (range+1)/2;
        mid_pt = x > mid_pt ? mid_pt + range: mid_pt - range;
        }
    }
    return out;
  }

  vector calc_logh(int N, int K, int G,
                       vector pi, vector s, vector mu, vector tau, real dtau,
                       vector stacked_log_M,
                       vector subset_length,
                       int subset_length_min, int subset_length_max,
                       vector M_offset,
                       int M_offset_min, int M_offset_max,
                       vector tau_offset,
                       int tau_offset_min, int tau_offset_max) {
    vector[N] logh; // Log-likelihood
    vector[G] logf; // log density
    vector[G] f; // density
    vector[K] logpi = log(pi); // log(pi), used for the log_sum_exp calculation
    real log_norm_factor;

    // Calculate the density as well as the normalization factor for the density
    logf = rep_vector(0,G);
    for (g in 1:G) {
      vector[K] lse = logpi;
      for (k in 1:K) {
        lse[k] += normal_lpdf(tau[g]|mu[k],s[k]);
      }
      logf[g] = log_sum_exp(lse);
    }
    f = exp(logf);
    //real log_norm_factor = log(sum(f)) + log(dtau);
    log_norm_factor = dtau;

    // Iterate over observations to calculate the log-likelihood, logh
    for (n in 1:N) {
        // Grid length for this sample
        int G_n = bin_search(subset_length[n], subset_length_min, subset_length_max);
        int M_offset_n = bin_search(M_offset[n], M_offset_min, M_offset_max);
        int tau_offset_n = bin_search(tau_offset[n], tau_offset_min, tau_offset_max);

        vector[G_n] lse = rep_vector(0, G_n); // Initialize the log-sum-exp vector
        // Iterate over grid points for this observation
        for (g in 1:G_n) {
          lse[g] = logf[tau_offset_n + g] + stacked_log_M[M_offset_n + g] - log_norm_factor;
        }
        logh[n] = log_sum_exp(lse);
    }
    return logh;
  }
}


data {
  int N; // Number of observations
  int G; // Number of grid pints
  int K; // Number of mixtures
  int stacked_log_M_length; // Length of stacked log M vector
  vector[stacked_log_M_length] stacked_log_M; // Stacked measurement matrix
  vector[G] tau; // Full calendar dates grid
  vector[N] subset_length; // Length of subset for each observation
  int subset_length_min; // Minimum value of subset_length
  int subset_length_max; // Maximum value of subset_length
  vector[N] M_offset; // Offset in stacked_log_M for each observation
  int M_offset_min; // Minimum M_offset value
  int M_offset_max; // Maximum M_offset value
  vector[N] tau_offset; // Offset in tau for each observation
  int tau_offset_min; // Minimum value of tau_offset
  int tau_offset_max; // Maximum value of tau_offset
  int tau_min; // Minimum calendar date
  int tau_max; // Maximum calendar date
  real dtau; // Grid spacing (not error checked for consistency with tau)
  real<lower=0> alpha_d; // Concentration parameter for pi prior
  real<lower=0> alpha_s; // shape parameter of gamma distribution for s
  real<lower=0> alpha_r; // rate parameter of gamma distribution for s
}

parameters {
  simplex[K] pi;
  ordered[K] mu;
  vector<lower=0>[K] s;
}

transformed parameters {
  vector[N] logh;
  logh =  calc_logh(N, K, G,
                    pi, s, mu, tau, dtau,
                    stacked_log_M,
                    subset_length,
                    subset_length_min, subset_length_max,
                    M_offset,
                    M_offset_min, M_offset_max,
                    tau_offset,
                    tau_offset_min, tau_offset_max);
}

model {
  pi ~ dirichlet(rep_vector(alpha_d,K));
  s ~ gamma(alpha_s,alpha_r);
  mu ~ uniform(tau_min,tau_max);
  target += sum(logh);
}
