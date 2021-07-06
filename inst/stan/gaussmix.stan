data {
  int N; // Number of observations
  int G; // Number of grid pints
  int K; // Number of mixtures
  matrix[G,N] Mt; // Measurement matrix (transposed)
  vector[G] tau; // Calendar grid points for the measurement matrix calculation
  real dtau; // Grid spacing (not error checked for consistency with tau)
  real<lower=0> alpha_d; // Concentration parameter for pi prior
  real tau_min; // Lower calendar date. Same as min(tau)
  real tau_max; // Upper calendar date. Same as max(tau)
  real<lower=0> alpha_s; // shape parameter of gamma distribution for s
  real<lower=0> alpha_r; // rate parameter of gamma distribution for s
}

parameters {
  simplex[K] pi;
  ordered[K] mu;
  vector<lower=0>[K] s;
}

// Do calculations inside a function so that variables are not created in
// transformed_parameters that are exported. If this is not done, the return
// object (and save files) can be very large.
//functions {
//  calc_h() {
//
//    row_vector[N] h;
//    row_vector[N] logh;
//    row_vector[G] f;
//    row_vector[G] logf;
//    vector[K] logpi = log(pi);
//
//    logf = to_row_vector(rep_vector(0,G));
//    for (g in 1:G) {
//      vector[K] lse = logpi;
//      for (k in 1:K) {
//        lse[k] += normal_lpdf(tau[g]|mu[k],s[k]);
//      }
//      logf[g] = log_sum_exp(lse);
//    }
//    f = exp(logf);
//
//    // The following line truncates the distribution by normalizing the density
//    // to integrate to 1 on the interval of tau, which is assumed to run from
//    tau_min to tau_max.
//    f = f / sum(f) / dtau;
//    h = f * Mt;
//    return abs_diff / avg_scale;
//  }
//}
//
//transformed parameters {
//  row_vector[N] h;
//  row_vector[N] logh;
//  h = calc_h();
//  logh = log(h);
//}

transformed parameters {
  row_vector[N] h;
  row_vector[N] logh;
  row_vector[G] f;
  row_vector[G] logf;
  vector[K] logpi = log(pi);

  logf = to_row_vector(rep_vector(0,G));
  for (g in 1:G) {
    vector[K] lse = logpi;
    for (k in 1:K) {
      lse[k] += normal_lpdf(tau[g]|mu[k],s[k]);
    }
    logf[g] = log_sum_exp(lse);
  }
  f = exp(logf);

  // The following line truncates the distribution by normalizing the density to
  // integrate to 1 on the interval of tau, which is assumed to run from tau_min
  // to tau_max.
  f = f / sum(f) / dtau;
  h = f * Mt;
  logh = log(h);
}

model {
  pi ~ dirichlet(rep_vector(alpha_d,K));
  s ~ gamma(alpha_s,alpha_r);
  mu ~ uniform(tau_min,tau_max);
  target += sum(logh);
}

generated quantities {
  vector[N] log_post_vect;
  for (n in 1:N) {
    log_post_vect[n] = logh[n]
                       + dirichlet_lpdf( pi | rep_vector(alpha_d,K))
                       + gamma_lpdf( s | alpha_s, alpha_r)
                       + uniform_lpdf( mu | tau_min, tau_max);
  }
}
