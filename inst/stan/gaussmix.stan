functions {
  row_vector calc_logh(int N, int G, int K, vector pi, vector s, vector mu, vector tau, real dtau, matrix Mt) {
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
    return logh;
  }
}

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

transformed parameters {
  row_vector[N] logh;
  logh = calc_logh(N, G, K, pi, s, mu, tau, dtau, Mt);
}

model {
  pi ~ dirichlet(rep_vector(alpha_d,K));
  s ~ gamma(alpha_s,alpha_r);
  mu ~ uniform(tau_min,tau_max);
  target += sum(logh);
}
