data {
  int N; // Number of observations
  int G; // Number of grid pints
  int K; // Number of mixtures
  matrix[G,N] Mt; // Measurement matrix (transposed)
  vector[G] tau; // Calendar grid points for the measurement matrix calculation
  real<lower=0> alpha_d; // Concentration parameter for pi prior
  real taumin; // Lower calendar date. Same as min(tau)
  real taumax; // Upper calendar date. Same as max(tau)
  real<lower=0> alpha_s; // shape parameter of gamma distribution for sigma
  real<lower=0> alpha_r; // rate parameter of gamma distribution for sigma
}

parameters {
  simplex[K] pi;
  ordered[K] mu;
  vector<lower=0>[K] sig;
}

model {
  row_vector[N] L;
  row_vector[G] f;
  row_vector[G] logf;
  vector[K] logpi = log(pi);

  logf = to_row_vector(rep_vector(0,G));
  for (g in 1:G) {
    vector[K] lse = logpi;
    for (k in 1:K) {
      lse[k] += normal_lpdf(tau[g]|mu[k],sig[k]);
    }
    logf[g] = log_sum_exp(lse);
  }
  f = exp(logf);

  // Since the target only needs to be proportional to the likelihood,
  // there is no need to multiply by dtau on the following line.
  f = f / sum(f);

  pi ~ dirichlet(rep_vector(alpha_d,K));
  sig ~ gamma(alpha_s,alpha_r);
  mu ~ uniform(taumin,taumax);
  L = f * Mt;
  target += sum(log(L));
}
