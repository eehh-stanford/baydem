data {
  int N; // Number of observations
  int G; // Number of grid pints
  int K; // Number of mixtures
  matrix[G,N] Mt; // Measurement matrix (transposed)
  vector[G] ygrid; // Grid points for the measurement matrix
  real<lower=0> dirichParam; // Parameter for the mixture distribution
  real ymin;
  real ymax;
  real<lower=0> sigAlpha;
  real<lower=0> sigBeta;
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
      lse[k] += normal_lpdf(ygrid[g]|mu[k],sig[k]);
    }
    logf[g] = log_sum_exp(lse);
  }
  f = exp(logf);

  
  // TODO: Use the trapezoid rule to determine the normalization rather than dividing by sum(f)
  f = f / sum(f);

  pi ~ dirichlet(rep_vector(dirichParam,K));
  sig ~ gamma(sigAlpha,sigBeta);
  mu ~ uniform(ymin,ymax);
  L = f * Mt;
  target += sum(log(L));
}
