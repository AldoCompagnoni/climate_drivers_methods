
data {
  int n_time;
  int K;
  vector[n_time] y;
  matrix[K,n_time] clim_yr;
}

parameters {
  real alpha;
  vector[K] beta;
  real<lower=0> y_sd;
}

transformed parameters {
  
  // transformed parameters for beta binomial regression
  real<lower=0,upper=1> mu[n_time]; // transformed linear predictor for mean of beta distribution
  real<lower=0> A[n_time];          // parameter for beta distn
  real<lower=0> B[n_time];          // parameter for beta distn

  for(n in 1:n_time){
    mu[n]  = inv_logit(alpha + sum(beta .* clim_yr[,n]));
    A[n]   = mu[n] * y_sd;
    B[n]   = (1.0 - mu[n]) * y_sd;
  }
}

model {
  // likelihood
  y ~ beta(A, B);
}

generated quantities {
  vector[n_time] log_lik;
  
  for (n in 1:n_time)
    log_lik[n] = beta_lpdf(y[n] | A[n], B[n]);
}
