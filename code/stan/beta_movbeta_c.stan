
data {
  int<lower=0> n_time; // number of data points, length(y)
  int<lower=0> n_lag; // number of monthly lags
  vector[n_time] y;    // response
  matrix[n_time,n_lag] clim;  // matrix of climate covariates
}

parameters {
  real alpha;
  real<lower=0> y_sd;
  real mu_beta;
  real<lower=0> sigma_beta;
  vector[n_lag] beta;
}

transformed parameters {
  
  // transformed parameters for beta binomial regression
  real<lower=0,upper=1> mu[n_time]; // transf. lin. pred. for mean of beta distribution
  real<lower=0> A[n_time];          // parameter for beta distn
  real<lower=0> B[n_time];          // parameter for beta distn

  // params for random beta
  vector[n_time] yhat;
  
  // linear predictor
  yhat = alpha + clim * beta;
  
  // beta reparameterization
  for(n in 1:n_time){
    mu[n]  = inv_logit(yhat[n]);
    A[n]   = mu[n] * y_sd;
    B[n]   = (1.0 - mu[n]) * y_sd;
  }
  
}

model {
  alpha ~ normal(0, 5);
  y_sd ~ gamma(1, 1);
  
  mu_beta ~ normal(0, 3);
  sigma_beta ~ normal(0, 3);
  
  beta ~ normal(mu_beta, sigma_beta);
  
  y ~ beta(A, B);
}

generated quantities {
  vector[n_time] log_lik;
  
  for (n in 1:n_time) {
    log_lik[n] = beta_lpdf(y[n] | yhat[n], y_sd);
  }
}
