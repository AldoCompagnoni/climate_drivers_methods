
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
  vector[n_time] yhat;
  
  // linear predictor
  yhat = alpha + clim * beta;
}

model {
  alpha ~ normal(0, 5);
  y_sd ~ gamma(1, 1);
  
  mu_beta ~ normal(0, 3);
  sigma_beta ~ normal(0, 3);
  
  beta ~ normal(mu_beta, sigma_beta);
  
  y ~ normal(yhat, y_sd);
}

generated quantities {
  vector[n_time] log_lik;
  
  for (n in 1:n_time) {
    log_lik[n] = normal_lpdf(y[n] | yhat[n], y_sd);
  }
}
