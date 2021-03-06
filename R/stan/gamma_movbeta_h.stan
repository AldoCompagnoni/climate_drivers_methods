
data {
  int<lower=0> n_time;
  int<lower=0> n_lag;       // number of monthly lags 
  vector[n_time] y;
  matrix[n_time,n_lag] clim;  // matrix of climate covariates
}

parameters {
  real alpha;
  real<lower=0> y_sd;
  real mu_beta;
  real<lower=0> sigma_beta;
  vector[n_lag] z;    // unit normal prior for non-centered term
}

transformed parameters {
  
  // params for random beta
  vector[n_time] yhat;
  vector[n_lag] beta;
  
  // non-centered parameterization
  beta = mu_beta + sigma_beta * z;
  
  // linear predictor
  yhat = exp(alpha + clim * beta);
  
}

model {
  
  // hyper-parameters to weight climate effects
  z ~ normal(0, 1);
  mu_beta ~ normal(0, 3);
  sigma_beta ~ normal(0, 3);
  
  // priors
  alpha ~ normal(0,1);
  y_sd  ~ gamma(1,1); 

  y ~ gamma(y_sd, y_sd ./ yhat);
}

generated quantities {
  vector[n_time] log_lik;
  
  for (n in 1:n_time)
    log_lik[n] = gamma_lpdf(y[n] | y_sd, (y_sd / yhat[n]) );
}
