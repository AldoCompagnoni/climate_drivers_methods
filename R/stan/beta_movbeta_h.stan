
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
  vector[n_lag] z;    // unit normal prior for non-centered term
}

transformed parameters {
  
  // transformed parameters for beta binomial regression
  real<lower=0,upper=1> yhat[n_time]; // transf. lin. pred. for mean of beta distribution
  real<lower=0> A[n_time];          // parameter for beta distn
  real<lower=0> B[n_time];          // parameter for beta distn

  // params for random beta
  vector[n_time] mu;
  vector[n_lag] beta;
  
  // non-centered parameterization
  beta = mu_beta + sigma_beta * z;
  
  // linear predictor
  mu = alpha + clim * beta;
  
  // beta reparameterization
  for(n in 1:n_time){
    yhat[n] = inv_logit(mu[n]);
    A[n]    = yhat[n] * y_sd;
    B[n]    = (1.0 - yhat[n]) * y_sd;
  }
  
}

model {
  
  // hyper-parameters to weight climate effects
  z ~ normal(0, 1);
  mu_beta ~ normal(0, 3);
  sigma_beta ~ normal(0, 3);
  
  // parameters of data model
  alpha ~ normal(0, 5);
  y_sd ~ gamma(1, 1);
  
  y ~ beta(A, B);
}

generated quantities {
  vector[n_time] log_lik;
  
  for (n in 1:n_time) {
    log_lik[n] = beta_lpdf(y[n] | A[n], B[n]);
  }
}
