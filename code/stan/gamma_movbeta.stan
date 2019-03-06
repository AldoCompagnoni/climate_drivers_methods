
data {
  int n_time;
  int<lower=0> n_lag;       // number of monthly lags
  vector[n_time] y;
  matrix[n_time,n_lag] clim;   // matrix of climate covariates
}

transformed data {
  real month[n_lag]; // vector of months, to create distance matrix
  for(k in 1:n_lag) month[k] = k;
}

parameters {
  real alpha;
  real mu_beta;
  real<lower=0> y_sd;
  real<lower=0> eta;  // maximum covariance for betas
  real<lower=0> rho;  // degree of temporal autocorrelation for betas
  vector[n_lag] z;    // unit normal prior for non-centered term
}

transformed parameters {
  
  // params for random beta
  vector[n_time] yhat;
  vector[n_lag] beta;
  matrix[n_lag,n_lag] sigma_beta; // covariance matrix
  matrix[n_lag,n_lag] L;     // cholesky of covariance matrix
  
  // covariance
  sigma_beta = cov_exp_quad(month, eta, rho) + diag_matrix(rep_vector(0.001, n_lag));
  L = cholesky_decompose(sigma_beta);
  
  // non-centered parameterization for beta
  beta = mu_beta + L * z;
  
  // linear predictor
  yhat = alpha + clim * beta;
  
}

model {
  // place holder  
  vector[n_time] mu; // transf. lin. pred. for mean
  
  // hyper-parameters to weight climate effects
  z ~ normal(0, 1);
  rho ~ normal(0, 5);
  eta ~ normal(0, 1);
  
  // parameters of data model
  alpha ~ normal(0, 5);
  mu_beta ~ normal(0, 5);
  y_sd ~ gamma(1,1); 

  for(n in 1:n_time)
    mu[n] = exp(yhat[n]);
    
  y ~ gamma(y_sd, y_sd ./ mu);
}

generated quantities {
  vector[n_time] log_lik;
  
  for (n in 1:n_time)
    log_lik[n] = gamma_lpdf(y[n] | y_sd, (y_sd / exp(yhat[n])) );
}
