
data {
  int n_train;
  int n_test;
  int n_lag;
  vector[n_train] y_train;
  vector[n_test]  y_test;
  matrix[n_train, n_lag] clim_train;
  matrix[n_test,  n_lag] clim_test;
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
  vector[n_train] yhat;
  vector[n_lag] beta;
  matrix[n_lag,n_lag] Sigma; // covariance matrix
  matrix[n_lag,n_lag] L;     // cholesky of covariance matrix
  
  // covariance
  Sigma = cov_exp_quad(month, eta, rho) + diag_matrix(rep_vector(0.001, n_lag));
  L = cholesky_decompose(Sigma);
  
  // non-centered parameterization for beta
  beta = mu_beta + L * z;
  
  // linear predictor
  yhat = alpha + clim_train * beta;
}

model {
  
  // priors
  z ~ normal(0, 1);
  
  rho ~ normal(0, 5);
  eta ~ normal(0, 1);
  
  alpha ~ normal(0, 5);
  mu_beta ~ normal(0, 5);
  y_sd ~ gamma(1, 1);

  // model
  y_train ~ normal(yhat, y_sd);
  
}

generated quantities {
  vector[n_train] log_lik;
  vector[n_test]  log_lik_test;
  
  // in sample prediction
  for (n in 1:n_train)
    log_lik[n] = normal_lpdf(y_train[n] | alpha + row(clim_train,n) * beta, y_sd);
  
  // out of sample prediction
  for(n in 1:n_test){
    log_lik_test[n] = normal_lpdf(y_test[n] | alpha + row(clim_test,n) * beta, y_sd);
  }
  
}
