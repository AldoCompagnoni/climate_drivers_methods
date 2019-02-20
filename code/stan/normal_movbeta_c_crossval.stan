
data {
  int n_train;
  int n_test;
  int n_lag;
  vector[n_train] y_train;
  vector[n_test]  y_test;
  matrix[n_train, n_lag] clim_train;
  matrix[n_test,  n_lag] clim_test;
}

parameters {
  real alpha;
  real<lower=0> y_sd;
  real mu_beta;
  real<lower=0> sigma_beta;
  vector[n_lag] beta;
}

transformed parameters {
  vector[n_train] yhat;
  
  // linear predictor
  yhat = alpha + clim_train * beta;
}

model {
  alpha ~ normal(0, 5);
  y_sd ~ gamma(1, 1);
  
  mu_beta ~ normal(0, 3);
  sigma_beta ~ normal(0, 3);
  
  beta ~ normal(mu_beta, sigma_beta);
  
  y_train ~ normal(yhat, y_sd);
}

generated quantities {
  vector[n_train] log_lik;
  vector[n_test]  pred_y;
  vector[n_test]  log_lik_test;
  
  
  // in sample prediction
  for (n in 1:n_train)
    log_lik[n] = normal_lpdf(y_train[n] | alpha + row(clim_train,n) * beta, y_sd);
  
  // out of sample prediction
  for(n in 1:n_test){
    pred_y[n]       = alpha + row(clim_test,n) * beta;
    log_lik_test[n] = normal_lpdf(y_test[n] | pred_y[n], y_sd);
  }
}
