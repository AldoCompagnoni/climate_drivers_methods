
data {
  int n_train;
  int n_test;
  int K;
  vector[n_train] y_train;
  vector[n_test] y_test;
  matrix[K,n_train] clim_yr_train;
  matrix[K,n_test] clim_yr_test;
}

parameters {
  real alpha;
  vector[K] beta;
  real<lower=0> y_sd;
}

model {
  real m[n_train];
  
  for(n in 1:n_train){
    m[n]  = alpha + sum(beta .* clim_yr_train[,n]);
  }
  
  // model
  y_train ~ normal(m, y_sd);
}

generated quantities {
  vector[n_train] log_lik;
  vector[n_test]  pred_y;
  vector[n_test]  log_lik_test;
  
  for (n in 1:n_train)
    log_lik[n] = normal_lpdf(y_train[n] | alpha + sum(beta .* clim_yr_train[,n]), y_sd);

  // out of sample prediction
  for(n in 1:n_test){
    pred_y[n]       = alpha + sum(beta .* clim_yr_test[,n]);
    log_lik_test[n] = normal_lpdf(y_test[n] | alpha + sum(beta .* clim_yr_test[,n]), y_sd);
  }
}
