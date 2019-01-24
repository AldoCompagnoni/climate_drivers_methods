
data {
  int n_train;
  int n_test;
  vector[n_train] y_train;
  vector[n_test] y_test;
  vector[n_train] clim_means_train;
  vector[n_test] clim_means_test;
}

parameters {
  real alpha;
  real beta;
  real<lower=0> y_sd;
}

model {
  // model
  y_train ~ normal(alpha + beta * clim_means_train, y_sd);
}

generated quantities {
  vector[n_train] log_lik;
  vector[n_test]  pred_y;
  vector[n_test]  log_lik_test;
  
  for (n in 1:n_train)
    log_lik[n] = normal_lpdf(y_train[n] | alpha + beta * clim_means_train[n], y_sd);
    
  // out of sample prediction
  for(n in 1:n_test){
    pred_y[n]       = alpha + beta * clim_means_test[n];
    log_lik_test[n] = normal_lpdf(y_test[n] | alpha + beta * clim_means_test[n], y_sd);
  }
}

