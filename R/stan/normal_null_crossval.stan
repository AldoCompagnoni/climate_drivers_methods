  
data {
  //int n_time;
  int n_test;
  int n_train;
  vector[n_train] y_train;
  vector[n_test] y_test;
}

parameters {
  real alpha;
  real<lower=0> y_sd;
}

model {
  
  // priors
  alpha ~ normal(0,0.5);
  y_sd  ~ gamma(0.01,0.01); 
  
  // model
  y_train ~ normal(alpha, y_sd);
}

generated quantities {
  vector[n_train] log_lik;
  vector[n_test]  pred_y;
  vector[n_test]  log_lik_test;
  
  for(n in 1:n_train)
    log_lik[n] = normal_lpdf(y_train[n] | alpha, y_sd);
    
  // out of sample prediction
  for(n in 1:n_test){
    pred_y[n] = alpha;
    log_lik_test[n] = normal_lpdf(y_test[n] | alpha, y_sd);
  }

}
