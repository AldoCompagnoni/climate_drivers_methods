
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
  real<lower=0.1> y_sd;       // Flower-to-fruit overdispersion parameter
  
}

model {
  
  // place holder  
  vector[n_train] mu;    // transformed linear predictor for mean of beta distribution

  // likelihood
  for(n in 1:n_train)
    mu[n]  = exp(alpha + clim_means_train[n] * beta);
    
  y_train ~ normal(mu, y_sd);
  
}

generated quantities {
  
  vector[n_train] log_lik;
  vector[n_test]  log_lik_test;
  vector[n_test]  pred_y;    // transformed linear predictor for mean of beta distribution

  for (n in 1:n_train)
    log_lik[n] = normal_lpdf(y_train[n] | exp(alpha + clim_means_train[n] * beta), y_sd);

  for(n in 1:n_test){
    pred_y[n]       = exp(alpha + clim_means_test[n] * beta);
    log_lik_test[n] = normal_lpdf(y_test[n] | exp(alpha + clim_means_test[n] * beta), y_sd);
  }

}
