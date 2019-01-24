
data {
  //int n_train;
  int n_test;
  int n_train;
  vector[n_train] y_train;
  vector[n_test] y_test;
}

parameters {
  real alpha;
  real<lower=0.1> y_sd;       // Flower-to-fruit overdispersion parameter
}

model {
  
  // place holder  
  vector[n_train] mu;    // transformed linear predictor for mean of beta distribution

  // model
  for(n in 1:n_train){
    mu[n] = exp( alpha );
  }
  y_train ~ gamma(y_sd, (y_sd ./ mu) );
  
}

generated quantities {
  
  vector[n_train]  log_lik;
  vector[n_test]   pred_y;
  vector[n_test]   log_lik_test;

  for (n in 1:n_train)
    log_lik[n] = gamma_lpdf(y_train[n] | y_sd, (y_sd ./ exp(alpha)) );
    
  // out of sample prediction
  for(n in 1:n_test){
    // out of sample predictions 
    pred_y[n]       = exp( alpha );
    // out of sample ELPD
    log_lik_test[n] = gamma_lpdf(y_test[n] | y_sd, (y_sd ./ exp(alpha)) );
  }

}
