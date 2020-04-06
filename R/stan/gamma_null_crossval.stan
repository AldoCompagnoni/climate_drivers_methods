
data {
  int n_test;
  int n_train;
  vector[n_train] y_train;
  vector[n_test] y_test;
}

parameters {
  real alpha;
  real<lower=0.1> y_sd;       // Flower-to-fruit overdispersion parameter
}

transformed parameters {
  // place holder (here for consistency with other STAN models) 
  real yhat; // transf. lin. pred. for mean
  
  yhat = exp(alpha);
  
}

model {
  
  // parameters of data model
  alpha ~ normal(0,1);
  y_sd  ~ gamma(1,1);
  
  y_train ~ gamma(y_sd, (y_sd ./ yhat) );
  
}

generated quantities {
  
  vector[n_train]  log_lik;
  vector[n_test]   pred_y;
  vector[n_test]   log_lik_test;

  for (n in 1:n_train)
    log_lik[n] = gamma_lpdf(y_train[n] | y_sd, (y_sd ./ yhat ) );
    
  // out of sample prediction
  for(n in 1:n_test){
    // out of sample predictions 
    pred_y[n]       = exp( alpha );
    // out of sample ELPD
    log_lik_test[n] = gamma_lpdf(y_test[n] | y_sd, (y_sd ./ exp( alpha ) ) );
  }

}
