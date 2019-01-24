
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

transformed parameters{
  
  real<lower=0,upper=1> mu;    // transformed linear predictor for mean of beta distribution
  real<lower=0> A;               // parameter for beta distn
  real<lower=0> B;               // parameter for beta distn
  
  mu  = inv_logit(alpha);
  A   = mu * y_sd;
  B   = (1.0 - mu) * y_sd;
  
}

model {
  // model
  y_train ~ beta(A, B);
}

generated quantities {
  vector[n_train]  log_lik;
  vector[n_test]   pred_y;
  vector[n_test]   log_lik_test;

  // in sample log pointwise log likelihood
  for(n in 1:n_train)
    log_lik[n] = beta_lpdf(y_train[n] | A, B);
    
  // out of sample prediction
  for(n in 1:n_test){
    // out of sample predictions 
    pred_y[n]       = inv_logit(alpha);
    // out of sample ELPD
    log_lik_test[n] = beta_lpdf(y_test[n] | A, B);
  }

}
