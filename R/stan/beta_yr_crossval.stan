
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

transformed parameters{
  
  real<lower=0,upper=1> mu[n_train];    // transformed linear predictor for mean of beta distribution
  real<lower=0> A[n_train];               // parameter for beta distn
  real<lower=0> B[n_train];               // parameter for beta distn
  
  for(n in 1:n_train){
    mu[n]  = inv_logit(alpha + clim_means_train[n] * beta);
    A[n]   = mu[n] * y_sd;
    B[n]   = (1.0 - mu[n]) * y_sd;
  }
}

model {
  
  // priors
  alpha ~ normal(0,1);
  beta  ~ normal(0,1);
  y_sd  ~ gamma(1,1); 
  
  // model
  y_train ~ beta(A, B);
}

generated quantities {
  vector[n_train] log_lik;
  vector[n_test]  pred_y;
  vector[n_test]  log_lik_test;
  real<lower=0> A_test[n_test];               // parameter for beta distn
  real<lower=0> B_test[n_test];               // parameter for beta distn

  for(n in 1:n_train)
    log_lik[n] = beta_lpdf(y_train[n] | A[n], B[n]);
    
  for(n in 1:n_test){
    pred_y[n]  = inv_logit(alpha + clim_means_test[n] * beta);
    A_test[n]   = pred_y[n] * y_sd;
    B_test[n]   = (1.0 - pred_y[n]) * y_sd;
    log_lik_test[n] = beta_lpdf(y_test[n] | A_test[n], B_test[n]);
  }
}
