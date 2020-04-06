
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

transformed parameters {
  
  // transformed parameters for beta binomial regression
  real<lower=0,upper=1> mu[n_train]; // transformed linear predictor for mean of beta distribution
  real<lower=0> A[n_train];          // parameter for beta distn
  real<lower=0> B[n_train];          // parameter for beta distn

  for(n in 1:n_train){
    mu[n]  = inv_logit(alpha + sum(beta .* clim_yr_train[,n]));
    A[n]   = mu[n] * y_sd;
    B[n]   = (1.0 - mu[n]) * y_sd;
  }
}

model {
  // likelihood
  y_train ~ beta(A, B);
}

generated quantities {
  vector[n_train] log_lik;
  vector[n_test]  log_lik_test;
  vector[n_test]  pred_y;
  real<lower=0>   A_test[n_test];               // parameter for beta distn
  real<lower=0>   B_test[n_test];               // parameter for beta distn
  
  for (n in 1:n_train)
    log_lik[n] = beta_lpdf(y_train[n] | A[n], B[n]);
  
  // out of sample prediction
  for(n in 1:n_test){
    pred_y[n]       = inv_logit(alpha + sum(beta .* clim_yr_test[,n]) );
    A_test[n]       = pred_y[n] * y_sd;
    B_test[n]       = (1.0 - pred_y[n]) * y_sd;
    log_lik_test[n] = beta_lpdf(y_test[n] | A_test[n], B_test[n]);
  }
}
