
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
  
  // transformed parameters for beta binomial regression
  real<lower=0,upper=1> mu[n_train]; // transf. lin. pred. for mean of beta distribution
  real<lower=0> A[n_train];          // parameter for beta distn
  real<lower=0> B[n_train];          // parameter for beta distn

  // params for random beta
  vector[n_train] yhat;
  
  // linear predictor
  yhat = alpha + clim_train * beta;
  
  // beta reparameterization
  for(n in 1:n_train){
    mu[n]  = inv_logit(yhat[n]);
    A[n]   = mu[n] * y_sd;
    B[n]   = (1.0 - mu[n]) * y_sd;
  }
  
}

model {
  alpha ~ normal(0, 5);
  y_sd ~ gamma(1, 1);
  
  mu_beta ~ normal(0, 3);
  sigma_beta ~ normal(0, 3);
  
  beta ~ normal(mu_beta, sigma_beta);
  
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
    pred_y[n]       = inv_logit(alpha + row(clim_test,n) * beta);
    A_test[n]       = pred_y[n] * y_sd;
    B_test[n]       = (1.0 - pred_y[n]) * y_sd;
    log_lik_test[n] = beta_lpdf(y_test[n] | A_test[n], B_test[n]);
  }
  
}
