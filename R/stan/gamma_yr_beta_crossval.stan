
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
  // place holder  
  vector[n_train] mu;    // transformed linear predictor for mean of beta distribution
  
  // likelihood
  for(n in 1:n_train)
    mu[n] = exp(alpha + sum(beta .* clim_yr_train[,n]) );
    
  y_train ~ gamma(y_sd, y_sd ./ mu);
}

generated quantities {
  vector[n_train] log_lik;
  vector[n_test]  log_lik_test;
  vector[n_test]  pred_y;    // transformed linear predictor for mean of beta distribution

  for(n in 1:n_train)
    log_lik[n] = gamma_lpdf(y_train[n] | y_sd, (y_sd / exp(alpha + sum(beta .* clim_yr_train[,n]))) );

  for(n in 1:n_test){
    pred_y[n]       = exp(alpha + sum(beta .* clim_yr_test[,n]));
    log_lik_test[n] = gamma_lpdf(y_test[n] | y_sd, (y_sd / exp(alpha + sum(beta .* clim_yr_test[,n]))) );
  }
}
