
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
  simplex[K] theta;
  real alpha;
  real beta;
  real<lower=0> y_sd;
}

transformed parameters {
  vector[n_train] x;
  
  for(i in 1:n_train)
    x[i] = sum(theta .* clim_yr_train[,i]);
}

model {
  // place holder  
  vector[n_train] mu; // transf. lin. pred. for mean
  
  // likelihood
  for(n in 1:n_train)
    mu[n] = exp(alpha + x[n] * beta);
    
  y_train ~ gamma(y_sd, y_sd ./ mu);
}

generated quantities {
  vector[n_train] log_lik;
  vector[n_test]  log_lik_test;
  vector[n_test]  pred_y;    // transformed linear predictor for mean of beta distribution
  vector[n_test]  pred_x;

  for (n in 1:n_train)
    log_lik[n] = gamma_lpdf(y_train[n] | y_sd, (y_sd / exp(alpha + x[n] * beta)) );

  for(n in 1:n_test){
    pred_x[n]       = sum(theta .* clim_yr_test[,n]);
    pred_y[n]       = exp(alpha + pred_x[n] * beta);
    log_lik_test[n] = gamma_lpdf(y_test[n] | y_sd, (y_sd / exp(alpha + pred_x[n] * beta)) );
  }
}
