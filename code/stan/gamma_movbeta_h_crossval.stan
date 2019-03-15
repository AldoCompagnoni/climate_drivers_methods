
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
  vector[n_lag] z;    // unit normal prior for non-centered term
}

transformed parameters {
  
  // params for random beta
  vector[n_train] yhat;
  vector[n_lag] beta;
  
  // non-centered parameterization
  beta = mu_beta + sigma_beta * z;
  
  // linear predictor
  yhat = exp(alpha + clim_train * beta);
  
}

model {
  
  // hyper-parameters to weight climate effects
  z ~ normal(0, 1);
  mu_beta ~ normal(0, 3);
  sigma_beta ~ normal(0, 3);
  
  // priors
  alpha ~ normal(0,1);
  y_sd  ~ gamma(1,1); 

  y_train ~ gamma(y_sd, y_sd ./ yhat);
}

generated quantities {
  vector[n_train] log_lik;
  vector[n_test]  log_lik_test;
  vector[n_test]  pred_y;    // transf. lin. pred. for mean of gamma distrib.

  for (n in 1:n_train)
    log_lik[n] = gamma_lpdf(y_train[n] | y_sd, (y_sd / yhat[n]) );

  for(n in 1:n_test){
    pred_y[n]       = exp(alpha + row(clim_test,n) * beta);
    log_lik_test[n] = gamma_lpdf(y_test[n] | y_sd, (y_sd / pred_y[n]) );
  }
  
}
