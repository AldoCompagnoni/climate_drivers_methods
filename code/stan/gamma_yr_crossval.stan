
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

model{
  // place holder  
  vector[n_time] mu; // transf. lin. pred. for mean
  
  // parameters of data model
  alpha ~ normal(0,1);
  beta  ~ normal(0,1);
  y_sd  ~ gamma(1,1);
  
  // likelihood
  for(n in 1:n_train)
    mu[n] = exp(alpha + clim_means_train[n] * beta);
    
  y_train ~ gamma(y_sd, y_sd ./ mu);
}

generated quantities{
  vector[n_train] log_lik;
  vector[n_test]  pred_y;
  vector[n_test]  log_lik_test;

  for(n in 1:n_train)
    log_lik[n] = gamma_lpdf(y_train[n] | y_sd, (y_sd / exp(alpha + clim_means_train[n] * beta)) );
    
  for(n in 1:n_test){
    pred_y[n]  = exp(alpha + clim_means_test[n] * beta);
    log_lik_test[n] = gamma_lpdf(y_test[n] | y_sd, (y_sd / exp(alpha + clim_means_test[n] * beta)) );
  }
}
