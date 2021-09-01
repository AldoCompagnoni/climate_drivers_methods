
data {
  int n_train;
  int n_test;
  int n_lag;
  vector[n_train] y_train;
  vector[n_test]  y_test;
  // real y_test;
  matrix[n_lag, n_train] clim_train;
  matrix[n_lag, n_test] clim_test;
  // row_vector[n_lag] clim_test;
}

parameters {
  simplex[n_lag] theta;
  real alpha;
  real beta;
  real<lower=0> y_sd;
}

transformed parameters {
  vector[n_train] x;
  
  for(i in 1:n_train)
    x[i] = sum(theta .* clim_train[,i]);
}

model {
  
  // priors  
  alpha ~ normal(0,0.5);
  beta  ~ normal(0,1);
  y_sd  ~ gamma(0.01,0.01); 
  theta ~ dirichlet(rep_vector(1.0, n_lag));
  
  // model
  y_train ~ normal(alpha + beta * x, y_sd);
}

generated quantities {
  vector[n_train] log_lik;
  vector[n_test]  pred_y;
  vector[n_test]  log_lik_test;
  real            pred_x;
  
  for(n in 1:n_train)
    log_lik[n] = normal_lpdf(y_train[n] | alpha + beta * x[n], y_sd);
  
  // out of sample prediction
  for(n in 1:n_test){
    pred_x          = sum(theta .* clim_test[,n]);
    pred_y[n]       = alpha + beta * pred_x;
    log_lik_test[n] = normal_lpdf(y_test[n] | alpha + beta * pred_x, y_sd);
  }
  
}
