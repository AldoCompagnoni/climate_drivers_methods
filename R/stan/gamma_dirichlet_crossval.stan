
data {
  int n_train;
  int n_test;
  int n_lag;
  vector[n_train] y_train;
  vector[n_test]  y_test;
  // real y_test;
  matrix[n_lag, n_train] clim_train;
  matrix[n_lag, n_test] clim_test;
}

parameters {
  simplex[n_lag] theta;
  real alpha;
  real beta;
  real<lower=0> y_sd;
}

transformed parameters{
  vector[n_train] x;
  
  for(i in 1:n_train)
    x[i] = sum(theta .* clim_train[,i]);
}

model{
  
  // place holder  
  vector[n_train] mu;    // transformed linear predictor for mean of beta distribution
  
  // priors
  alpha ~ normal(0,0.5);
  beta  ~ normal(0,1);
  y_sd  ~ gamma(0.01,0.01); 
  theta ~ dirichlet(rep_vector(1.0, n_lag));
  
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
    pred_x[n]       = sum(theta .* clim_test[,n]);
    pred_y[n]       = exp(alpha + pred_x[n] * beta);
    log_lik_test[n] = gamma_lpdf(y_test[n] | y_sd, (y_sd / exp(alpha + pred_x[n] * beta)) );
  }
}
