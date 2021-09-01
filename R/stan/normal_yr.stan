
data {
  int n_time;
  vector[n_time] y;
  vector[n_time] clim_means;
}

parameters {
  real alpha;
  real beta;
  real<lower=0> y_sd;
}

transformed parameters{
  vector[n_time] yhat;
  
  yhat = alpha + beta * clim_means;
}

model {
  
  // priors
  alpha ~ normal(0,0.5);
  beta  ~ normal(0,1);
  y_sd  ~ gamma(0.01,0.01); 
  
  // model
  y ~ normal(yhat, y_sd);
}

generated quantities {
  vector[n_time] log_lik;
  
  for (n in 1:n_time)
  log_lik[n] = normal_lpdf(y[n] | yhat[n], y_sd );

}
