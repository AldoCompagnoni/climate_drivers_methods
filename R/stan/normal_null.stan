
data {
  int n_time;
  vector[n_time] y;
  
}

parameters {
  real alpha;
  real<lower=0> y_sd;
}

transformed parameters {
  // just to facilitate pipeline for plots
  real yhat;
  yhat = alpha;
}


model {
  
  real mu[n_time];
  
  for(ny in 1:n_time){
    mu[ny] = alpha;
  }
  
  // priors
  alpha ~ normal(0,0.5);
  y_sd  ~ gamma(0.01,0.01); 
  
  // model
  y ~ normal(mu, y_sd);
}

generated quantities {
  vector[n_time] log_lik;
  
  for (n in 1:n_time)
  log_lik[n] = normal_lpdf(y[n] | alpha, y_sd );

}
