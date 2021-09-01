
data {
  int n_time;
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
  
  // priors
  alpha ~ normal(0,1);
  y_sd  ~ gamma(0.01,0.01); 
  
}

generated quantities {
  real y_sim[n_time];

  for(n in 1:n_time)
    y_sim[n] = normal_rng(alpha, y_sd);
  
}
