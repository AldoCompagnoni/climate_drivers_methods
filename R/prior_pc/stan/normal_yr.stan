
data {
  int n_time;
  vector[n_time] clim_means;
  real gamma_shape;
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
  alpha ~ normal(0,1);
  beta  ~ normal(0,1);
  y_sd  ~ gamma(gamma_shape,0.01); 
  
}

generated quantities {
  real y_sim[n_time];

  y_sim = normal_rng(yhat, y_sd);
}

