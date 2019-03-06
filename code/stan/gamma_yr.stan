
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

model {
  // place holder  
  vector[n_time] mu; // transf. lin. pred. for mean
  
  // parameters of data model
  alpha ~ normal(0,1);
  beta  ~ normal(0,1);
  y_sd  ~ gamma(1,1); 

  for(n in 1:n_time)
    mu[n] = exp(alpha + clim_means[n] * beta);
    
  y ~ gamma(y_sd, y_sd ./ mu);
}

generated quantities {
  vector[n_time] log_lik;
  
  for (n in 1:n_time)
    log_lik[n] = gamma_lpdf(y[n] | y_sd, (y_sd / exp(alpha + clim_means[n] * beta)) );
}
