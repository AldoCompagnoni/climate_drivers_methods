
data {
  
  int<lower=1> n_time;
  vector<lower=0>[n_time] y;
  vector[n_time] clim_means;
  
}

parameters {
  
  real alpha;
  real beta;
  real<lower=0.1> y_sd;       // Flower-to-fruit overdispersion parameter
  
}

model {
  
  // place holder  
  vector[n_time] mu;    // transformed linear predictor for mean of beta distribution

  // likelihood
  for(n in 1:n_time)
    mu[n] = exp(alpha + clim_means[n] * beta);
    
  y ~ normal(mu, y_sd);
  
}

generated quantities {
  vector[n_time] log_lik;
  
  for (n in 1:n_time)
    log_lik[n] = normal_lpdf(y[n] | exp(alpha + clim_means[n] * beta), y_sd);

}
