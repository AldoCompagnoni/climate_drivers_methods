
data {
  int n_time;
  int K;
  vector[n_time] y;
  matrix[K,n_time] clim_yr;
}

parameters {
  real alpha;
  vector[K] beta;
  real<lower=0> y_sd;
}

model {
  // place holder  
  vector[n_time] mu;    // transformed linear predictor for mean of beta distribution
  
  // likelihood
  for(n in 1:n_time)
    mu[n] = exp(alpha + sum(beta .* clim_yr[,n]) );
    
  y ~ gamma(y_sd, y_sd ./ mu);
}

generated quantities {
  vector[n_time] log_lik;
  
  for (n in 1:n_time)
    log_lik[n] = gamma_lpdf(y[n] | y_sd, (y_sd / exp(alpha + sum(beta .* clim_yr[,n])) ) );
}
