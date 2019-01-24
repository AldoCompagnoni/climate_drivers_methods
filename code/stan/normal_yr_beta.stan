
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
  real m[n_time];
  
  for(n in 1:n_time){
    m[n]  = alpha + sum(beta .* clim_yr[,n]);
  }
  
  // model
  y ~ normal(m, y_sd);
}

generated quantities {
  vector[n_time] log_lik;

  for (n in 1:n_time)
    log_lik[n] = normal_lpdf(y[n] | alpha + sum(beta .* clim_yr[,n]), y_sd);
}
