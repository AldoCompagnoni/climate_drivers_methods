
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
  // model
  y ~ normal(alpha + beta * clim_means, y_sd);
}

generated quantities {
  vector[n_time] log_lik;
  
  for (n in 1:n_time)
    log_lik[n] = normal_lpdf(y[n] | alpha + beta * clim_means[n], y_sd);
}

