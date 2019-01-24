
data {
  int n_time;
  int K;
  vector[n_time] y;
  matrix[K,n_time] clim_yr;
}

parameters {
  simplex[K] theta;
  real alpha;
  real beta;
  real<lower=0> y_sd;
}

transformed parameters {
  vector[n_time] x;
  
  for(i in 1:n_time)
    x[i] = sum(theta .* clim_yr[,i]);
}

model {
  // model
  y ~ normal(alpha + beta * x, y_sd);
}

generated quantities {
  vector[n_time] log_lik;

  for (n in 1:n_time)
    log_lik[n] = normal_lpdf(y[n] | alpha + beta * x[n], y_sd);
}
