
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
  // place holder  
  vector[n_time] mu;    // transformed linear predictor for mean of beta distribution
  
  // likelihood
  for(n in 1:n_time)
    mu[n] = exp(alpha + x[n] * beta);
    
  y ~ gamma(y_sd, y_sd ./ mu);
}

generated quantities {
  vector[n_time] log_lik;
  
  for (n in 1:n_time)
    log_lik[n] = gamma_lpdf(y[n] | y_sd, (y_sd / exp(alpha + x[n] * beta)) );
}
