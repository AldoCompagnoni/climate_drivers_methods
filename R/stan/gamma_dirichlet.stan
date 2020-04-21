
data {
  int n_time;
  int n_lag;
  vector[n_time] y;
  matrix[n_lag,n_time] clim;
}

parameters {
  simplex[n_lag] theta;
  real alpha;
  real beta;
  real<lower=0> y_sd;
}

transformed parameters {
  vector[n_time] x;
  vector[n_time] yhat;    // transformed linear predictor for mean of beta distribution
  
  for(i in 1:n_time)
    x[i] = sum(theta .* clim[,i]);
    
  yhat = exp(alpha + x * beta);
    
}

model {
  
  // priors
  alpha ~ normal(0,1);
  beta  ~ normal(0,1);
  y_sd  ~ gamma(1,1);
  theta ~ dirichlet(rep_vector(1.0, n_lag));
  
  // likelihood
  y ~ gamma(y_sd, y_sd ./ yhat);
}

generated quantities {
  vector[n_time] log_lik;
  
  for (n in 1:n_time)
    log_lik[n] = gamma_lpdf(y[n] | y_sd, (y_sd / exp(alpha + x[n] * beta)) );
}
