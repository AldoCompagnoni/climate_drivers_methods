
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
  vector[n_time] yhat;
  
  for(i in 1:n_time)
    x[i] = sum(theta .* clim[,i]);
    
  yhat = alpha + beta * x;
    
}

model {
  
  // priors  
  alpha ~ normal(0,0.5);
  beta  ~ normal(0,1);
  y_sd  ~ gamma(0.01,0.01); 
  theta ~ dirichlet(rep_vector(1.0, n_lag));
   
  // model
  y ~ normal(yhat, y_sd); 
   
}

generated quantities {
  vector[n_time] log_lik;
  
  for (n in 1:n_time)
  log_lik[n] = normal_lpdf(y[n] | yhat[n], y_sd );

}
