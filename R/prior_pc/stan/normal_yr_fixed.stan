
data {
  int n_time;
  vector[n_time] clim_means;
  real gamma_shape;
}

generated quantities {
  
  real alpha = normal_rng(0,1);
  real beta  = normal_rng(0,1);
  real y_sd  = gamma_rng(gamma_shape,0.01); 
  
  real yhat[n_time];
  real y_sim[n_time];
  
  for(n in 1:n_time)
    yhat[n] = alpha + beta * clim_means[n];
  
  for(n in 1:n_time)
    y_sim[n] = normal_rng(yhat[n], y_sd);
  
}

