
data {
  int n_time;
  int n_lag;
  real alpha_dir;
  real gamma_shape;
  matrix[n_lag,n_time] clim;
}

transformed parameters {
  
    
}

generated quantities {
  
  real alpha           = normal_rng(0,1);
  real beta            = normal_rng(0,1);
  real<lower=0> y_sd   = gamma_rng(gamma_shape,0.01);
  simplex[n_lag] theta = dirichlet_rng(rep_vector(alpha_dir, n_lag));
  
  vector[n_time] x;
  vector[n_time] yhat;
  real y_sim[n_time];
  
  for(i in 1:n_time)
    x[i] = sum(theta .* clim[,i]);
    
  yhat = alpha + beta * x;
  
  y_sim = normal_rng(yhat, y_sd);
  
}
