
data {
  int n_time;
  int n_lag;
  real alpha_dir;
  real gamma_shape;
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
  alpha ~ normal(0,1);
  beta  ~ normal(0,1);
  y_sd  ~ gamma(gamma_shape,0.01); 
  theta ~ dirichlet(rep_vector(alpha_dir, n_lag));
   
}

generated quantities {
  real y_sim[n_time];

  y_sim = normal_rng(yhat, y_sd);
  
}
