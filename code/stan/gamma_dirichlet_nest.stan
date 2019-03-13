
data {
  int n_time;
  int K;
  int M;
  vector[n_time] y;
  matrix[M,n_time] clim1;
  matrix[M,n_time] clim2;
  matrix[M,n_time] clim3;
}

parameters {
  simplex[K] theta_y;
  simplex[M] theta_m;
  real alpha;
  real beta;
  real<lower=0> y_sd;
}

transformed parameters {
  vector[n_time] x;
  vector[n_time] yhat;
  matrix[K,n_time] x_m;
  
  for(i in 1:n_time){
    x_m[1,i] = sum(theta_m .* clim1[,i]); 
    x_m[2,i] = sum(theta_m .* clim2[,i]);
    x_m[3,i] = sum(theta_m .* clim3[,i]);
  }
  
  for(i in 1:n_time)
    x[i] = sum(theta_y .* x_m[,i]);
    
  yhat = exp(alpha + x * beta);
  
}

model {
  // hyper-parameters to weight climate effects
  theta_y ~ dirichlet(rep_vector(1.0, K));
  theta_m ~ dirichlet(rep_vector(1.0, M));

  // priors
  alpha ~ normal(0,1);
  beta  ~ normal(0,1);
  y_sd  ~ gamma(1,1); 
  
  y ~ gamma(y_sd, y_sd ./ yhat);
  
}

generated quantities {
  vector[n_time] log_lik;
  
  for (n in 1:n_time)
    log_lik[n] = gamma_lpdf(y[n] | y_sd, (y_sd / yhat[n]) );
}
