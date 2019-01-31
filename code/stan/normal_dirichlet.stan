
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
  
  for(i in 1:n_time)
    x[i] = sum(theta .* clim[,i]);
}

model {
  
  // priors  
  alpha ~ normal(0,1);
  beta  ~ normal(0,1);
  y_sd  ~ gamma(1,1);
   
  // model
  y ~ normal(alpha + beta * x, y_sd);
}

// generated quantities {
//   vector[n_time] log_lik;
// 
//   for (n in 1:n_time)
//     log_lik[n] = normal_lpdf(y[n] | alpha + beta * x[n], y_sd);
// }
