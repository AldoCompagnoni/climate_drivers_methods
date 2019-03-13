
data {
  int n_time;
  int K;
  int M;
  vector[n_time] y;
  vector[n_time] x0;
  matrix[M,n_time] clim1;
  matrix[M,n_time] clim2;
  matrix[M,n_time] clim3;
}

parameters {
  simplex[K] theta_y;
  simplex[M] theta_m;
  real alpha;
  real beta;
  real dd;
  real<lower=0> y_sd;
}

transformed parameters {
  matrix[K,n_time] x_m;
  vector[n_time] x;
  
  for(i in 1:n_time){
    x_m[1,i] = sum(theta_m .* clim1[,i]); 
    x_m[2,i] = sum(theta_m .* clim2[,i]);
    x_m[3,i] = sum(theta_m .* clim3[,i]);
  }
  
  for(i in 1:n_time)
    x[i] = sum(theta_y .* x_m[,i]);
}

model {
  
  // priors
  alpha   ~ normal(0,1);
  beta    ~ normal(0,1);
  dd      ~ normal(1,1);
  y_sd    ~ gamma(1,1);
  theta_y ~ dirichlet(rep_vector(1.0, K));
  theta_m ~ dirichlet(rep_vector(1.0, M));

  // model
  y ~ normal(alpha + beta * x + dd * x0, y_sd);
}

// generated quantities {
//   vector[n_time] log_lik;
// 
//   for (n in 1:n_time)
//     log_lik[n] = normal_lpdf(y[n] | alpha + beta * x[n], y_sd);
// }
// 
