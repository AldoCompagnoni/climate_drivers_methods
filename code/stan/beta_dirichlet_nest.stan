
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
  
  // transformed parameters for beta binomial regression
  real<lower=0,upper=1> mu[n_time]; // transformed linear predictor for mean of beta distribution
  real<lower=0> A[n_time];          // parameter for beta distn
  real<lower=0> B[n_time];          // parameter for beta distn

  matrix[K,n_time] x_m;
  vector[n_time] x;
  
  for(i in 1:n_time){
    x_m[1,i] = sum(theta_m .* clim1[,i]); 
    x_m[2,i] = sum(theta_m .* clim2[,i]);
    x_m[3,i] = sum(theta_m .* clim2[,i]);
  }
  
  for(i in 1:n_time)
    x[i] = sum(theta_y .* x_m[,i]);
 
  for(n in 1:n_time){
    mu[n]  = inv_logit(alpha + x[n] * beta);
    A[n]   = mu[n] * y_sd;
    B[n]   = (1.0 - mu[n]) * y_sd;
  }
}

model {
  // likelihood
  y ~ beta(A, B);
}

generated quantities {
  vector[n_time] log_lik;
  
  for (n in 1:n_time)
    log_lik[n] = beta_lpdf(y[n] | A[n], B[n]);
}

