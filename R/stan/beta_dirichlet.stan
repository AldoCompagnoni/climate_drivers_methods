
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
  
  // transformed parameters for beta binomial regression
  real<lower=0,upper=1> yhat[n_time]; // transformed linear predictor for mean of beta distribution
  real<lower=0> A[n_time];          // parameter for beta distn
  real<lower=0> B[n_time];          // parameter for beta distn

  vector[n_time] x;
  
  for(i in 1:n_time)
    x[i] = sum(theta .* clim[,i]);
    
  for(n in 1:n_time){
    yhat[n]  = inv_logit(alpha + x[n] * beta);
    A[n]   = yhat[n] * y_sd;
    B[n]   = (1.0 - yhat[n]) * y_sd;
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
