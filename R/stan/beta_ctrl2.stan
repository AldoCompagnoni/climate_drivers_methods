
data {
  
  int<lower=1> n_time;
  vector<lower=0, upper=1>[n_time] y;
  vector[n_time] clim_means;
  
}

parameters {
  
  real alpha;
  real beta;
  real<lower=0.1> y_sd;       // Flower-to-fruit overdispersion parameter
  
}

transformed parameters{
  
  real<lower=0,upper=1> mu[n_time];    // transformed linear predictor for mean of beta distribution
  real<lower=0> A[n_time];               // parameter for beta distn
  real<lower=0> B[n_time];               // parameter for beta distn
  
  for(n in 1:n_time){
    mu[n]  = inv_logit(alpha + clim_means[n] * beta);
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
