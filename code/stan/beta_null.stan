
data {
  int n_time;
  vector[n_time] y;
}

parameters {
  real alpha;
  real<lower=0.1> y_sd;       // Flower-to-fruit overdispersion parameter
}

transformed parameters{
  
  real<lower=0,upper=1> mu[n_time];    // transformed linear predictor for mean of beta distribution
  real<lower=0> A[n_time];               // parameter for beta distn
  real<lower=0> B[n_time];               // parameter for beta distn
  
  for(n in 1:n_time){
    mu[n]  = inv_logit( alpha );
    A[n]   = mu[n] * y_sd;
    B[n]   = (1.0 - mu[n]) * y_sd;
  }
  
}

model {
  // model
  y ~ beta(A, B);
}

generated quantities {
  vector[n_time] log_lik;
  
  for (n in 1:n_time)
    log_lik[n] = beta_lpdf(y[n] | A, B);

}
