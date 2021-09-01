
data {
  int n_time;
  vector[n_time] y;
}

parameters {
  real alpha;
  real<lower=0.1> y_sd;       // dispersion parameter
}

transformed parameters{
  
  real<lower=0,upper=1> yhat[n_time]; // transf. lin. pred. for mean of beta distrib.
  real<lower=0> A[n_time];          // parameter for beta distn
  real<lower=0> B[n_time];          // parameter for beta distn
  
  for(n in 1:n_time){
    yhat[n]  = inv_logit( alpha );
    A[n]     = yhat[n] * y_sd;
    B[n]     = (1.0 - yhat[n]) * y_sd;
  }
  
}

model {
  
  // priors
  alpha ~ normal(0,0.5);
  y_sd  ~ gamma(1,1); 
  
  // model
  y ~ beta(A, B);
}

generated quantities {
  vector[n_time] log_lik;
  
  for (n in 1:n_time)
    log_lik[n] = beta_lpdf(y[n] | A, B);

}
