
data {
  int n_time;
  vector[n_time] y;
  vector[n_time] clim_means;
}

parameters {
  real alpha;
  real beta;
  real<lower=0> sigma2;
}

transformed parameters {
  // place holder  
  vector[n_time] yhat; // transf. lin. pred. for mean
  real<lower=0> A[n_time];    // parameter for beta distn
  real<lower=0> B[n_time];    // parameter for beta distn
  
  yhat = exp(alpha + clim_means * beta);
  
  // beta reparameterization
  for(n in 1:n_time){
    // Parameters A and B 
    
    // BRMS parameterization
    A[n]      = yhat[n] * sigma2;
    B[n]      = sigma2 ;
    
  }
  
}


model {
  // parameters of data model
  alpha   ~ normal(0, 2);
  beta    ~ normal(0, 2);
  sigma2  ~ gamma(0.01,0.01); 

  y ~ gamma(A, B);
}

generated quantities {
  vector[n_time] log_lik;

  for (n in 1:n_time)
    log_lik[n] = gamma_lpdf( y[n] | A[n], B[n] );
}
