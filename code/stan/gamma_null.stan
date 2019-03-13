
data {
  int<lower=1> n_time;
  vector<lower=0>[n_time] y;
}

parameters {
  real alpha;
  real<lower=0.1> y_sd;
}

transformed parameters{
  // place holder  
  real yhat; // transf. lin. pred. for mean

  yhat = exp( alpha );

}


model {
  
  // priors
  alpha ~ normal(0,1);
  y_sd  ~ gamma(1,1); 

  y ~ gamma(y_sd, (y_sd ./ yhat) );
  
}

generated quantities {
  vector[n_time] log_lik;
  
  for (n in 1:n_time)
    log_lik[n] = gamma_lpdf(y[n] | y_sd, (y_sd ./ exp(alpha)) );
    
}
