
data {
  int<lower=1> n_time;
  vector<lower=0>[n_time] y;
}

parameters {
  real alpha;
  real<lower=0.1> y_sd;       // Flower-to-fruit overdispersion parameter
}

model {
  
  // place holder  
  vector[n_time] mu;    // transformed linear predictor for mean of beta distribution

  // model
  for(n in 1:n_time){
    mu[n] = exp( alpha );
  }
  y ~ gamma(y_sd, (y_sd ./ mu) );
  
}

generated quantities {
  vector[n_time] log_lik;
  
  for (n in 1:n_time)
    log_lik[n] = gamma_lpdf(y[n] | y_sd, (y_sd ./ exp(alpha)) );
    
}
