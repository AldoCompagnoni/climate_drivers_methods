
data {
  
  int<lower=1> n;
  vector<lower=0>[n] y;
  vector[n] x1;
  
}

parameters {

  real beta[2];
  real<lower=0.1> y_sd;       // Flower-to-fruit overdispersion parameter
  
}

model {
  
  // place holder  
  vector[n] mu;    // transformed linear predictor for mean of beta distribution

  // likelihood
  for(ny in 1:n){
    mu[ny] = exp( beta[1] + beta[2] * x1[ny]);
  }
  y ~ gamma(y_sd, (y_sd ./ mu) );
  
}

generated quantities {
  vector[n] log_lik;
  
  for (ny in 1:n)
    log_lik[ny] = gamma_lpdf(y[ny] | y_sd, (y_sd / exp(beta[1] + beta[2] * x1[ny])) );

}
