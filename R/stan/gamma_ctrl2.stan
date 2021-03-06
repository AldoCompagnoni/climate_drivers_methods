
data {
  
  int<lower=1> n_time;
  vector<lower=0>[n_time] y;
  vector[n_time] clim_means;
  
}

parameters {
  
  real alpha;
  real beta;
  real<lower=0.1> y_sd;
  
}

model {
  
  // place holder  
  vector[n_time] mu;    // transformed linear predictor for mean of beta distribution

  // likelihood
  for(n in 1:n_time){
    mu[n] = exp( alpha + beta * clim_means[n] );
  }
  y ~ gamma(y_sd, (y_sd ./ mu) );
  
}

generated quantities {
  vector[n_time] log_lik;
  
  for (n in 1:n_time)
    log_lik[n] = gamma_lpdf(y[n] | y_sd, (y_sd / exp(alpha + clim_means[n] * beta)) );

}
