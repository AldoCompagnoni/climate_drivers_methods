
functions {
  real dnorm(real x, real mu, real sigma) {
    return((1 / sqrt(2 *pi()*pow(sigma, 2))) * exp(-((x-mu)^2) / (2*pow(sigma, 2))));
  }
}

data {
  
  int n_train;
  int n_test;
  int n_lag;
  vector[n_train] y_train;
  vector[n_test]  y_test;
  matrix[n_train, n_lag] clim_train;
  matrix[n_test,  n_lag] clim_test;
  
}

parameters {
  
  real<lower=0,upper=n_lag> sens_mu;
  real<lower=0,upper=n_lag*2> sens_sd;
  real alpha;
  real beta;
  real<lower=0> y_sd;
  
}

transformed parameters{
  
  // transformed parameters for moving window
  vector[n_train] x;
  row_vector[n_lag] sens;
  
  for(i in 1:n_lag) { sens[i] = dnorm(i, sens_mu, sens_sd); }
  
  sens = sens / sum(sens);
  
  for(n in 1:n_train){ x[n] = sum(sens .* row(clim_train, n)); }
    
}

model {
  
  // place holder  
  vector[n_train] mu;    // transformed linear predictor for mean of beta distribution
  
  // likelihood
  for(n in 1:n_train)
    mu[n] = exp(alpha + x[n] * beta);
    
  y_train ~ gamma(y_sd, y_sd ./ mu);
  
}

generated quantities {
  
  vector[n_train] log_lik;
  vector[n_test]  log_lik_test;
  vector[n_test]  pred_y;    // transformed linear predictor for mean of beta distribution
  vector[n_test]  pred_x;

  for (n in 1:n_train)
    log_lik[n] = gamma_lpdf(y_train[n] | y_sd, (y_sd / exp(alpha + x[n] * beta)) );

  for(n in 1:n_test){
    pred_x[n]       = sum(sens .* row(clim_test,n));
    pred_y[n]       = exp(alpha + pred_x[n] * beta);
    log_lik_test[n] = gamma_lpdf(y_test[n] | y_sd, (y_sd / exp(alpha + pred_x[n] * beta)) );
  }

}
