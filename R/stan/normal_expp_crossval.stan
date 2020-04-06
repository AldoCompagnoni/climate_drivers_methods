
functions {
  real dexppow(real x, real mu, real sigma, real beta) {
    return((beta / (2 * sigma * tgamma(1.0/beta)) ) * exp(-(fabs(x - mu)/sigma)^beta));
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
  real expp_beta;
}

parameters {
  real<lower=0,upper=n_lag> sens_mu;
  real<lower=1> sens_sd;
  real alpha;
  real beta;
  real<lower=0> y_sd;
}

transformed parameters {
  vector[n_train] x;
  row_vector[n_lag] sens;
  
  for(i in 1:n_lag)
    sens[i] = dexppow(i, sens_mu, sens_sd, expp_beta);
  
  sens = sens / sum(sens);
  
  for(i in 1:n_train)
    x[i] = sum(sens .* row(clim_train, i));
}

model {
  
  // priors
  alpha ~ normal(0,1);
  beta  ~ normal(0,1);
  y_sd  ~ gamma(1,1);
  sens_sd ~ normal(0.5, 12);
  sens_mu ~ normal(18.5, 36); 

  // model
  y_train ~ normal(alpha + beta * x, y_sd);
}

generated quantities {
  vector[n_train] log_lik;
  vector[n_test]  pred_y;
  vector[n_test]  log_lik_test;
  real            pred_x;
  
  for (n in 1:n_train)
    log_lik[n] = normal_lpdf(y_train[n] | alpha + beta * x[n], y_sd);
    
  // out of sample prediction
  for(n in 1:n_test){
    pred_x          = sum(sens .* row(clim_test,n) );
    pred_y[n]       = alpha + beta * pred_x;
    log_lik_test[n] = normal_lpdf(y_test[n] | alpha + beta * pred_x, y_sd);
  }

}
