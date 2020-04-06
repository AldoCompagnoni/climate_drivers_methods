
functions {
  real dexppow(real x, real mu, real sigma, real beta) {
    return((beta / (2 * sigma * tgamma(1.0/beta)) ) * exp(-(fabs(x - mu)/sigma)^beta));
  }
}

data {
  
  int n_time;
  int n_lag;
  vector[n_time] y;
  matrix[n_time, n_lag] clim;
  real expp_beta;
  
}

parameters {
  
  real<lower=0,upper=n_lag> sens_mu;
  real<lower=1> sens_sd;
  real alpha;
  real beta;
  real<lower=0> y_sd;
  
}

transformed parameters{
  
  // transformed parameters for beta binomial regression
  real<lower=0,upper=1> yhat[n_time];    // transformed linear predictor for mean of beta distribution
  real<lower=0> A[n_time];               // parameter for beta distn
  real<lower=0> B[n_time];               // parameter for beta distn
  
  // transformed parameters for moving window
  vector[n_time] x;
  row_vector[n_lag] sens;
  
  
  for(i in 1:n_lag)
    sens[i] = dexppow(i, sens_mu, sens_sd, expp_beta);
  
  sens = sens / sum(sens);
  
  for(i in 1:n_time)
    x[i] = sum(sens .* row(clim, i));
    
  for(n in 1:n_time){
    yhat[n] = inv_logit(alpha + x[n] * beta);
    A[n]    = yhat[n] * y_sd;
    B[n]    = (1.0 - yhat[n]) * y_sd;
  }

}

model {
  
  // priors
  alpha ~ normal(0,1);
  beta  ~ normal(0,1);
  y_sd  ~ gamma(1,1);
  sens_sd ~ normal(0.5, 12);
  sens_mu ~ normal(18.5, 36);
  
  // likelihood
  y ~ beta(A, B);
  
}

generated quantities {
  vector[n_time] log_lik;
  
  for (n in 1:n_time)
    log_lik[n] = beta_lpdf(y[n] | A[n], B[n]);

}
