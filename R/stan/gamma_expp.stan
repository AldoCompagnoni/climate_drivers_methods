
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
  
  // transformed parameters for moving window
  vector[n_time] x;
  row_vector[n_lag] sens;
  vector[n_time] yhat;
  
  for(i in 1:n_lag)
    sens[i] = dexppow(i, sens_mu, sens_sd, expp_beta);
  
  sens = sens / sum(sens);
  
  for(n in 1:n_time)
    x[n] = sum(sens .* row(clim, n));
   
  yhat = exp(alpha + x * beta);
  
}

model {
  
  // hyper-parameters to weight climate effects
  sens_sd ~ normal(0.5, 12);
  sens_mu ~ normal(6.5, 12);

  // priors
  alpha ~ normal(0,1);
  beta  ~ normal(0,1);
  y_sd  ~ gamma(1,1);
  
  y ~ gamma(y_sd, y_sd ./ yhat);
  
}

generated quantities {
  vector[n_time] log_lik;
  
  for (n in 1:n_time)
    log_lik[n] = gamma_lpdf(y[n] | y_sd, (y_sd / exp(alpha + x[n] * beta)) );

}
