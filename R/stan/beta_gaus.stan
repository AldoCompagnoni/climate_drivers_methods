
functions {
  real dnorm(real x, real mu, real sigma) {
    return((1 / sqrt(2 *pi()*pow(sigma, 2))) * exp(-((x-mu)^2) / (2*pow(sigma, 2))));
  }
}

data {
  
  int n_time;
  int n_lag;
  vector<lower=0, upper=1>[n_time] y;
  matrix[n_time, n_lag] clim;
  
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
    sens[i] = dnorm(i, sens_mu, sens_sd);
  
  sens = sens / sum(sens);
  
  for(n in 1:n_time)
    x[n] = sum(sens .* row(clim, n)); 
    
  for(n in 1:n_time){
    yhat[n] = inv_logit(alpha + x[n] * beta);
    A[n]    = yhat[n] * y_sd;
    B[n]    = (1.0 - yhat[n]) * y_sd;
  }

}

model {
  
  // hyper-parameters to weight climate effects 
  sens_sd ~ normal(0.5, 12);
  sens_mu ~ normal(18.5, 36);
  
  // parameters of data model
  alpha ~ normal(0,1);
  beta  ~ normal(0,1);
  y_sd  ~ gamma(1,1);

  // likelihood
  y ~ beta(A, B);
  
}

generated quantities {
  vector[n_time] log_lik;
  
  for (n in 1:n_time)
    log_lik[n] = beta_lpdf(y[n] | A[n], B[n]);

}
