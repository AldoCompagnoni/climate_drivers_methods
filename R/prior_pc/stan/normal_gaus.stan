
functions {
  real dnorm(real x, real mu, real sigma) {
    return((1 / sqrt(2 *pi()*pow(sigma, 2))) * exp(-((x-mu)^2) / (2*pow(sigma, 2))));
  }
}

data {
  int n_time;
  int n_lag;
  matrix[n_time, n_lag] clim;
  real gamma_shape;
}

parameters {
  real<lower=0,upper=n_lag> sens_mu;
  real<lower=0> sens_sd; //
  real alpha;
  real beta;
  real<lower=0> y_sd;
}

transformed parameters {
  vector[n_time] x;
  row_vector[n_lag] sens;
  vector[n_time] yhat;
  
  for(i in 1:n_lag)
    sens[i] = dnorm(i, sens_mu, sens_sd);
  
  sens = sens / sum(sens);
  
  for(i in 1:n_time)
    x[i] = sum(sens .* row(clim, i));
  
  yhat = alpha + beta * x;  
    
}

model {
  
  // hyper-parameters to weight climate effects
  sens_sd ~ normal(0.5, 12);
  sens_mu ~ normal(6.5, 12);

  // priors
  alpha ~ normal(0,1);
  beta  ~ normal(0,1);
  y_sd  ~ gamma(gamma_shape,0.01); 
 
}

generated quantities {
  
  real y_sim[n_time]; 
  
  y_sim = normal_rng(yhat, y_sd);

}

