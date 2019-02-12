
functions{
  // GEV log pdf
  real dgev(real y,real mu,real sigma,real xi){
      real z;
      real logp;
      real zi;
      
      z  = 1 + xi*(y - mu)/sigma;
      zi = z^(-1/xi);
      
      logp = log(sigma) + (1 + 1/xi)*log(z) + zi;
      return exp(logp); 
  }
}

data {
  int n_time;
  int n_lag;
  vector[n_time] y;
  matrix[n_time, n_lag] clim;
}
  
parameters {
  real<lower=-2,upper=2> xi;
  real<lower=1> sigma; //
  real<lower=0,upper=n_lag> mu;
  real alpha;
  real beta;
  real<lower=0> y_sd;
}

transformed parameters {
  vector[n_time] x;
  row_vector[n_lag] sens;
  
  for(i in 1:n_lag)
    sens[i] = dgev(i, mu, sigma, xi);
  
  sens = sens / sum(sens);
  
  for(i in 1:n_time)
    x[i] = sum(sens .* row(clim, i));
}

model {
  
  alpha ~ normal(0,1);
  beta  ~ normal(0,1);
  y_sd  ~ gamma(1,1);
  sigma ~ normal(1, 12);
  
  // model
  y ~ normal(alpha + beta * x, y_sd);
}

// generated quantities {
//   vector[n_time] log_lik;
//   
//   for (n in 1:n_time)
//     log_lik[n] = normal_lpdf(y[n] | alpha + beta * x[n], y_sd);
// }
