
functions {
  real dnorm(real x, real alphaw, real sigma) {
    return( (alphaw / sigma) * (x / sigma)^(alphaw-1) * exp(-((x/sigma)^alphaw)) );
  }
}


data {
  int n_time;
  int n_lag;
  vector[n_time] y;
  matrix[n_time, n_lag] clim;
}

parameters {
  real<lower=0> sens_mu;
  real<lower=0> sens_sd;
  real alpha;
  real beta;
  real<lower=0> y_sd;
}

transformed parameters {
  vector[n_time] x;
  row_vector[n_lag] sens;
  
  for(i in 1:n_lag)
    sens[i] = dnorm(i, sens_mu, sens_sd);
  
  sens = sens / sum(sens);
  
  for(i in 1:n_time)
    x[i] = sum(sens .* row(clim, i));
}

model {
  // model
  y ~ normal(alpha + beta * x, y_sd);
}

// generated quantities {
//   vector[n_time] log_lik;
//   
//   for (n in 1:n_time)
//     log_lik[n] = normal_lpdf(y[n] | alpha + beta * x[n], y_sd);
// }

