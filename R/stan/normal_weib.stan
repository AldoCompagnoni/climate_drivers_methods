
functions {
  real my_weib(real y, real lambda, real sigma) {
    real pp;
    real out;
    real one;
    real two;
    real three;
    
    one   = lambda / sigma;
    two   = pow((y / sigma),(lambda-1));
    three = -pow((y/sigma),lambda);
    
    pp    = one * two * exp(three);
    return pp;
  }
}

data {
  int n_time;
  int n_lag;
  vector[n_time] y;
  matrix[n_time, n_lag] clim;
}

parameters {
  real<lower=0> lambda;
  real<lower=0,upper=50> sigma;
  real alpha;
  real beta;
  real<lower=0> y_sd;
}

transformed parameters {
  vector[n_time] x;
  row_vector[n_lag] sens;
  
  for(i in 1:n_lag)
    sens[i] = my_weib(i, lambda, sigma);
  
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

