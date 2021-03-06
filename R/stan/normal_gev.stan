
functions{
  real dgev(real y, real loc, real scale, real shape){
    
    real inv_scale;
    real neg_inv_shape; 
    real inv_shape_p1; 
    real x;
    real xx;
    real lp;

    x             = (y - loc)/scale;
    inv_scale     = 1/scale;
    neg_inv_shape = -1/shape;
    inv_shape_p1  = (1/shape) + 1;

    xx = 1 + shape * x;
    lp = log(inv_scale) - pow(xx,neg_inv_shape) - (inv_shape_p1 * log(xx));
     
    return exp(lp);
    
  }
  
}

data {
  int n_time;
  int n_lag;
  vector[n_time] y;
  matrix[n_time, n_lag] clim;
  real shape;
}

parameters {
  // real<lower=-0.5,upper=0.5> shape;
  real<lower=1> scale; //
  real<lower=0,upper=n_lag> loc;
  real alpha;
  real beta;
  real<lower=0> y_sd;
}

transformed parameters {
  vector[n_time] x;
  row_vector[n_lag] sens;
  
  for(i in 1:n_lag)
    sens[i] = dgev(i, loc, scale, shape);
  
  sens = sens / sum(sens);
  
  for(i in 1:n_time)
    x[i] = sum(sens .* row(clim, i));
}

model {
  
  alpha ~ normal(0,1);
  beta  ~ normal(0,1);
  y_sd  ~ gamma(1,1);
  // scale ~ normal(0.5, 12);
  
  // model
  y ~ normal(alpha + beta * x, y_sd);
}

// generated quantities {
//   vector[n_time] log_lik;
//   
//   for (n in 1:n_time)
//     log_lik[n] = normal_lpdf(y[n] | alpha + beta * x[n], y_sd);
// }
