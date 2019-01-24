
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
  vector<lower=0>[n_time] y;
  matrix[n_time, n_lag] clim;
  
}

parameters {
  real<lower=-2,upper=2> shape;
  real<lower=0,upper=n_lag> scale;
  real<lower=-(n_lag*2),upper=n_lag*2> loc;
  real alpha;
  real beta;
  real<lower=0> y_sd;  
}

transformed parameters{
  
  // transformed parameters for moving window
  vector[n_time] x;
  row_vector[n_lag] sens;
  
  for(i in 1:n_lag)
    sens[i] = dgev(i, loc, scale, shape);
  
  sens = sens / sum(sens);
  
  for(n in 1:n_time)
    x[n] = sum(sens .* row(clim, n));
    
}

model {
  
  // place holder  
  vector[n_time] mu;    // transformed linear predictor for mean of beta distribution
  
  // likelihood
  for(n in 1:n_time)
    mu[n] = exp(alpha + x[n] * beta);
    
  y ~ gamma(y_sd, y_sd ./ mu);
  
}

generated quantities {
  vector[n_time] log_lik;
  
  for (n in 1:n_time)
    log_lik[n] = gamma_lpdf(y[n] | y_sd, (y_sd / exp(alpha + x[n] * beta)) );

}
