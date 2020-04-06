
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
  int n_train;
  int n_test;
  int n_lag;
  vector[n_train] y_train;
  vector[n_test]  y_test;
  // real y_test;
  matrix[n_train, n_lag] clim_train;
  matrix[n_test,  n_lag] clim_test;
}

parameters {
  real<lower=-2,upper=2> shape;
  real<lower=0,upper=n_lag> scale;
  real<lower=-(n_lag*2),upper=n_lag*2> loc;
  real alpha;
  real beta;
  real<lower=0> y_sd;
}

transformed parameters {
  vector[n_train] x;
  row_vector[n_lag] sens;
  
  for(i in 1:n_lag)
    sens[i] = dgev(i, loc, scale, shape);
  
  sens = sens / sum(sens);
  
  for(i in 1:n_train)
    x[i] = sum(sens .* row(clim_train, i));
}

model {
  // model
  y_train ~ normal(alpha + beta * x, y_sd);
}

generated quantities {
  vector[n_train] log_lik;
  vector[n_test]  pred_y;
  vector[n_test]  log_lik_test;
  real            pred_x;
  
  for (tr in 1:n_train)
    log_lik[tr] = normal_lpdf(y_train[tr] | alpha + beta * x[tr], y_sd);
  
  // out of sample prediction
  for(te in 1:n_test){
    pred_x           = sum(sens .* row(clim_test,te));
    pred_y[te]       = alpha + beta * pred_x;
    log_lik_test[te] = normal_lpdf(y_test[te] | alpha + beta * pred_x, y_sd);
  }
  
}
