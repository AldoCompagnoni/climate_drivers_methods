
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
  int K;
  int M;
 
  vector[n_train] y_train;
  vector[n_test]  y_test;

  matrix[M, n_train] clim1_train;
  matrix[M, n_train] clim2_train;
  matrix[M, n_train] clim3_train;
  matrix[M, n_test] clim1_test;
  matrix[M, n_test] clim2_test;
  matrix[M, n_test] clim3_test;
}

parameters {
  real<lower=-2,upper=2> shape;
  real<lower=0,upper=M> scale;
  real<lower=-(M*2),upper=M*2> loc;
  simplex[K] theta_y; 
  real alpha;
  real beta;
  real<lower=0> y_sd;
}

transformed parameters {
  
  vector[n_train] x;
  vector[M] sens_m;
  matrix[K,n_train] x_m;
  
  for(i in 1:M)
    sens_m[i] = dgev(i, loc, scale, shape);
  
  sens_m = sens_m / sum(sens_m);
  
  for(i in 1:n_train) {
    x_m[1,i] = sum(sens_m .* clim1_train[,i]); 
    x_m[2,i] = sum(sens_m .* clim2_train[,i]);
    x_m[3,i] = sum(sens_m .* clim3_train[,i]);
  }

  for(i in 1:n_train)
    x[i] = sum(theta_y .* x_m[,i]);
}

model {
  // place holder  
  vector[n_train] mu;    // transformed linear predictor for mean of beta distribution
  
  // likelihood
  for(n in 1:n_train)
    mu[n] = exp(alpha + x[n] * beta);
    
  y_train ~ gamma(y_sd, y_sd ./ mu);
}

generated quantities {
  vector[n_train]  log_lik;
  vector[n_test]   log_lik_test;
  vector[n_test]   pred_y;    // transformed linear predictor for mean of beta distribution
  vector[n_test]   pred_x;
  matrix[K,n_test] pred_x_m;

  for (n in 1:n_train)
    log_lik[n] = gamma_lpdf(y_train[n] | y_sd, (y_sd / exp(alpha + x[n] * beta)) );

  // out of sample prediction
  for(n in 1:n_test){
    pred_x_m[1,n] = sum(sens_m  .* clim1_test[,n]); 
    pred_x_m[2,n] = sum(sens_m  .* clim2_test[,n]);
    pred_x_m[3,n] = sum(sens_m  .* clim3_test[,n]);
    pred_x[n]       = sum(theta_y .* pred_x_m[,n]);
    pred_y[n]       = exp(alpha + pred_x[n] * beta);
    log_lik_test[n] = gamma_lpdf(y_test[n] | y_sd, (y_sd / exp(alpha + pred_x[n] * beta)) );
  }
}
