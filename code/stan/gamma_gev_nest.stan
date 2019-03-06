
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
  int M; // Number of months
  int K; // Number of years
  vector[n_time] y;
  matrix[M,n_time] clim1;
  matrix[M,n_time] clim2;
  matrix[M,n_time] clim3;
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
  vector[n_time] x;
  vector[M] sens_m;
  matrix[K,n_time] x_m;
  
  for(i in 1:M)
    sens_m[i] = dgev(i, loc, scale, shape);
  
  sens_m = sens_m / sum(sens_m);
  
  for(i in 1:n_time) {
    x_m[1,i] = sum(sens_m .* clim1[,i]);
    x_m[2,i] = sum(sens_m .* clim2[,i]);
    x_m[3,i] = sum(sens_m .* clim3[,i]);
  }

  for(i in 1:n_time)
    x[i] = sum(theta_y .* x_m[,i]);

}

model {
  // place holder  
  vector[n_train] mu; // transf. lin. pred. for mean

  // priors
  alpha ~ normal(0,1);
  beta  ~ normal(0,1);
  y_sd  ~ gamma(1,1); 
  
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
