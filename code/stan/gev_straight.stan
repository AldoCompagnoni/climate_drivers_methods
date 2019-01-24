
functions {

  real gev_lpdf(real y, real loc, real scale, real shape){
    
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

    return lp;  
  }

}

data {
  int<lower=0> N;
  vector[N] y;
}

transformed data {
  real min_y;
  real max_y;
  real sd_y;
  min_y = min(y);
  max_y = max(y);
  sd_y = sd(y);
}

parameters {
  real shape;
  real<lower=0> scale;
  // location has upper/lower bounds depending on the value of xi
  real<lower=if_else( shape < 0, min_y + scale / shape, negative_infinity() ),
       upper=if_else( shape > 0, positive_infinity(), max_y + scale / shape )> loc;
}

model {

  //priors
  for(ii in 1:N)
    target += gev_lpdf(y[ii]| loc, scale, shape);
    
}
