

functions {

  real gev_lpdf (vector y, real mu, real sigma, real xi){
    vector[rows(y)] z;
    real inv_xi;
    real neg_inv_xi; 
    real inv_xi_p1; 
    vector[rows(y)+1] lp;
    int N;

    N = rows(y);

    z = 1 + (y - mu) * xi / sigma;
    inv_xi = 1/xi;
    neg_inv_xi = -inv_xi; 
    inv_xi_p1 = 1 + inv_xi; 

    for(n in 1:N){
      lp[n] = inv_xi_p1*log(z[n]) + pow(z[n],neg_inv_xi);
    }

    lp[N+1] = N * log(sigma);
    #print("2: ",lp);
    return -sum(lp);  
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
  real<lower=-0.5,upper=0.5> xi;
  real<lower=0> sigma;
  // location has upper/lower bounds depending on the value of xi
  real<lower=if_else( xi > 0, min_y, negative_infinity()),
       upper=if_else( xi > 0, positive_infinity(), max_y )> mu;
}
model {

  # priors
  sigma ~ normal(sd_y, 100);
  #mu ~ uniform
  #xi ~ uniform(-.5,.5);

  increment_log_prob(gev_lpdf(y, mu, sigma, xi)); 
}