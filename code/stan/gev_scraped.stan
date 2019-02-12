
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
  int<lower=0> N;
  vector[N] y; //vector containing data
}
parameters {
  real<lower=0> sigma;
  real xi;
  real<lower=0> mu;
}

model {
  sigma ~ normal(0.1,1);
  xi ~ normal(0.05,1);
  mu ~ normal(1.2,1);
  
  y ~ gev(mu,sigma,xi);
}

