functions{
  real gev_lpdf(vector y,real mu,real sigma,real xi){
  vector[rows(y)] z;
  vector[rows(y)] logp;
  vector[rows(y)] zi;
  
  z = 1 + xi*(y - mu)/sigma;
  for(i in 1:rows(y)){
      zi[i] = z[i]^(-1/xi);
  }
  
  logp = log(sigma) + (1 + 1/xi)*log(z) + zi;
  return -sum(logp); 
  }
}

data {
 int<lower=0> N;
 vector[N] y; //vector containing data
}

parameters {
  real<lower=0> sigma;
  real<lower=-2,upper=2> xi;
  real<lower=0> mu;
}

model {
  sigma ~ normal(1.05,1);
  xi ~ normal(0.05,1);
  mu ~ normal(1.2,1);
  
  y ~ gev(mu,sigma,xi);
}
