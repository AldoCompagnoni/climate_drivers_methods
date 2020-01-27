  
data {
  int<lower=0> n_time;
  int<lower=0> n_lag;       // number of monthly lags 
  vector[n_time] y;
  matrix[n_time,n_lag] clim;  // matrix of climate covariates
}

parameters {
  real alpha;
  vector[n_lag] beta; 
  real<lower=0> sigma2;
  //hyperparameters prior
  real<lower = 0> lambda; //penalty parameter
}

transformed parameters {
  
  // place holder  
  vector[n_time] yhat;     // transf. lin. pred. for mean
  real<lower=0> A[n_time]; // parameter for beta distn
  real<lower=0> B[n_time]; // parameter for beta distn
  real<lower=0> tau2;      // prior variance
  real<lower = 0> y_sd;    // error sd 
  
  yhat = exp(alpha + clim * beta);
  
  // beta reparameterization
  for(n in 1:n_time){
    // Parameters A and B 
    A[n]      = square(yhat[n]) / sigma2;
    B[n]      = yhat[n] / sigma2;
  }
  
  // Ridge regression parameters
  tau2   = sigma2 / lambda;
  y_sd   = sqrt(sigma2);
  
}

model {
  
  // prior regression coefficients: ridge
  beta   ~ normal(0, sqrt(tau2));
  lambda ~ cauchy(0, 1);
 
  //priors nuisance parameters: uniform on log(sigma^2) & mu
	target += -2 * log(y_sd); 
 
  // priors
  sigma2  ~ inv_gamma(0.01,0.01); // this is a variance

  y ~ gamma(A, B);
  
}
