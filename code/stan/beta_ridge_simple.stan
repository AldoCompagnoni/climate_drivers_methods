
// https://stats.stackexchange.com/questions/12232/calculating-the-parameters-of-a-beta-distribution-using-the-mean-and-variance

data {
  int<lower=0> n_time; // number of data points, length(y)
  int<lower=0> n_lag; // number of monthly lags
  vector[n_time] y;    // response
  matrix[n_time,n_lag] clim;  // matrix of climate covariates
}

parameters {
  real alpha;
  real<lower=0,upper=0.25> sigma2;
  vector[n_lag] beta;     // regression parameters
  //hyperparameters prior
  real<lower = 0> lambda; //penalty parameter
  
}

transformed parameters {
  
  // transformed parameters for beta binomial regression
  real<lower=0,upper=1> yhat[n_time]; // transf. lin. pred. for mean of beta distribution
  real<lower=0> A[n_time];    // parameter for beta distn
  real<lower=0> B[n_time];    // parameter for beta distn
  vector[n_time] mu;          // placeholder to create parameters
  real<lower = 0> y_sd;    // error sd 
  
  // parameters for ridge regression
  real<lower=0> tau2;     

  // linear predictor
  mu = alpha + clim * beta;
  
  // beta reparameterization
  for(n in 1:n_time){
    
    // Parameters A and B 
    yhat[n]   = inv_logit(mu[n]);
    A[n]      = (yhat[n] * ( sigma2 + square(yhat[n]) - yhat[n]) ) / sigma2;
    B[n]      = ( (sigma2 + square(yhat[n]) - yhat[n]) * (yhat[n] - 1)) / sigma2;
  }
  
  // Ridge regression parameters
  tau2   = sigma2 / lambda;
  y_sd   = sqrt(sigma2);
 
}

model {
  
  // prior regression coefficients: ridge
  beta   ~ normal(0, sqrt(tau2) );
  lambda ~ cauchy(0, 1);

  //priors nuisance parameters: uniform on log(sigma^2) & mu
	target += -2 * log(y_sd); 

  // Prior on variance 
  sigma2 ~ inv_gamma(0.001, 0.001); 
  
  y ~ beta(A, B);
}
