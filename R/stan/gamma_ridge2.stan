  
data {
  int<lower=0> n_time;
  int<lower=0> n_lag;       // number of monthly lags 
  vector[n_time] y;
  matrix[n_time,n_lag] clim;  // matrix of climate covariates
  real<lower=0> sigma;
}

parameters {
  real alpha;
  vector[n_lag] beta; 
  //hyperparameters prior
  real<lower = 0> lambda; //penalty parameter
  real<lower = 0> y_sd;    // error sd 
}

transformed parameters {
  
  // place holder  
  vector[n_time] yhat;     // transf. lin. pred. for mean
  real<lower=0> A[n_time]; // parameter for beta distn
  real<lower=0> B[n_time]; // parameter for beta distn
  // real<lower=0> tau2;      // prior variance
  
  yhat = exp(alpha + clim * beta);
  
  // beta reparameterization
  for(n in 1:n_time){
    // Parameters A and B 
    A[n]      = y_sd;
    B[n]      = y_sd / yhat[n];
  }
  
  // // Ridge regression parameters
  // tau2   = sigma2 / lambda;
  // y_sd   = sqrt(sigma2);
  
}

model {
  
  // prior regression coefficients: ridge
  beta ~ double_exponential(0, sigma/lambda);
	lambda ~ cauchy(0, 1);
 
 //priors nuisance parameters: uniform on log(sigma^2) & mu
	// target += -2 * log(sigma); 
 
  // priors
  y_sd  ~ gamma(0.01,0.01); // this is a variance

  y ~ gamma(A, B);
  
}
