
data {
  int<lower=0> n_time; // number of data points, length(y)
  int<lower=0> n_lag; // number of monthly lags
  vector[n_time] y;    // response
  matrix[n_time,n_lag] clim;  // matrix of climate covariates
}

parameters {
  real alpha;
  real<lower=0> phi;
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
  
  // parameters for ridge regression
  real<lower=0> sigma2[n_time];
  real<lower=0> tau2[n_time];     

  // linear predictor
  mu = alpha + clim * beta;
  
  // beta reparameterization
  for(n in 1:n_time){
    // Parameters A and B 
    yhat[n]   = inv_logit(mu[n]);
    A[n]      = yhat[n] * phi;
    B[n]      = (1.0 - yhat[n]) * phi;
    
    // Ridge regression parameters
    sigma2[n] = (yhat[n] * (1 - yhat[n])) / (1 + phi);
    tau2[n]   = sigma2[n] / lambda;
 }
  
 
}

model {
  
  // prior regression coefficients: ridge
  beta   ~ normal(0, sqrt(tau2) );
  lambda ~ cauchy(0, 1);

  // parameters of data model
  alpha  ~ normal(0, 5);
  // Prior when running Beta(link='logit') in BRMS
  phi    ~ gamma(0.001, 0.001); 
  
  y ~ beta(A, B);
}
