
data {
  int<lower=0> n_time;  // number of data points, length(y)
  int<lower=0> n_lag;       // number of monthly lags
  vector[n_time] y;          // response
  matrix[n_time,n_lag] clim;   // matrix of climate covariates
}

transformed data {
  real month[n_lag]; // vector of months, to create distance matrix
  for(k in 1:n_lag) month[k] = k;
}

parameters {
  real alpha;
  real mu_beta;
  real<lower=0> y_sd;
  real<lower=0> eta;  // maximum covariance for betas
  real<lower=0> rho;  // degree of temporal autocorrelation for betas
  vector[n_lag] z;    // unit normal prior for non-centered term
}

transformed parameters {
 
  // transformed parameters for beta binomial regression
  real<lower=0,upper=1> yhat[n_time]; // transf. lin. pred. for mean of beta distribution
  real<lower=0> A[n_time];          // parameter for beta distn
  real<lower=0> B[n_time];          // parameter for beta distn
 
  // params for random beta
  vector[n_time] mu;
  vector[n_lag] beta;
  matrix[n_lag,n_lag] sigma_beta; // covariance matrix
  matrix[n_lag,n_lag] L;     // cholesky of covariance matrix
  
  // covariance
  sigma_beta = cov_exp_quad(month, eta, rho) + diag_matrix(rep_vector(0.001, n_lag));
  L = cholesky_decompose(sigma_beta);
  
  // non-centered parameterization for beta
  beta = mu_beta + L * z;
  
  // linear predictor
  mu = alpha + clim * beta;
  
  // beta reparameterization
  for(n in 1:n_time){
    yhat[n]  = inv_logit(mu[n]);
    A[n]   = yhat[n] * y_sd;
    B[n]   = (1.0 - yhat[n]) * y_sd;
  }
  
}

model {
  
  // hyper-parameters to weight climate effects
  z ~ normal(0, 1);
  rho ~ normal(0, 5);
  eta ~ normal(0, 1);
  
  // parameters of data model
  alpha ~ normal(0, 5);
  mu_beta ~ normal(0, 5);
  y_sd ~ gamma(1,1); 

  // model
  y ~ beta(A, B);
}

generated quantities {
  vector[n_time] log_lik;
  
  for (n in 1:n_time) {
    log_lik[n] = beta_lpdf(y[n] | A[n], B[n]);
  }
}
