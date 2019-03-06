
data {
  int<lower=0> n_time;  // number of data points, length(y)
  int<lower=0> n_lag;       // number of monthly lags
  vector[n_time] y;
  matrix[n_time,n_lag] clim;   // matrix of climate covariates
}

transformed data {
  int m1[M];
  int m2[M];
  int m3[M];
  real month[M];
  
  // month-year indices
  for (i in 1:M) {
    m1[i] = i;
    m2[i] = i + M;
    m3[i] = i + 2*M;
  }
  
  // vector of months, to create distance matrix
  for(m in 1:M) month[m] = m;
}

parameters {
  real alpha;
  real mu_beta;
  real<lower=0> y_sd;
  real<lower=0> eta;  // maximum covariance for betas
  real<lower=0> rho;  // degree of temporal autocorrelation for betas
  vector[n_lag] z;    // unit normal prior for non-centered term
  simplex[K] theta_y;
}

transformed parameters {
  
  // params for random beta
  vector[n_time] yhat;
  vector[n_lag] beta;
  vector[M*K] beta_wt;
  matrix[n_lag,n_lag] sigma_beta; // covariance matrix
  matrix[n_lag,n_lag] L;     // cholesky of covariance matrix
    
  // covariance
  sigma_beta = cov_exp_quad(month, eta, rho) + diag_matrix(rep_vector(0.001, n_lag));
  L = cholesky_decompose(sigma_beta);
  
  // non-centered parameterization for beta
  beta = mu_beta + L * z;
  
  // apply yearly weights
  beta_wt[m1] = theta_y[1] * beta;
  beta_wt[m2] = theta_y[2] * beta;
  beta_wt[m3] = theta_y[3] * beta;
  
  // linear predictor
  yhat = alpha + clim * beta_wt;
  
}

model {
  // place holder  
  vector[n_time] mu; // transf. lin. pred. for mean
  
  // hyper-parameters to weight climate effects
  z ~ normal(0, 1);
  rho ~ normal(0, 5);
  eta ~ normal(0, 1);
  theta_y ~ dirichlet(rep_vector(1.0, K));
  
  // parameters of data model
  alpha ~ normal(0, 5);
  mu_beta ~ normal(0, 5);
  y_sd ~ gamma(1,1); 

  for(n in 1:n_time)
    mu[n] = exp(yhat[n]);
    
  y ~ gamma(y_sd, y_sd ./ mu);
}

generated quantities {
  vector[n_time] log_lik;
  
  for (n in 1:n_time)
    log_lik[n] = gamma_lpdf(y[n] | y_sd, (y_sd / exp(yhat[n])) );
}
