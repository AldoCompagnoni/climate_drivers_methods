
data {
  int n_time;        // number of data points, length(y)
  int M;             // number of months-within-years
  int K;             // number of years
  vector[n_time] y;  // response
  matrix[n_time,M*K] clim; // matrix of climate covariates
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
  real<lower=0> eta; // maximum covariance for betas
  real<lower=0> rho; // degree of temporal autocorrelation for betas
  vector[M] z;       // unit normal prior for non-centered term
  simplex[K] theta_y;
}

transformed parameters {
  vector[n_time] yhat;
  vector[M] beta;
  vector[M*K] beta_wt;
  matrix[M,M] sigma_beta; // covariance matrix
  matrix[M,M] L;     // cholesky of covariance matrix
  
  // covariance
  sigma_beta = cov_exp_quad(month, eta, rho) + diag_matrix(rep_vector(0.001, M));
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
  z ~ normal(0, 1);
  
  rho ~ normal(0, 5);
  eta ~ normal(0, 1);
  
  alpha ~ normal(0, 5);
  mu_beta ~ normal(0, 5);
  y_sd ~ normal(0, 3);
  
  theta_y ~ dirichlet(rep_vector(1.0, K));
  
  y ~ normal(yhat, y_sd);
}

generated quantities {
  vector[n_time] log_lik;
  
  for (n in 1:n_time) {
    log_lik[n] = normal_lpdf(y[n] | yhat[n], y_sd);
  }
}
