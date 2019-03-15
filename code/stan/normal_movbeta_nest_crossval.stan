
data {
  int n_train;
  int n_test;
  int K;
  int M;
  vector[n_train] y_train;
  vector[n_test]  y_test;
  matrix[n_train, M*K] clim_train;
  matrix[n_test,  M*K] clim_test;
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
  vector[n_train] yhat;
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
  yhat = alpha + clim_train * beta_wt;
}

model {
  
  // Priors of Hyper-parmeters
  z ~ normal(0, 1);
  rho ~ normal(0, 5);
  eta ~ normal(0, 1);
  theta_y ~ dirichlet(rep_vector(1.0, K));
  
  // Priors of parmeters
  alpha ~ normal(0, 5);
  mu_beta ~ normal(0, 5);
  y_sd ~ normal(0, 3);
  
  y_train ~ normal(yhat, y_sd);
}

generated quantities {
  vector[n_train] log_lik;
  vector[n_test]  pred_y;
  vector[n_test]  log_lik_test;
  
  // in sample prediction
  for (n in 1:n_train)
    log_lik[n] = normal_lpdf(y_train[n] | yhat[n], y_sd);

  // out of sample prediction
  for(n in 1:n_test){
    pred_y[n]       = alpha + row(clim_test,n) * beta_wt;
    log_lik_test[n] = normal_lpdf(y_test[n] | pred_y[n], y_sd);
  }
}
