
data {
  int n_time;
  vector[n_time] y;
  matrix[n_time,n_lag] clim;  // matrix of climate covariates
}

transformed data {
  int m1[M];   // month-year indices
  int m2[M];
  int m3[M];
  
  for (i in 1:M) {
    m1[i] = i;
    m2[i] = i + M;
    m3[i] = i + 2*M;
  }
}

parameters {
  real alpha;
  real<lower=0> y_sd;
  real mu_beta;
  real<lower=0> sigma_beta;
  simplex[K] theta_y;
}

transformed parameters {
  
  // params for random beta
  vector[n_time] yhat;
  vector[n_lag] beta;
  vector[M*K] beta_wt;
  
  // non-centered parameterization for beta
  beta = mu_beta + L * z;
  
  // apply yearly weights
  beta_wt[m1] = theta_y[1] * beta;
  beta_wt[m2] = theta_y[2] * beta;
  beta_wt[m3] = theta_y[3] * beta;
  
  // linear predictor
  yhat = alpha + clim * beta;
  
  
}

model {
  // place holder  
  vector[n_time] mu; // transf. lin. pred. for mean
  
  // hyper-parameters to weight climate effects
  z ~ normal(0, 1);
  mu_beta ~ normal(0, 3);
  sigma_beta ~ normal(0, 3);
  theta_y ~ dirichlet(rep_vector(1.0, K));
  
  // priors
  alpha ~ normal(0,1);
  y_sd  ~ gamma(1,1); 

  for(n in 1:n_time)
    mu[n] = exp(yhat[n]);
    
  y ~ gamma(y_sd, y_sd ./ mu);
}

generated quantities {
  vector[n_time] log_lik;
  
  for (n in 1:n_time)
    log_lik[n] = gamma_lpdf(y[n] | y_sd, (y_sd / exp(yhat[n])) );
}
