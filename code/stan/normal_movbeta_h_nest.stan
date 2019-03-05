
data {
  int n_time;        // number of data points, length(y)
  int M;        // number of months-within-years
  int K;        // number of years
  vector[n_time] y;  // response
  matrix[n_time,M*K] clim; // matrix of climate covariates
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
  vector[M] z;   // unit normal prior for non-centered term
  simplex[K] theta_y;
}

transformed parameters {
  vector[n_time] yhat;
  vector[M] beta;
  vector[M*K] beta_wt;
  
  // non-centered parameterization
  beta = mu_beta + sigma_beta * z;
  
  // apply yearly weights
  beta_wt[m1] = theta_y[1] * beta;
  beta_wt[m2] = theta_y[2] * beta;
  beta_wt[m3] = theta_y[3] * beta;
  
  // linear predictor
  yhat = alpha + clim * beta_wt;
}

model {
  
  // hyper-parameters to weight climate effects
  theta_y ~ dirichlet(rep_vector(1.0, K));
  z ~ normal(0, 1);
  mu_beta ~ normal(0, 3);
  sigma_beta ~ normal(0, 3);
  
  // parameters of data model
  alpha ~ normal(0, 5);
  y_sd ~ normal(0, 3);
  
  y ~ normal(yhat, y_sd);
}

generated quantities {
  vector[n_time] log_lik;
  
  for (n in 1:n_time) {
    log_lik[n] = normal_lpdf(y[n] | yhat[n], y_sd);
  }
}
