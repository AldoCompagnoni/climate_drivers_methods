  
data {
  int n_time;
  int M;        // number of months-within-years
  int K;        // number of years
  vector[n_time] y;
  matrix[n_time,M*K] clim;  // matrix of climate covariates
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
  
  // params for random beta
  vector[n_time] yhat;
  vector[M] beta;
  vector[M*K] beta_wt;
  
  // non-centered parameterization for beta
  beta = mu_beta + sigma_beta * z;
  
  // apply yearly weights
  beta_wt[m1] = theta_y[1] * beta;
  beta_wt[m2] = theta_y[2] * beta;
  beta_wt[m3] = theta_y[3] * beta;
  
  // linear predictor
  yhat = exp(alpha + clim * beta_wt);
  
}

model {
  // hyper-parameters to weight climate effects
  z ~ normal(0, 1);
  mu_beta ~ normal(0, 3);
  sigma_beta ~ normal(0, 3);
  theta_y ~ dirichlet(rep_vector(1.0, K));
  
  // priors
  alpha ~ normal(0,1);
  y_sd  ~ gamma(1,1); 

  y ~ gamma(y_sd, y_sd ./ yhat);
}

generated quantities {
  vector[n_time] log_lik;
  
  for (n in 1:n_time)
    log_lik[n] = gamma_lpdf(y[n] | y_sd, (y_sd / yhat[n]) );
}
