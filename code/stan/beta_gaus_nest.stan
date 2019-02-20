
functions {
  real dnorm(real x, real mu, real sigma) {
    return((1 / sqrt(2 *pi()*pow(sigma, 2))) * exp(-((x-mu)^2) / (2*pow(sigma, 2))));
  }
}

data {
  int n_time;
  int M; // Number of months
  int K; // Number of years
  vector[n_time] y;
  matrix[M,n_time] clim1;
  matrix[M,n_time] clim2;
  matrix[M,n_time] clim3;
}

parameters {
  real<lower=0,upper=M> sens_mu;
  real<lower=1> sens_sd;
  simplex[K] theta_y;
  real alpha;
  real beta;
  real<lower=0> y_sd;
}

transformed parameters {
  
  // transformed parameters for beta binomial regression
  real<lower=0,upper=1> mu[n_time]; // transformed linear predictor for mean of beta distribution
  real<lower=0> A[n_time];          // parameter for beta distn
  real<lower=0> B[n_time];          // parameter for beta distn

  vector[n_time] x;
  vector[M] sens_m;
  matrix[K,n_time] x_m;
  
  for(n in 1:M)
    sens_m[n] = dnorm(n, sens_mu, sens_sd);
  
  sens_m = sens_m / sum(sens_m);
  
  for(n in 1:n_time) {
    x_m[1,n] = sum(sens_m .* clim1[,n]); 
    x_m[2,n] = sum(sens_m .* clim2[,n]);
    x_m[3,n] = sum(sens_m .* clim3[,n]);
  }

  for(n in 1:n_time)
    x[n] = sum(theta_y .* x_m[,n]);
    
  for(n in 1:n_time){
    mu[n]  = inv_logit(alpha + x[n] * beta);
    A[n]   = mu[n] * y_sd;
    B[n]   = (1.0 - mu[n]) * y_sd;
  }
}

model {
  
  // priors
  alpha ~ normal(0,1);
  beta  ~ normal(0,1);
  y_sd  ~ gamma(1,1);
  sens_sd ~ normal(0.5, 12);
  sens_mu ~ normal(6.5, 12);
  theta_y ~ dirichlet(rep_vector(1.0, K));
  
  // likelihood
  y ~ beta(A, B);
  
}

generated quantities {
  vector[n_time] log_lik;
  
  for (n in 1:n_time)
    log_lik[n] = beta_lpdf(y[n] | A[n], B[n]);
}
