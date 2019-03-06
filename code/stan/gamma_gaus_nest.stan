
functions {
  real dnorm(real x, real mu, real sigma) {
    return((1 / sqrt(2 *pi()*pow(sigma, 2))) * exp(-((x-mu)^2) / (2*pow(sigma, 2))));
  }
}

data {
  int n_time;
  int M; // Number of months
  int K; // Number of years
  vector<lower=0>[n_time] y;
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

transformed parameters{
  
  // transformed parameters for moving window
  vector[n_time] x;
  vector[M] sens_m;
  matrix[K,n_time] x_m;
  
  for(i in 1:M)
    sens_m[i] = dnorm(i, sens_mu, sens_sd, expp_beta);
  
  sens_m = sens_m / sum(sens_m);
  
  for(i in 1:n_time) {
    x_m[1,i] = sum(sens_m .* clim1[,i]); 
    x_m[2,i] = sum(sens_m .* clim2[,i]);
    x_m[3,i] = sum(sens_m .* clim3[,i]);
  }

  for(i in 1:n_time)
    x[i] = sum(theta_y .* x_m[,i]);

}

model {
  // place holder  
  vector[n_time] mu; // transf. lin. pred. for mean
  
  // hyper-parameters to weight climate effects
  sens_sd ~ normal(0.5, 12);
  sens_mu ~ normal(6.5, 12);
  theta_y ~ dirichlet(rep_vector(1.0, K));
  
  // priors
  alpha ~ normal(0,1);
  beta  ~ normal(0,1);
  y_sd  ~ gamma(1,1); 
  
  for(n in 1:n_time)
    mu[n] = exp(alpha + x[n] * beta);
    
  y ~ gamma(y_sd, y_sd ./ mu);
  
}

generated quantities {
  vector[n_time] log_lik;
  
  for (n in 1:n_time)
    log_lik[n] = gamma_lpdf(y[n] | y_sd, (y_sd / exp(alpha + x[n] * beta)) );

}
