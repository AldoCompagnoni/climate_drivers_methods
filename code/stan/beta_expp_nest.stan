
functions {
  real dexppow(real x, real mu, real sigma, real beta) {
    return((beta / (2 * sigma * tgamma(1.0/beta)) ) * exp(-(fabs(x - mu)/sigma)^beta));
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
  real expp_beta;
}

parameters {
  real<lower=0,upper=M> sens_mu;
  real<lower=0,upper=M*2> sens_sd;
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
  
  for(i in 1:M)
    sens_m[i] = dexppow(i, sens_mu, sens_sd, expp_beta);
  
  sens_m = sens_m / sum(sens_m);
  
  for(i in 1:n_time) {
    x_m[1,i] = sum(sens_m .* clim1[,i]); 
    x_m[2,i] = sum(sens_m .* clim2[,i]);
    x_m[3,i] = sum(sens_m .* clim3[,i]);
  }

  for(i in 1:n_time)
    x[i] = sum(theta_y .* x_m[,i]);
    
  for(n in 1:n_time){
    mu[n]  = inv_logit(alpha + x[n] * beta);
    A[n]   = mu[n] * y_sd;
    B[n]   = (1.0 - mu[n]) * y_sd;
  }
}

model {
  // likelihood
  y ~ beta(A, B);
}

generated quantities {
  vector[n_time] log_lik;
  
  for (n in 1:n_time)
    log_lik[n] = beta_lpdf(y[n] | A[n], B[n]);
}
