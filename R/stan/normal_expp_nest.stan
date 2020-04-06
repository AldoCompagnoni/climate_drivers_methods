
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
  real<lower=1> sens_sd;
  simplex[K] theta_y;
  real alpha;
  real beta;
  real<lower=0> y_sd;
}

transformed parameters {
  vector[n_time] x;
  vector[M] sens_m;
  matrix[K,n_time] x_m;
  vector[n_time] yhat;
  
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

  yhat = alpha + beta * x;
  
}

model {
  
  // priors
  alpha ~ normal(0,1);
  beta  ~ normal(0,1);
  y_sd  ~ gamma(1,1);
  sens_sd ~ normal(0.5, 12);
  sens_sd ~ normal(6.5, 12);
  theta_y ~ dirichlet(rep_vector(1.0, K));
  
  // model
  y ~ normal(yhat, y_sd);
}

generated quantities {
  vector[n_time] log_lik;

  for (n in 1:n_time)
    log_lik[n] = normal_lpdf(y[n] | yhat[n], y_sd);
}

