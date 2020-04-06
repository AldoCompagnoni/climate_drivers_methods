
functions {
  real dexppow(real x, real mu, real sigma, real beta) {
    return((beta / (2 * sigma * tgamma(1.0/beta)) ) * exp(-(fabs(x - mu)/sigma)^beta));
  }
}

data {
  int n_train;
  int n_test;
  int K;
  int M;
 
  vector[n_train] y_train;
  vector[n_test]  y_test;

  matrix[M, n_train] clim1_train;
  matrix[M, n_train] clim2_train;
  matrix[M, n_train] clim3_train;
  matrix[M, n_test] clim1_test;
  matrix[M, n_test] clim2_test;
  matrix[M, n_test] clim3_test;
  
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
  
  vector[n_train] x;
  vector[M] sens_m;
  matrix[K,n_train] x_m;
  
  for(i in 1:M)
    sens_m[i] = dexppow(i, sens_mu, sens_sd, expp_beta);
  
  sens_m = sens_m / sum(sens_m);
  
  for(i in 1:n_train) {
    x_m[1,i] = sum(sens_m .* clim1_train[,i]); 
    x_m[2,i] = sum(sens_m .* clim2_train[,i]);
    x_m[3,i] = sum(sens_m .* clim3_train[,i]);
  }

  for(i in 1:n_train)
    x[i] = sum(theta_y .* x_m[,i]);
  
}

model {
  
  // priors
  alpha ~ normal(0,1);
  beta  ~ normal(0,1);
  y_sd  ~ gamma(1,1);
  sens_sd ~ normal(0.5, 12);
  sens_mu ~ normal(6.5, 12); 
  theta_y ~ dirichlet(rep_vector(1.0, K));
  
  // model
  y_train ~ normal(alpha + beta * x, y_sd);
}

generated quantities {
   
  vector[n_train]  log_lik;
  vector[n_test]   pred_y;
  vector[n_test]   log_lik_test;
  matrix[K,n_test] pred_x_m;
  real             pred_x;
  
  for(n in 1:n_train)
    log_lik[n] = normal_lpdf(y_train[n] | alpha + beta * x[n], y_sd);
  
  // out of sample prediction
  for(n in 1:n_test){
    pred_x_m[1,n] = sum(sens_m  .* clim1_test[,n]); 
    pred_x_m[2,n] = sum(sens_m  .* clim2_test[,n]);
    pred_x_m[3,n] = sum(sens_m  .* clim3_test[,n]);
    pred_x        = sum(theta_y .* pred_x_m[,n]);

    pred_y[n]       = alpha + beta * pred_x;
    log_lik_test[n] = normal_lpdf(y_test[n] | alpha + beta * pred_x, y_sd);
  }
}
