
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
  real<lower=0,upper=M*2> sens_sd;
  simplex[K] theta_y; 
  real alpha;
  real beta;
  real<lower=0> y_sd;
}

transformed parameters {
  
  real<lower=0,upper=1> mu[n_train];    // transformed linear predictor for mean of beta distribution
  real<lower=0> A[n_train];             // parameter for beta distn
  real<lower=0> B[n_train];             // parameter for beta distn

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
  
  for(n in 1:n_train){
    mu[n]  = inv_logit(alpha + x[n] * beta);
    A[n]   = mu[n] * y_sd;
    B[n]   = (1.0 - mu[n]) * y_sd;
  }
}

model {
  // model
  y_train ~ beta(A, B);
}

generated quantities {
   
  vector[n_train]  log_lik;
  vector[n_test]   log_lik_test;
  vector[n_test]   pred_y;
  matrix[K,n_test] pred_x_m;
  real             pred_x;
  real<lower=0>    A_test[n_test]; // parameter for beta distn
  real<lower=0>    B_test[n_test]; // parameter for beta distn

  for(n in 1:n_train)
    log_lik[n] = beta_lpdf(y_train[n] | A[n], B[n]);
  
  // out of sample prediction
  for(n in 1:n_test){
    pred_x_m[1,n] = sum(sens_m  .* clim1_test[,n]); 
    pred_x_m[2,n] = sum(sens_m  .* clim2_test[,n]);
    pred_x_m[3,n] = sum(sens_m  .* clim3_test[,n]);
    pred_x          = sum(theta_y .* pred_x_m[,n]);
    pred_y[n]       = inv_logit(alpha + pred_x * beta);
    A_test[n]       = pred_y[n] * y_sd;
    B_test[n]       = (1.0 - pred_y[n]) * y_sd;
    log_lik_test[n] = beta_lpdf(y_test[n] | A_test[n], B_test[n]);
  }
}
