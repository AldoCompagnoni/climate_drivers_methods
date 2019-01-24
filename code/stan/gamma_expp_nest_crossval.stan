
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
  // place holder  
  vector[n_train] mu;    // transformed linear predictor for mean of beta distribution
  
  // likelihood
  for(n in 1:n_train)
    mu[n] = exp(alpha + x[n] * beta);
    
  y_train ~ gamma(y_sd, y_sd ./ mu);
}

generated quantities {
  vector[n_train]  log_lik;
  vector[n_test]   log_lik_test;
  matrix[K,n_test] pred_x_m;
  vector[n_test]   pred_y;
  vector[n_test]   pred_x;

  for(n in 1:n_train)
    log_lik[n] = gamma_lpdf(y_train[n] | y_sd, (y_sd / exp(alpha + x[n] * beta)) );
  
  // out of sample prediction
  for(n in 1:n_test){
    pred_x_m[1,n] = sum(sens_m  .* clim1_test[,n]); 
    pred_x_m[2,n] = sum(sens_m  .* clim2_test[,n]);
    pred_x_m[3,n] = sum(sens_m  .* clim3_test[,n]);
    pred_x[n]       = sum(theta_y .* pred_x_m[,n]);
    pred_y[n]       = exp(alpha + pred_x[n] * beta);
    log_lik_test[n] = gamma_lpdf(y_test[n] | y_sd, (y_sd / exp(alpha + pred_x[n] * beta)) );
  }
}
