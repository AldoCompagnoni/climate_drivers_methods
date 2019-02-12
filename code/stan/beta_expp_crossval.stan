
functions {
  real dexppow(real x, real mu, real sigma, real beta) {
    return((beta / (2 * sigma * tgamma(1.0/beta)) ) * exp(-(fabs(x - mu)/sigma)^beta));
  }
}

data {
  int n_train;
  int n_test;
  int n_lag;
  vector[n_train] y_train;
  vector[n_test]  y_test;
  matrix[n_train, n_lag] clim_train;
  matrix[n_test,  n_lag] clim_test;
  real expp_beta;
}

parameters {
  real<lower=0,upper=n_lag> sens_mu;
  real<lower=0,upper=n_lag*2> sens_sd;
  real alpha;
  real beta;
  real<lower=0> y_sd;
}

transformed parameters {
  vector[n_train] x;
  row_vector[n_lag] sens;
  
  real<lower=0,upper=1> mu[n_train];    // transformed linear predictor for mean of beta distribution
  real<lower=0> A[n_train];             // parameter for beta distn
  real<lower=0> B[n_train];             // parameter for beta distn

  for(i in 1:n_lag)
    sens[i] = dexppow(i, sens_mu, sens_sd, expp_beta);
  
  sens = sens / sum(sens);
  
  for(i in 1:n_train)
    x[i] = sum(sens .* row(clim_train, i));
    
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
  
  vector[n_train] log_lik;
  vector[n_test]  log_lik_test;
  vector[n_test]  pred_y;
  real            pred_x;
  real<lower=0>   A_test[n_test];               // parameter for beta distn
  real<lower=0>   B_test[n_test];               // parameter for beta distn

  for (n in 1:n_train)
    log_lik[n] = beta_lpdf(y_train[n] | A[n], B[n]);
  
  // out of sample prediction
  for(n in 1:n_test){
    pred_x          = sum(sens .* row(clim_test,n));
    pred_y[n]       = inv_logit(alpha + pred_x * beta);
    A_test[n]       = pred_y[n] * y_sd;
    B_test[n]       = (1.0 - pred_y[n]) * y_sd;
    log_lik_test[n] = beta_lpdf(y_test[n] | A_test[n], B_test[n]);
  }

}