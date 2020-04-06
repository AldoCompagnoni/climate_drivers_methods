
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
  real<lower=1> sens_sd;
  real alpha;
  real beta;
  real<lower=0> y_sd;
}

transformed parameters{
  
  // transformed parameters for moving window
  vector[n_train] x;
  row_vector[n_lag] sens;
  vector[n_train] yhat;
  
  for(i in 1:n_lag) 
    sens[i] = dexppow(i, sens_mu, sens_sd, expp_beta);
  
  sens = sens / sum(sens);
  
  for(n in 1:n_train)
    x[n] = sum(sens .* row(clim_train, n));
    
  yhat = exp(alpha + x * beta);
    
}

model {
  
  // hyper-parameters to weight climate effects
  sens_sd ~ normal(0.5, 12);
  sens_mu ~ normal(6.5, 12);
  
  // parameters of data model
  alpha ~ normal(0,1);
  beta  ~ normal(0,1);
  y_sd  ~ gamma(1,1); 
  
  y_train ~ gamma(y_sd, y_sd ./ yhat);
  
}

generated quantities {
  
  vector[n_train] log_lik;
  vector[n_test]  log_lik_test;
  vector[n_test]  pred_y; // transf. lin. pred. for mean of beta distrib.

  for (n in 1:n_train)
    log_lik[n] = gamma_lpdf(y_train[n] | y_sd, (y_sd / yhat[n]) );

  for(n in 1:n_test){
    pred_y[n]       = exp(alpha + sum(sens .* row(clim_test,n)) * beta);
    log_lik_test[n] = gamma_lpdf(y_test[n] | y_sd, (y_sd / pred_y[n]) );
  }

}
