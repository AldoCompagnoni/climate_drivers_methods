
data {
  int n_time;
  int K;
  vector[n_time] y;
  vector[n_time] clim_means;
}

parameters {
  real alpha;
  real beta;
  real<lower=0> y_sd;
}

transformed parameters {
  // place holder  
  vector[n_time] yhat; // transformed linear predictor for mean of beta distribution
  vector<lower=0>[n_time] B;          // parameter for gamma distn
  
  yhat  = exp(alpha + clim_means * beta);
  B     = y_sd ./ yhat;
  
}

model {

  // parameters of data model
  alpha ~ normal(0,0.5);
  beta  ~ normal(0,1);
  y_sd  ~ gamma(0.01,0.01); 

  y ~ gamma(y_sd, B);
  
}

generated quantities {
  vector[n_time] log_lik;
  
  for (n in 1:n_time)
    log_lik[n] = gamma_lpdf(y[n] | y_sd, B[n]);
}
