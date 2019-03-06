
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
}

parameters {
  simplex[K] theta_y;
  simplex[M] theta_m;
  real alpha;
  real beta;
  real<lower=0> y_sd;
}

transformed parameters {

  matrix[K,n_train] x_m;
  vector[n_train] x;
  
  for(i in 1:n_train){
    x_m[1,i] = sum(theta_m .* clim1_train[,i]); 
    x_m[2,i] = sum(theta_m .* clim2_train[,i]);
    x_m[3,i] = sum(theta_m .* clim3_train[,i]);
  }
  
  for(i in 1:n_train)
    x[i] = sum(theta_y .* x_m[,i]);
}

model{
  // place holder  
  vector[n_train] mu; // transf. lin. pred. for mean
  
  // hyper-parameters to weight climate effects
  theta_m ~ dirichlet(rep_vector(1.0, M));
  theta_y ~ dirichlet(rep_vector(1.0, K));

  // priors
  alpha ~ normal(0,1);
  beta  ~ normal(0,1);
  y_sd  ~ gamma(1,1); 

  // likelihood
  for(n in 1:n_train)
    mu[n] = exp(alpha + x[n] * beta);
    
  y_train ~ gamma(y_sd, y_sd ./ mu);
}

generated quantities{
  vector[n_train]  log_lik;
  vector[n_test]   log_lik_test;
  vector[n_test]   pred_y;
  vector[n_test]   pred_x;
  matrix[K,n_test] pred_x_m;

  for(n in 1:n_train)
    log_lik[n] = gamma_lpdf(y_train[n] | y_sd, (y_sd / exp(alpha + x[n] * beta)) );
  
  // out of sample prediction
  for(n in 1:n_test){
    pred_x_m[1,n]   = sum(theta_m .* clim1_test[,n]); 
    pred_x_m[2,n]   = sum(theta_m .* clim2_test[,n]);
    pred_x_m[3,n]   = sum(theta_m .* clim3_test[,n]);
    pred_x[n]       = sum(theta_y .* pred_x_m[,n]);
    pred_y[n]       = exp(alpha + pred_x[n] * beta);
    log_lik_test[n] = gamma_lpdf(y_test[n] | y_sd, (y_sd / exp(alpha + pred_x[n] * beta)) );
  }
}
