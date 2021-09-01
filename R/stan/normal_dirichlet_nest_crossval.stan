
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

model {
  
  // priors
  alpha   ~ normal(0,1);
  beta    ~ normal(0,1);
  y_sd    ~ gamma(0.01,0.01); 
  theta_y ~ dirichlet(rep_vector(1.0, K));
  theta_m ~ dirichlet(rep_vector(1.0, M));

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
    
    pred_x_m[1,n] = sum(theta_m .* clim1_test[,n]); 
    pred_x_m[2,n] = sum(theta_m .* clim2_test[,n]);
    pred_x_m[3,n] = sum(theta_m .* clim3_test[,n]);
    pred_x        = sum(theta_y .* pred_x_m[,n]);

    pred_y[n]       = alpha + beta * pred_x;
    log_lik_test[n] = normal_lpdf(y_test[n] | alpha + beta * pred_x, y_sd);
  }

}

