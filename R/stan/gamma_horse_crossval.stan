
functions {

  /* Efficient computation of the horseshoe prior
   * Args:
   *   zb: standardized population-level coefficients
   *   global: global horseshoe parameters
   *   local: local horseshoe parameters
   *   scale_global: global scale of the horseshoe prior
   *   c2: positive real number for regularization
   * Returns:
   *   population-level coefficients following the horseshoe prior
   */
  vector horseshoe(vector zb, vector[] local, real[] global,
                   real scale_global, real c2) {
    int n_lag = rows(zb);
    vector[n_lag] lambda = local[1] .* sqrt(local[2]);
    vector[n_lag] lambda2 = square(lambda);
    real tau = global[1] * sqrt(global[2]) * scale_global;
    vector[n_lag] lambda_tilde = sqrt(c2 * lambda2 ./ (c2 + tau^2 * lambda2));
    return zb .* lambda_tilde * tau;
  }
}
data {
  int<lower=1> n_train;  // number of observations
  int<lower=1> n_test;  // number of observations
  int n_lag;
  vector[n_train] y_train;
  vector[n_test]  y_test;
  matrix[n_train, n_lag] clim_train;
  matrix[n_test,  n_lag] clim_test;
  // data for the horseshoe prior
  real<lower=0> hs_df;  // local degrees of freedom
  real<lower=0> hs_df_global;  // global degrees of freedom
  real<lower=0> hs_df_slab;  // slab degrees of freedom
  real<lower=0> hs_scale_global;  // global prior scale
  real<lower=0> hs_scale_slab;  // slab prior scale
}

parameters {
  // local parameters for horseshoe prior
  vector[n_lag] zb;
  vector<lower=0>[n_lag] hs_local[2];
  real alpha;  // temporary intercept for centered predictors
  // horseshoe shrinkage parameters
  real<lower=0> hs_global[2];  // global shrinkage parameters
  real<lower=0> hs_c2;  // slab regularization parameter
  real<lower=0> y_sd;  // shape parameter
}

transformed parameters {
  vector[n_lag] beta;  // population-level effects
  vector[n_train] yhat;
  vector[n_train] mu_aux;
  vector[n_train] mu;
  
  // compute actual regression coefficients
  beta   = horseshoe(zb, hs_local, hs_global, hs_scale_global, hs_scale_slab^2 * hs_c2);
  yhat   = exp(alpha + clim_train * beta); // store means
  mu_aux = alpha + clim_train * beta;      // auxiliary variable
  
  // initialize linear predictor term
  for (n in 1:n_train) {
    // apply the inverse link function
    mu[n] = y_sd * exp(-mu_aux[n]);
  }
  
}
model {
  
  // priors including all constants
  target += normal_lpdf(zb | 0, 1);
  target += normal_lpdf(hs_local[1] | 0, 1)
    - 36 * log(0.5);
  target += inv_gamma_lpdf(hs_local[2] | 0.5 * hs_df, 0.5 * hs_df);
  target += normal_lpdf(alpha | 0, 0.5);
  target += normal_lpdf(hs_global[1] | 0, 1)
    - 1 * log(0.5);
  target += inv_gamma_lpdf(hs_global[2] | 0.5 * hs_df_global, 0.5 * hs_df_global);
  target += inv_gamma_lpdf(hs_c2 | 0.5 * hs_df_slab, 0.5 * hs_df_slab);
  target += gamma_lpdf(y_sd | 0.01, 0.01);
  // likelihood including all constants
  target += gamma_lpdf(y_train | y_sd, mu);
}

generated quantities {
  vector[n_train] log_lik;
  vector[n_lag] pred_x;
  vector[n_test] pred_y;
  vector[n_test] log_lik_test;
  vector[n_test] mu_aux_test;
  vector[n_test] mu_test;
  
  for(n in 1:n_train)
    log_lik[n] = gamma_lpdf(y_train[n] | y_sd, mu[n] );
  
  // out of sample prediction
  for(n in 1:n_test){
    pred_y[n]       = exp(alpha + sum(clim_test[n,]' .* beta) ); // store means
    mu_aux_test[n]  = alpha + sum(clim_test[n,]' .* beta);      // auxiliary variable
    // apply the inverse link function
    mu_test[n]      = y_sd * exp(-mu_aux_test[n]);
    log_lik_test[n] = gamma_lpdf(y_test[n] | y_sd, mu_test[n]);
  }

}
