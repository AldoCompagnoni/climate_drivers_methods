
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
  int<lower=1> n_time;  // number of observations
  vector[n_time] y;  // response variable
  int<lower=1> n_lag;  // number of population-level effects
  matrix[n_time, n_lag] clim;  // population-level design matrix
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
  real<lower=0> y_sd;  // residual SD
}

transformed parameters {
  vector[n_lag] beta;  // population-level effects
  vector[n_time] yhat;
  // compute actual regression coefficients
  beta = horseshoe(zb, hs_local, hs_global, hs_scale_global * y_sd, hs_scale_slab^2 * hs_c2);
  yhat = alpha + clim * beta;
}
model {
  // priors including all constants
  target += normal_lpdf(zb | 0, 1);
  target += normal_lpdf(hs_local[1] | 0, 1)
    - 36 * log(0.5);
  target += inv_gamma_lpdf(hs_local[2] | 0.5 * hs_df, 0.5 * hs_df);
  target += normal_lpdf(alpha | 0, 2);
  target += normal_lpdf(hs_global[1] | 0, 1)
    - 1 * log(0.5);
  target += inv_gamma_lpdf(hs_global[2] | 0.5 * hs_df_global, 0.5 * hs_df_global);
  target += inv_gamma_lpdf(hs_c2 | 0.5 * hs_df_slab, 0.5 * hs_df_slab);
  target += student_t_lpdf(y_sd | 3, 0, 10)
    - 1 * student_t_lccdf(0 | 3, 0, 10);
  // likelihood including all constants
  y ~ normal( yhat, y_sd);
  
}

generated quantities {
  vector[n_time] log_lik;
  
  for (n in 1:n_time)
  log_lik[n] = normal_lpdf(y[n] | yhat[n], y_sd );

}
