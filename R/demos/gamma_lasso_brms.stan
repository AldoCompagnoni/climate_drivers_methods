
// generated with brms 2.11.0
functions {
}
data {
  int<lower=1> N;  // number of observations
  vector[N] Y;  // response variable
  int<lower=1> K;  // number of population-level effects
  matrix[N, K] X;  // population-level design matrix
  // data for the lasso prior
  real<lower=0> lasso_df;  // prior degrees of freedom
  real<lower=0> lasso_scale;  // prior scale
  int prior_only;  // should the likelihood be ignored?
}
transformed data {
  int Kc = K - 1;
  matrix[N, Kc] Xc;  // centered version of X without an intercept
  vector[Kc] means_X;  // column means of X before centering
  for (i in 2:K) {
    means_X[i - 1] = mean(X[, i]);
    Xc[, i - 1] = X[, i] - means_X[i - 1];
  }
}
parameters {
  vector[Kc] b;  // population-level effects
  // temporary intercept for centered predictors
  real Intercept;
  // lasso shrinkage parameter
  real<lower=0> lasso_inv_lambda;
  real<lower=0> shape;  // shape parameter
}
transformed parameters {
}
model {
  // initialize linear predictor term
  vector[N] mu = Intercept + Xc * b;
  for (n in 1:N) {
    // apply the inverse link function
    mu[n] = shape * exp(-(mu[n]));
  }
  // priors including all constants
  target += double_exponential_lpdf(b | 0, lasso_scale * lasso_inv_lambda);
  target += normal_lpdf(Intercept | 0, 2);
  target += chi_square_lpdf(lasso_inv_lambda | lasso_df);
  target += gamma_lpdf(shape | 0.01, 0.01);
  // likelihood including all constants
  if (!prior_only) {
    target += gamma_lpdf(Y | shape, mu);
  }
}
generated quantities {
  // actual population-level intercept
  real b_Intercept = Intercept - dot_product(means_X, b);
}
