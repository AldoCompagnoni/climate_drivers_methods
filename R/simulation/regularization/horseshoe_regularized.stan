data {
  int < lower =0> n; // number of observations
  int < lower =0> d; // number of predictors
  vector [n] y; // outputs
  matrix [n,d] x; // inputs
  real < lower =0> scale_icept ; // prior std for the intercept
  real < lower =0> scale_global ; // scale for the half -t prior for tau
  real < lower =1> nu_global ; // degrees of freedom for the half -t prior
  // for tau
  real < lower =1> nu_local ; // degrees of freedom for the half -t priors
  // for lambdas
  real < lower =0> slab_scale ; // slab scale for the regularized horseshoe
  real < lower =0> slab_df ; // slab degrees of freedom for the regularized
  // horseshoe
}
parameters {
  real logsigma ;
  real Intercept ;
  vector [d] z;
  real < lower =0> tau; // global shrinkage parameter
  vector < lower =0 >[d] lambda ; // local shrinkage parameter
  real < lower =0> caux ;
}

transformed parameters {
  real < lower =0> sigma ; // noise std
  vector < lower =0 >[d] lambda_tilde ; // 'truncated ' local shrinkage parameter
  real < lower =0> c; // slab scale
  vector [d] b ; // regression coefficients
  vector [n] f; // latent function values
  sigma = exp ( logsigma );
  c = slab_scale * sqrt ( caux );
  lambda_tilde = sqrt ( c^2 * square ( lambda ) ./ (c^2 + tau ^2* square ( lambda )) );
  b = z .* lambda_tilde *tau;
  f = Intercept + x * b ;
}

model {
  // half -t priors for lambdas and tau , and inverse - gamma for c^2
  z ~ normal (0, 1);
  lambda ~ student_t ( nu_local , 0, 1); // Note mean 0, sigma 1
  tau ~ student_t ( nu_global , 0, scale_global * sigma );
  caux ~ inv_gamma (0.5* slab_df , 0.5* slab_df );
  
  Intercept ~ normal (0, scale_icept );
  
  y ~ normal (f, sigma);
}
