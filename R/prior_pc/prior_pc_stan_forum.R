library(rstan)
library(MASS)
library(tidyverse)

# STAN models 
suggested_ppc <- '
  
  data {
    int N;
    vector[N] x;
  }
  
  generated quantities {
    real alpha = normal_rng(0,1);
    real beta  = normal_rng(0,1);
    real y_sd  = gamma_rng(0.01,0.01); 
    
    real yhat[N];
    real y_sim[N];
    
    for(n in 1:N)
      yhat[n] = alpha + beta * x[n];
    
    for(n in 1:N)
      y_sim[n] = normal_rng(yhat[n], y_sd);
    
  }

'

brms_ppc      <- '
  
  data {
    int N;
    vector[N] x;
  }
  
  parameters {
    real alpha;
    real beta;
    real<lower=0> y_sd;
  }
  
  transformed parameters{
    vector[N] yhat;
    yhat = alpha + beta * x;
  }
  
  model {
    alpha ~ normal(0,1);
    beta  ~ normal(0,1);
    y_sd  ~ gamma(0.01,0.01); 
  }
  
  generated quantities {
    real y_sim[N];
    for(n in 1:N)
      y_sim[n] = normal_rng(yhat[n], y_sd);
  }

'

# Create stanmodel objects
suggested_ppc_mod  <- stan_model( model_code = suggested_ppc  )
brms_ppc_mod       <- stan_model( model_code = brms_ppc )

# create fictitious predictor data
set.seed( 2013 + 702)
dat <- list( x = rnorm( 10 ),
             N = 10 )

# fit stan models 
suggested_ppc_fit <- sampling( object = suggested_ppc_mod, 
                               data   = dat,
                               iter   = 4000,
                               warmup = 1000,
                               thin   = 2,
                               chains = 3,
                               algorithm = 'Fixed_param'
                              )

brms_ppc_fit      <- sampling( object = brms_ppc_mod, 
                               data   = dat,
                               iter   = 4000,
                               warmup = 1000,
                               thin   = 2,
                               chains = 3 
                            )

# The variance of output suggested_ppc always > brms_ppc
suggested_ppc_fit %>% 
  rstan::extract() %>% 
  .$y_sim %>% 
  as.numeric %>% 
  sd(na.rm=T) # Did not understand why there are NAs!

brms_ppc_fit %>% 
  rstan::extract() %>% 
  .$y_sim %>% 
  as.numeric %>% 
  sd()
