library(rstan)
library(tidyverse)
library(moments)

stud_t <- ' 

data {
  real df;
  real scale;
}

generated quantities {
  
  real y_sim = student_t_rng(df,0,scale);
  
}

'

stud_mod    <- stan_model( model_code = stud_t )

stud_fit    <- sampling( object = stud_mod, 
                         data   = list( df = 25, scale = 1),
                         iter   = 4000,
                         warmup = 1000,
                         thin   = 2,
                         chains = 1,
                         algorithm = 'Fixed_param'
                      )

stud_fit %>% rstan::extract() %>% .$y_sim %>% hist
stud_fit %>% rstan::extract() %>% .$y_sim %>% sd
stud_fit %>% rstan::extract() %>% .$y_sim %>% kurtosis
