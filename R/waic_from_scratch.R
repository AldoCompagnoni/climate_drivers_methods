# CHECK: example of how to calculate WAIC from skratch
setwd("C:/CODE/climate_drivers_methods")
library(rstan)
library(loo)
library(dplyr)
library(testthat)

# organize data into list to pass to stan
dat_stan <- list(
  n_time  = 30,
  y       = rnorm(30, 5, 1)
)

# simulation parameters
sim_pars <- list(
  warmup = 1000, 
  iter = 4000, 
  thin = 2, 
  chains = 3
)

# NULL model (model of the mean)
fit_ctrl1 <- stan(
  file = paste0("code/stan/normal_null.stan"),
  data = dat_stan,
  pars = c('alpha', 'y_sd', 'log_lik'),
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains
  #control = list(adapt_delta = 0.999, stepsize = 0.001, max_treedepth = 20)
)

# lppd
lppd <- fit_ctrl1 %>% 
  extract() %>% 
  as.data.frame %>% 
  select( grep('log_lik',names(.)) ) %>% 
  apply(2,exp) %>% 
  apply(2,mean) %>% 
  log %>% sum

# variance of the lppd (effective n. of params in WAIC formula)
lppd_var <- fit_ctrl1 %>% 
  extract() %>% 
  as.data.frame %>% 
  select( grep('log_lik',names(.)) ) %>% 
  apply(2,var) %>% 
  sum

# here is my estimate
waic_skratch <- lppd - lppd_var

# waic through loo versus waic from skratch
loo_waic <- extract_log_lik(fit_ctrl1) %>% 
  waic %>% 
  # get estimates table
  .$estimates %>% 
  # get estimates
  .[,'Estimate'] %>% 
  # get WAIC
  .[1] %>% 
  # remove names (otherwise error from expect_equal)
  setNames(NULL) %>% 
  # this need be equal to 
  expect_equal( waic_skratch )

