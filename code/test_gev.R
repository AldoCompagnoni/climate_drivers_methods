# script to test that my code for GEV is correct
library(tidyverse)
library(mgcv)
library(testthat)
library(rstan)
library(evd)
library(dplyr)

# set rstan options
rstan_options( auto_write = TRUE )
options( mc.cores = parallel::detectCores() )

# organize data into list to pass to stan
dat_stan <- list(
  N  = 100,
  y  = rgev(100,5,1,0)
)

# simulation parameters
sim_pars <- list(
  warmup = 1000, 
  iter = 4000, 
  thin = 2, 
  chains = 3
)

fit_gev <- stan(
  file = paste0("code/stan/gev.stan"),
  data = dat_stan,
  pars = c('mu', 'sigma', 'xi'), #, "shape", 'log_lik'
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  init= 2,
  chains = sim_pars$chains#,
  # control = list(adapt_delta = 0.999, stepsize = 0.001, max_treedepth = 20)
)

library(shinystan)
launch_shinystan(fit_gev)


fit_gev_s <- stan(
  file = paste0("code/stan/gev_straight.stan"),
  data = dat_stan,
  pars = c('loc', 'scale', 'shape'), #, "shape", 'log_lik'
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  init = 2,
  chains = sim_pars$chains#,
  # control = list(adapt_delta = 0.999, stepsize = 0.001, max_treedepth = 20)
)

launch_shinystan(fit_gev_s)

fit_gev_sc <- stan(
  file = paste0("code/stan/gev_scraped.stan"),
  data = dat_stan,
  pars = c('sigma', 'xi', 'mu'), #, "shape", 'log_lik'
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  init_r = 2,
  chains = sim_pars$chains#,
  # control = list(adapt_delta = 0.999, stepsize = 0.001, max_treedepth = 20)
)

fit_gev_sc

launch_shinystan(fit_gev_sc)
