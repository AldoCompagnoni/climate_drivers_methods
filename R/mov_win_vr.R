rm(list=ls())
source("R/format_data.R")
library(tidyverse)
library(mgcv)
library(testthat)
library(rstan)
library(loo)

# set rstan options
rstan_options( auto_write = TRUE )
options( mc.cores = parallel::detectCores() )

# climate predictor, response, months back, max. number of knots
response  <- "grow"
clim_var  <- "airt"
m_back    <- 36    
st_dev    <- FALSE

# read data -----------------------------------------------------------------------------------------
lam       <- read.csv("data/all_demog_updt.csv", stringsAsFactors = F)
# m_info    <- read.csv("C:/cloud/Dropbox/sApropos/MatrixEndMonth_information.csv", stringsAsFactors = F)
clim      <- data.table::fread(paste0('data/',clim_var,"_chelsa_prism_hays_2014.csv"),  stringsAsFactors = F)
spp       <- lam$SpeciesAuthor %>% unique

# format data --------------------------------------------------------------------------------------

# set up model "family" based on response
if( response == "surv" | response == "grow" )             family = "beta" 
if( response == "fec" )                                   family = "gamma"
if( grepl("PreRep", response) | grepl("Rep", response) )  family = "beta"
if( response == "rho" | response == "react_fsa" )         family = "gamma"
if( response == "log_lambda" )                             family = "normal"

expp_beta     <- 20

# set species (I pick Sphaeraclea_coccinea)
ii            <- 2
spp_name      <- spp[ii]

# lambda data
spp_resp      <- format_species(spp_name, lam, response)

# climate data
clim_separate <- clim_list(spp_name, clim, spp_resp)
#clim_detrnded <- lapply(clim_separate, clim_detrend, clim_var)
clim_detrnded <- lapply(clim_separate, clim_detrend, clim_var)
clim_mats     <- Map(clim_long, clim_detrnded, spp_resp, m_back)

# model data
mod_data          <- lambda_plus_clim(spp_resp, clim_mats, response)
mod_data$climate  <- mod_data$climate #/ diff(range(mod_data$climate))

# throw error if not enough data
if( nrow(mod_data$resp) < 6 ) stop( paste0("not enough temporal replication for '", 
                                              spp_name, "' and response variable '" , response, "'") )


# Transform response variables (if needed) ------------------------------------------------------------------

# replace 0 and 1s to allow fitting beta models
if( response == "surv" | response == "grow" | grepl("PreRep", response) | grepl("Rep", response) ){
  
  raw_x <- mod_data$resp[,response]
  raw_x <- replace( raw_x, raw_x == 1, 0.99999 )
  raw_x <- replace( raw_x, raw_x == 9, 0.00001 )
  
  # replace numbers
  mod_data$resp[,response] <- new_x

}

# avoid absolute zeros
if( response == "fec" ){
  # transform from [0, infinity) to (0, infinity) 
  # I add quantity 2 orders of mag. lower than lowest obs value.
  mod_data$resp[,response] <- mod_data$resp[,response] + 1.54e-12 
} 

if( response == "rho" | response == "react_fsa" ){
  # bound responses to (0,infinity) instead of [1, infinity) 
  mod_data$resp[,response] <- mod_data$resp[,response] - 0.99999
} 


# Fit models ----------------------------------------------------------------------------------------

# organize data into list to pass to stan
dat_stan <- list(
  n_time  = nrow(mod_data$climate),
  n_lag   = ncol(mod_data$climate),
  y       = mod_data$resp[,response],
  clim    = mod_data$climate,
  clim_means = rowMeans(mod_data$climate),
  clim_yr = list( rowMeans(mod_data$climate[, 1:12]),
                  rowMeans(mod_data$climate[,13:24]),
                  rowMeans(mod_data$climate[,25:36]) ) %>% do.call(rbind,.),
  M       = 12,    # number of months in a year
  K       = ncol(mod_data$climate) / 12, # number of years
  S       = mod_data$resp$population %>% unique %>% length,
  site_i  = mod_data$resp$population %>% as.factor %>% as.numeric,
  expp_beta = expp_beta,
  # parameters for horseshoe models
  hs_df           = 1,   # variance of 
  hs_df_global    = 1,   # 
  hs_df_slab      = 25,   # slab degrees of freedom
  hs_scale_global = (4 / (36-4)) / sqrt(nrow(mod_data$climate)), # global prior scale
  hs_scale_slab   = 2    # slab prior scale
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
  file = paste0("R/stan/",family,"_null.stan"),
  data = dat_stan,
  pars = c('alpha', 'y_sd', 
           'yhat',  'log_lik'),
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains
  #control = list(adapt_delta = 0.999, stepsize = 0.001, max_treedepth = 20)
)

# year t
dat_stan$clim_means  <- rowMeans(mod_data$climate[,1:12 ])
fit_yr1 <- stan(
  file = paste0("R/stan/",family,"_yr.stan"),
  data = dat_stan,
  pars = c('alpha', 'beta', 'y_sd', 
           'yhat','log_lik'),
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains,
  control = list(adapt_delta = 0.8)
)

# year t-1
dat_stan$clim_means  <- rowMeans(mod_data$climate[,13:24])
fit_yr2 <- stan(
  file = paste0("R/stan/",family,"_yr.stan"),
  data = dat_stan,
  pars = c('alpha', 'beta', 'y_sd', 
           'yhat','log_lik'),
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains,
  control = list(adapt_delta = 0.8)
)

# year t-2
dat_stan$clim_means  <- rowMeans(mod_data$climate[,25:36])
fit_yr3 <- stan(
  file = paste0("R/stan/",family,"_yr.stan"),
  data = dat_stan,
  pars = c('alpha', 'beta', 'y_sd', 
           'yhat','log_lik'),
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains,
  control = list(adapt_delta = 0.8)
)
dat_stan$clim_means  <- rowMeans(mod_data$climate)


# gaussian moving window in year t
dat_stan$clim   <- mod_data$climate[,1:12]
dat_stan$n_lag  <- 12
fit_gaus1 <- stan(
  file = paste0("R/stan/",family,"_gaus.stan"),
  data = dat_stan,
  pars = c('sens_mu', 'sens_sd', 'alpha', 'beta', 'y_sd', 
           'yhat','log_lik'), #
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains,
  control = list(adapt_delta = 0.8)
)

# gaussian moving window in year t-1
dat_stan$clim   <- mod_data$climate[,13:24]
fit_gaus2 <- stan(
  file = paste0("R/stan/",family,"_gaus.stan"),
  data = dat_stan,
  pars = c('sens_mu', 'sens_sd', 'alpha', 'beta', 'y_sd', 
           'yhat','log_lik'), #
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains,
  control = list(adapt_delta = 0.8)
)

# gaussian moving window in year t-2
dat_stan$clim   <- mod_data$climate[,25:36]
fit_gaus3 <- stan(
  file = paste0("R/stan/",family,"_gaus.stan"),
  data = dat_stan,
  pars = c('sens_mu', 'sens_sd', 'alpha', 'beta', 'y_sd', 
           'yhat','log_lik'), #
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains,
  control = list(adapt_delta = 0.8)
)

# simplex year t
dat_stan$clim   <- t(mod_data$climate[,1:12])
fit_simpl1 <- stan(
  file = paste0("R/stan/",family,"_dirichlet.stan"),
  data = dat_stan,
  pars = c('theta', 'alpha', 'beta', 'y_sd', 
           'yhat',  'log_lik'),
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains#,
  #control = list(adapt_delta = 0.999, stepsize = 0.001, max_treedepth = 20)
)

# simplex year t-1
dat_stan$clim   <- t(mod_data$climate[,13:24])
fit_simpl2 <- stan(
  file = paste0("R/stan/",family,"_dirichlet.stan"),
  data = dat_stan,
  pars = c('theta', 'alpha', 'beta', 'y_sd', 
           'yhat',  'log_lik'),
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains#,
  #control = list(adapt_delta = 0.999, stepsize = 0.001, max_treedepth = 20)
)

# simplex year t-2
dat_stan$clim   <- t(mod_data$climate[,25:36])
fit_simpl3 <- stan(
  file = paste0("R/stan/",family,"_dirichlet.stan"),
  data = dat_stan,
  pars = c('theta', 'alpha', 'beta', 'y_sd', 
           'yhat',  'log_lik'),
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains#,
  #control = list(adapt_delta = 0.999, stepsize = 0.001, max_treedepth = 20)
)

# Ridge regression at time t
dat_stan$clim   <- mod_data$climate[,1:12]
fit_ridge1 <- stan(
  file = paste0("R/stan/",family,"_horse.stan"),
  data = dat_stan,
  pars = c('alpha', 'beta', 'y_sd', 'yhat', 'log_lik'),
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains,
  control = list(adapt_delta = 0.99, stepsize = 0.001, max_treedepth = 20)
)

# Ridge regression at time t-1
dat_stan$clim   <- mod_data$climate[,13:24]
fit_ridge2 <- stan(
  file = paste0("R/stan/",family,"_horse.stan"),
  data = dat_stan,
  pars = c('alpha', 'beta', 'y_sd', 'yhat', 'log_lik'),
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains,
  control = list(adapt_delta = 0.99, stepsize = 0.001, max_treedepth = 20)
)

# Ridge regression at time t-2
dat_stan$clim   <- mod_data$climate[,25:36]
fit_ridge3 <- stan(
  file = paste0("R/stan/",family,"_horse.stan"),
  data = dat_stan,
  pars = c('alpha', 'beta', 'y_sd', 'yhat', 'log_lik'),
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains,
  control = list(adapt_delta = 0.99, stepsize = 0.001, max_treedepth = 20)
)


# gaussian moving window ALL THREE YEARS
dat_stan$clim   <- mod_data$climate
dat_stan$n_lag  <- 36
fit_gaus_all <- stan(
  file = paste0("R/stan/",family,"_gaus.stan"),
  data = dat_stan,
  pars = c('sens_mu', 'sens_sd', 'alpha', 'beta', 'y_sd', 
           'yhat','log_lik'),
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains,
  control = list(adapt_delta = 0.8)
)

# Ridge ALL YEARS
fit_ridge_all <- stan(
  file = paste0("R/stan/",family,"_horse.stan"),
  data = dat_stan,
  pars = c('alpha', 'beta', 'y_sd', 'yhat', 'log_lik'),
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains,
  control = list(adapt_delta = 0.99, stepsize = 0.001, max_treedepth = 20)
)

# Nested simplex models 
dat_stan$clim         <- t(mod_data$climate)
dat_stan$clim1        <- t(mod_data$climate)[1:12 ,]
dat_stan$clim2        <- t(mod_data$climate)[13:24,]
dat_stan$clim3        <- t(mod_data$climate)[25:36,]
fit_simpl_all <- stan(
  file = paste0("R/stan/",family,"_dirichlet_nest.stan"),
  data = dat_stan,
  pars = c('theta_y', 'theta_m', 'alpha', 'beta', 'y_sd', 
           'yhat',    'log_lik'),
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains#,
  #control = list(adapt_delta = 0.999, stepsize = 0.001, max_treedepth = 20)
)


# parameter values and diagnostics ----------------------------------------------------------------

# list of model fits
mod_fit   <- list( ctrl1    = fit_ctrl1,     
                   yr1      = fit_yr1,     
                   yr2      = fit_yr2,
                   yr3      = fit_yr3,
                   gaus1    = fit_gaus1,
                   gaus2    = fit_gaus2,
                   gaus3    = fit_gaus3,
                   ridge1   = fit_ridge1,
                   ridge2   = fit_ridge2,
                   ridge3   = fit_ridge3,
                   simpl1   = fit_simpl1, 
                   simpl2   = fit_simpl2, 
                   simpl3   = fit_simpl3, 
                   gaus     = fit_gaus_all,
                   ridge    = fit_ridge_all,
                   simpl_n  = fit_simpl_all )
                   

# get central tendencies
pars_diag_extract <- function(x){
  
  # central tendencies
  tmp         <- rstan::extract(x) %>% as.data.frame
  tmp         <- setNames(tmp, gsub("\\.","_",names(tmp)))
  par_means   <- sapply(tmp, function(x) mean(x)) %>%
                    setNames( paste0(names(tmp),"_mean") )
  par_medians <- sapply(tmp, function(x) median(x)) %>%
                    setNames( paste0(names(tmp),"_median") )
  central_tend<- c(par_means, par_medians)
  
  # diagnostics
  diverg      <- do.call(rbind, args = get_sampler_params(x, inc_warmup = F))[,5]
  n_diverg    <- length(which(diverg == 1))
  df_summ     <- as.data.frame(summary(x)$summary)
  rhat_high   <- length(which(df_summ$Rhat > 1.1))
  n_eff       <- df_summ$n_eff / length(diverg)
  n_eff_low   <- length(which(n_eff < 0.1))
  mcse_high   <- length(which(df_summ$se_mean / df_summ$sd > 0.1))
  diagnostics <- c(n_diverg = n_diverg, rhat_high = rhat_high,
                   n_eff_low = n_eff_low, mcse_high = mcse_high)
  out         <- c( central_tend, diagnostics ) %>% t %>% as.data.frame
  
  rm(tmp) ; return(out)
  
}

# extract rhat~n_eff for all parameters
all_diag_extract <- function(x,y){
  
  as.data.frame(summary(x)$summary) %>%
    dplyr::select(n_eff, Rhat, se_mean, sd ) %>%
    tibble::add_column(.before=1, model = y)
  
}

# store posteriors
posterior_extract <- function(model_fit, model_name){
  
  # central tendencies
  tmp     <- rstan::extract(model_fit) %>% as.data.frame
  post_df <- setNames(tmp, gsub("\\.","_",names(tmp)))
  post_df <- tibble::add_column(post_df,
                                model = model_name, .before=1)
  
  rm(tmp) ; return(post_df)
  
}

# calculate central tendencies
pars_diag_l   <- lapply(mod_fit, pars_diag_extract)
mod_pars_diag <- Reduce(function(...) bind_rows(...), pars_diag_l) %>%
                    tibble::add_column(model = names(mod_fit), .before = 1)

# store posteriors
posts_l       <- Map(posterior_extract, mod_fit, names(mod_fit) )
posteriors    <- bind_rows(posts_l)

# extract rhat/neff
diag_df       <- Map( all_diag_extract,
                      mod_fit,
                      names(mod_fit) ) %>% bind_rows


# WAIC model comparison --------------------------------------------------------------------

# wAIC model selection using loo approximation (from library 'loo')
log_liks   <- lapply(mod_fit, extract_log_lik)

# data frame to "name" models
mod_names  <- data.frame( model = c("ctrl1",  
                                    "yr1",     "yr2",   "yr3", 
                                    "gaus1",   "gaus2", "gaus3", 
                                    "simpl1",  "simpl2","simpl3", 
                                    "ridge1",  "ridge2","ridge3", 
                                    "gaus",    "ridge", "simpl_n"),
                          mod_n = c(1:16),
                          stringsAsFactors = F )

# leave-one-out estimates
loo_l      <- lapply(log_liks, loo) %>%
                  setNames( c("loo_ctrl1",   
                              "loo_yr1",     "loo_yr2",   "loo_yr3", 
                              "loo_gaus1",   "loo_gaus2", "loo_gaus3", 
                              "loo_simpl1",  "loo_simpl2","loo_simpl3", 
                              "loo_ridge1",  "loo_ridge2","loo_ridge3", 
                              "loo_gaus",    "loo_ridge", "loo_simpl_n") )
loo_df     <- loo_compare(loo_l$loo_ctrl1,   
                          loo_l$loo_yr1,     loo_l$loo_yr2,    loo_l$loo_yr3, 
                          loo_l$loo_gaus1,   loo_l$loo_gaus2,  loo_l$loo_gaus3, 
                          loo_l$loo_simpl1,  loo_l$loo_simpl2, loo_l$loo_simpl3, 
                          loo_l$loo_ridge1,  loo_l$loo_ridge2, loo_l$loo_ridge3, 
                          loo_l$loo_gaus,    loo_l$loo_ridge,  loo_l$loo_simpl_n
                           ) %>%
                as.data.frame %>%
                tibble::add_column(model = gsub("loo_l\\$loo_","",row.names(.) ), 
                                   .before = 1) %>% 
                rename( se_diff_loo = se_diff ) %>% 
                mutate( model = gsub('model','',model) %>% as.numeric ) %>% 
                rename( mod_n = model ) %>%   
                left_join( mod_names )  

# WAIC estimates
waic_l    <- lapply(log_liks, waic) %>%
                setNames(c("waic_ctrl1",   
                           "waic_yr1",     "waic_yr2",   "waic_yr3", 
                           "waic_gaus1",   "waic_gaus2", "waic_gaus3", 
                           "waic_simpl1",  "waic_simpl2","waic_simpl3", 
                           "waic_ridge1",  "waic_ridge2","waic_ridge3", 
                           "waic_gaus",    "waic_ridge", "waic_simpl_n") )
waic_df   <- loo_compare(waic_l$waic_ctrl1,   
                          waic_l$waic_yr1,     waic_l$waic_yr2,    waic_l$waic_yr3, 
                          waic_l$waic_gaus1,   waic_l$waic_gaus2,  waic_l$waic_gaus3, 
                          waic_l$waic_simpl1,  waic_l$waic_simpl2, waic_l$waic_simpl3, 
                          waic_l$waic_ridge1,  waic_l$waic_ridge2, waic_l$waic_ridge3, 
                          waic_l$waic_gaus,    waic_l$waic_ridge,  waic_l$waic_simpl_n) %>%
                as.data.frame %>%
                tibble::add_column(model = gsub("waic_l\\$waic_","",row.names(.) ), 
                                   .before = 1) %>% 
                # this is useless to me, causes a conflict down the line
                dplyr::select(-elpd_diff) %>% 
                rename( se_diff_waic = se_diff ) %>% 
                mutate( model = gsub('model','',model) %>% as.numeric ) %>% 
                rename( mod_n = model ) %>%   
                left_join( mod_names )  


# leave-one-out crossvalidation ------------------------------------------------------------------------

# crossvalidation function
CrossVal <- function(i, mod_data, response) {       # i is index for row to leave out
  
  # identify years
  uniq_yr           <- mod_data$resp$year %>% unique 
  test_i            <- which(mod_data$resp$year == uniq_yr[i])
  
  # put all in matrix form 
  x_clim            <- mod_data$climate
  x_clim_means      <- rowMeans(mod_data$climate)   # climate averages over entire window (for control model #2)
  
  # response variable
  y_train           <- mod_data$resp[-test_i, response]
  y_test            <- mod_data$resp[test_i, response]
  
  # climate variable
  clim_train        <- x_clim[-test_i,]
  clim_test         <- x_clim[test_i,] 
  
  # climate averages over full 24-month window (for control model #2)
  clim_means_train  <- x_clim_means[-test_i]
  clim_means_test   <- x_clim_means[test_i]
  
  # organize data into list to pass to stan
  dat_stan_crossval <- list(
    n_train = length(y_train),  # number of data points in train set (length of response var)
    n_test  = length(y_test),   # number of data points in test set
    n_lag   = ncol(clim_train), # maximum lag
    y_train = array(y_train),
    y_test  = array(y_test),
    clim    = clim_train,
    clim_train = array(clim_train),
    clim_test  = array(clim_test),
    clim_means_train = array(clim_means_train), # climate averages over full 24-month window (for control model #2)
    clim_means_test  = array(clim_means_test),   # climate averages over full 24-month window (for control model #2)
    clim_yr_train    = list( rowMeans(clim_train[, 1:12]),
                             rowMeans(clim_train[,13:24]),
                             rowMeans(clim_train[,25:36]) ) %>% do.call(rbind,.),
    clim_yr_test     = list( rowMeans(clim_test[, 1:12]),
                             rowMeans(clim_test[,13:24]),
                             rowMeans(clim_test[,25:36]) ) %>% do.call(rbind,.),
    expp_beta = expp_beta, # beta paramater for exponential power distribution
    M         = 12,  # number of months in a year
    K         = ncol(clim_train) / 12,
    # parameters for horseshoe models
    hs_df           = 1,   # variance of 
    hs_df_global    = 1,   # 
    hs_df_slab      = 25,  # slab degrees of freedom
    hs_scale_global = (4 / (36-4)) / sqrt(nrow(mod_data$climate)), # global prior scale
    hs_scale_slab   = 2    # slab prior scale
  )
  
  # fit control 1 (intercept only)
  fit_ctrl1_crossval <- sampling(
    object = readRDS(paste0("R/stan/",family,"_null_crossval.RDS")),
    data   = dat_stan_crossval,
    pars   = c('alpha',  'y_sd', 
               'pred_y', 'log_lik_test'),
    warmup = sim_pars$warmup,
    iter = sim_pars$iter,
    thin = sim_pars$thin,
    chains = sim_pars$chains
  )
  
  # year t
  dat_stan_crossval$clim_means_test   <- rowMeans( mod_data$climate[test_i, 1:12,drop=F] ) %>% array
  dat_stan_crossval$clim_means_train  <- rowMeans( mod_data$climate[-test_i,1:12,drop=F] ) %>% array
  fit_yr1_crossval <- sampling(
    object = readRDS(paste0("R/stan/",family,"_yr_crossval.RDS")),
    data = dat_stan_crossval,
    pars = c('alpha', 'beta', 'y_sd', 
             'pred_y', 'log_lik_test'),
    warmup = sim_pars$warmup,
    iter = sim_pars$iter,
    thin = sim_pars$thin,
    chains = sim_pars$chains
  )
  
  # year t - 1
  dat_stan_crossval$clim_means_test   <- rowMeans( mod_data$climate[test_i, 13:24,drop=F] ) %>% array
  dat_stan_crossval$clim_means_train  <- rowMeans( mod_data$climate[-test_i,13:24,drop=F] ) %>% array
  fit_yr2_crossval <- sampling(
    object = readRDS(paste0("R/stan/",family,"_yr_crossval.RDS")),
    data = dat_stan_crossval,
    pars = c('alpha', 'beta', 'y_sd', 
             'pred_y', 'log_lik_test'),
    warmup = sim_pars$warmup,
    iter = sim_pars$iter,
    thin = sim_pars$thin,
    chains = sim_pars$chains
  )
  
  # year t-2
  dat_stan_crossval$clim_means_test   <- rowMeans( mod_data$climate[test_i, 25:36,drop=F] ) %>% array
  dat_stan_crossval$clim_means_train  <- rowMeans( mod_data$climate[-test_i,25:36,drop=F] ) %>% array
  fit_yr3_crossval <- sampling(
    object = readRDS(paste0("R/stan/",family,"_yr_crossval.RDS")),
    data = dat_stan_crossval,
    pars = c('alpha', 'beta', 'y_sd', 
             'pred_y', 'log_lik_test'),
    warmup = sim_pars$warmup,
    iter = sim_pars$iter,
    thin = sim_pars$thin,
    chains = sim_pars$chains
  )
  
  # fit moving window, gaussian year t
  dat_stan_crossval$n_lag       <- 12
  dat_stan_crossval$clim_train  <- array(clim_train[,1:12])
  dat_stan_crossval$clim_test   <- array(clim_test[,1:12])
  fit_gaus1_crossval <- sampling(
    object = readRDS(paste0("R/stan/",family,"_gaus_crossval.RDS")),
    data = dat_stan_crossval,
    pars = c('sens_mu', 'sens_sd', 'alpha', 'beta', 'y_sd', 
             'pred_y', 'log_lik_test'),
    warmup = sim_pars$warmup,
    iter = sim_pars$iter,
    thin = sim_pars$thin,
    chains = sim_pars$chains,
    control = list(adapt_delta = 0.99)
  )
  
  # fit moving window, gaussian year t - 1
  dat_stan_crossval$clim_train  <- array(clim_train[,13:24])
  dat_stan_crossval$clim_test   <- array(clim_test[,13:24])
  fit_gaus2_crossval <- sampling(
    object = readRDS(paste0("R/stan/",family,"_gaus_crossval.RDS")),
    data = dat_stan_crossval,
    pars = c('sens_mu', 'sens_sd', 'alpha', 'beta', 'y_sd', 
             'pred_y', 'log_lik_test'),
    warmup = sim_pars$warmup,
    iter = sim_pars$iter,
    thin = sim_pars$thin,
    chains = sim_pars$chains,
    control = list(adapt_delta = 0.99)
  )
  
  # fit moving window, gaussian year t - 1
  dat_stan_crossval$clim_train  <- array(clim_train[,25:36])
  dat_stan_crossval$clim_test   <- array(clim_test[,25:36])
  fit_gaus3_crossval <- sampling(
    object = readRDS(paste0("R/stan/",family,"_gaus_crossval.RDS")),
    data = dat_stan_crossval,
    pars = c('sens_mu', 'sens_sd', 'alpha', 'beta', 'y_sd', 
             'pred_y', 'log_lik_test'),
    warmup = sim_pars$warmup,
    iter = sim_pars$iter,
    thin = sim_pars$thin,
    chains = sim_pars$chains,
    control = list(adapt_delta = 0.99)
  )
  
  # Simplex year t
  dat_stan_crossval$clim_train  <- array(clim_train[,1:12]) %>% t
  dat_stan_crossval$clim_test   <- array(clim_test[,1:12]) %>% t
  fit_simpl1_crossval <- sampling(
    object = readRDS(paste0("R/stan/",family,"_dirichlet_crossval.RDS")),
    data = dat_stan_crossval,
    pars = c('theta',  'alpha', 'beta', 'y_sd', 
             'pred_y', 'log_lik_test'),
    warmup = sim_pars$warmup,
    iter = sim_pars$iter,
    thin = sim_pars$thin,
    chains = sim_pars$chains,
    control = list(adapt_delta = 0.99)
  )
  
  # Simplex year t - 1
  dat_stan_crossval$clim_train  <- array(clim_train[,13:24]) %>% t
  dat_stan_crossval$clim_test   <- array(clim_test[,13:24]) %>% t
  fit_simpl2_crossval <- sampling(
    object = readRDS(paste0("R/stan/",family,"_dirichlet_crossval.RDS")),
    data = dat_stan_crossval,
    pars = c('theta',  'alpha', 'beta', 'y_sd', 
             'pred_y', 'log_lik_test'),
    warmup = sim_pars$warmup,
    iter = sim_pars$iter,
    thin = sim_pars$thin,
    chains = sim_pars$chains,
    control = list(adapt_delta = 0.99)
  )
  
  # Simplex year t - 2
  dat_stan_crossval$clim_train  <- array(clim_train[,25:36]) %>% t
  dat_stan_crossval$clim_test   <- array(clim_test[,25:36]) %>% t
  fit_simpl3_crossval <- sampling(
    object = readRDS(paste0("R/stan/",family,"_dirichlet_crossval.RDS")),
    data = dat_stan_crossval,
    pars = c('theta',  'alpha', 'beta', 'y_sd', 
             'pred_y', 'log_lik_test'),
    warmup = sim_pars$warmup,
    iter = sim_pars$iter,
    thin = sim_pars$thin,
    chains = sim_pars$chains,
    control = list(adapt_delta = 0.99)
  )
  
  # Ridge year t 
  dat_stan_crossval$clim_train  <- array(clim_train[,1:12])
  dat_stan_crossval$clim_test   <- array(clim_test[,1:12])
  fit_ridge1_crossval <- sampling(
    object = readRDS(paste0("R/stan/",family,"_horse_crossval.RDS")),
    data = dat_stan_crossval,
    pars = c('alpha',  'beta', 'y_sd', 
             'pred_y', 'log_lik_test'),
    warmup = sim_pars$warmup,
    iter = sim_pars$iter,
    thin = sim_pars$thin,
    chains = sim_pars$chains,
    control = list(adapt_delta = 0.99)
  )
  
  # Ridge year t - 1
  dat_stan_crossval$clim_train  <- array(clim_train[,13:24])
  dat_stan_crossval$clim_test   <- array(clim_test[,13:24])
  fit_ridge2_crossval <- sampling(
    object = readRDS(paste0("R/stan/",family,"_horse_crossval.RDS")),
    data = dat_stan_crossval,
    pars = c('alpha',  'beta', 'y_sd', 
             'pred_y', 'log_lik_test'),
    warmup = sim_pars$warmup,
    iter = sim_pars$iter,
    thin = sim_pars$thin,
    chains = sim_pars$chains,
    control = list(adapt_delta = 0.99)
  )
  
  # Ridge year t - 2
  dat_stan_crossval$clim_train  <- array(clim_train[,25:36])
  dat_stan_crossval$clim_test   <- array(clim_test[,25:36])
  fit_ridge3_crossval <- sampling(
    object = readRDS(paste0("R/stan/",family,"_horse_crossval.RDS")),
    data = dat_stan_crossval,
    pars = c('alpha',  'beta', 'y_sd', 
             'pred_y', 'log_lik_test'),
    warmup = sim_pars$warmup,
    iter = sim_pars$iter,
    thin = sim_pars$thin,
    chains = sim_pars$chains,
    control = list(adapt_delta = 0.99)
  )
  
  # fit moving window, gaussian year t through t-2
  dat_stan_crossval$n_lag       <- 36
  dat_stan_crossval$clim_train  <- array(clim_train)
  dat_stan_crossval$clim_test   <- array(clim_test)
  fit_gaus_all_crossval <- sampling(
    object = readRDS(paste0("R/stan/",family,"_gaus_crossval.RDS")),
    data = dat_stan_crossval,
    pars = c('sens_mu', 'sens_sd', 'alpha', 'beta', 'y_sd', 
             'pred_y', 'log_lik_test'),
    warmup = sim_pars$warmup,
    iter = sim_pars$iter,
    thin = sim_pars$thin,
    chains = sim_pars$chains,
    control = list(adapt_delta = 0.99)
  )
  
  # Ridge year t through t-2
  fit_ridge_all_crossval <- sampling(
    object = readRDS(paste0("R/stan/",family,"_horse_crossval.RDS")),
    data = dat_stan_crossval,
    pars = c('alpha',  'beta', 'y_sd', 
             'pred_y', 'log_lik_test'),
    warmup = sim_pars$warmup,
    iter = sim_pars$iter,
    thin = sim_pars$thin,
    chains = sim_pars$chains,
    control = list(adapt_delta = 0.99)
  )
  
  
  # update data list
  dat_stan_crossval$clim_train        <- t(dat_stan_crossval$clim_train)
  dat_stan_crossval$clim_test         <- t(dat_stan_crossval$clim_test)
  
  dat_stan_crossval$clim1_train       <- dat_stan_crossval$clim_train[1:12,]
  dat_stan_crossval$clim1_test        <- dat_stan_crossval$clim_test[ 1:12, ,drop=F]
  
  dat_stan_crossval$clim2_train       <- dat_stan_crossval$clim_train[13:24,]
  dat_stan_crossval$clim2_test        <- dat_stan_crossval$clim_test[ 13:24, ,drop=F]
  
  dat_stan_crossval$clim3_train       <- dat_stan_crossval$clim_train[25:36,]
  dat_stan_crossval$clim3_test        <- dat_stan_crossval$clim_test[ 25:36, ,drop=F]
  
  dat_stan_crossval$clim1_means_train <- colMeans(dat_stan_crossval$clim_train[1:12,])
  dat_stan_crossval$clim1_means_test  <- colMeans(dat_stan_crossval$clim_test[ 1:12, ,drop=F]) %>% array
  
  dat_stan_crossval$clim2_means_train <- colMeans(dat_stan_crossval$clim_train[13:24,])
  dat_stan_crossval$clim2_means_test  <- colMeans(dat_stan_crossval$clim_test[ 13:24, ,drop=F]) %>% array
  
  dat_stan_crossval$clim3_means_train <- colMeans(dat_stan_crossval$clim_train[25:36,])
  dat_stan_crossval$clim3_means_test  <- colMeans(dat_stan_crossval$clim_test[ 25:36, ,drop=F]) %>% array
  
  # fit simplex nested within year
  fit_36_nest_crossval <- sampling(
    object = readRDS(paste0("R/stan/",family,"_dirichlet_nest_crossval.RDS")),
    data = dat_stan_crossval,
    pars = c('theta_m', "theta_y",'alpha', 'beta', 'y_sd', 
             'pred_y', 'log_lik_test'),
    warmup = sim_pars$warmup,
    iter = sim_pars$iter,
    thin = sim_pars$thin,
    chains = sim_pars$chains,
    control = list(adapt_delta = 0.99)
  )
  
  # posterior mean prediction for the out-of-sample value
  crossval_mods <- list( ctrl1   = fit_ctrl1_crossval, 
                         yr1     = fit_yr1_crossval, 
                         yr2     = fit_yr2_crossval, 
                         yr3     = fit_yr3_crossval, 
                         gaus1   = fit_gaus1_crossval, 
                         gaus2   = fit_gaus2_crossval, 
                         gaus3   = fit_gaus3_crossval, 
                         simpl1  = fit_simpl1_crossval,
                         simpl2  = fit_simpl2_crossval,
                         simpl3  = fit_simpl3_crossval,
                         ridge1  = fit_ridge1_crossval,
                         ridge2  = fit_ridge2_crossval,
                         ridge3  = fit_ridge3_crossval,  
                         
                         gaus    = fit_gaus_all_crossval,  
                         ridge   = fit_ridge_all_crossval,  
                         simpl_n = fit_36_nest_crossval )
  
  
  # mean predictions
  mod_preds <- lapply(crossval_mods, function(x) rstan::extract(x, 'pred_y')$pred_y %>% apply(2,mean) )
  
  # Expected Log Predictive Density
  mod_elpds <- lapply(crossval_mods, function(x){
    rstan::extract(x, 'log_lik_test')$log_lik_test %>% 
      exp %>%
      apply(2,mean) %>%
      log 
  } )
  
  # mean posterior of predictions
  pred_post <- function(x, test_i ){
    
    # first 
    post_raw <- rstan::extract(x, 'pred_y') %>% 
      .$pred_y %>% 
      as.data.frame %>% 
      stack
    
    # update index values
    updt_df  <- data.frame( ind    = unique(post_raw$ind) %>% sort,
                            test_i = unlist(test_i) )
    
    suppressMessages( left_join( post_raw, updt_df ) ) %>%
      dplyr::select(test_i, values )
    
  }
  
  # posterior of predictions, with formatting
  post_raw_l <- lapply(crossval_mods, pred_post, list(test_i) ) 
  post_df    <- lapply( 1:length(post_raw_l), function(ii) post_raw_l[[ii]] %>% 
                          mutate( mod = names(post_raw_l)[i] )
  ) %>% 
    bind_rows
  
  # diagnostics 
  diagnostics <- function(fit_obj, name_mod){
    
    diverg      <- do.call(rbind, args = get_sampler_params(fit_obj, inc_warmup = F))[,5]
    n_diverg    <- length(which(diverg == 1))
    df_summ     <- as.data.frame(summary(fit_obj)$summary)
    rhat_high   <- length(which(df_summ$Rhat > 1.1))
    n_eff       <- df_summ$n_eff / length(diverg)
    n_eff_low   <- length(which(n_eff < 0.1))
    mcse_high   <- length(which(df_summ$se_mean / df_summ$sd > 0.1))
    out         <- data.frame(n_diverg,  rhat_high, 
                              n_eff_low, mcse_high) 
    # out         <- setNames(out, paste0(names(out),"_",name_mod) )
    return(out)
    
  }
  
  # store diagnostics
  # gaus_expp   <- crossval_mods[c("gaus","expp","gev","simpl")]
  diagnost_l  <- Map(diagnostics, crossval_mods, names(crossval_mods))
  diagnost_df <- do.call(cbind, diagnost_l) %>%
    bind_cols( unique( dplyr::select(mod_data$resp[test_i,],year) ) ) %>% 
    setNames( gsub('\\.','_',names(.)) )
  
  # function
  pred_elpd_df<- mod_data$resp[test_i,] %>%
    mutate( # predictions
      ctrl1_pred   = mod_preds$ctrl1,
      yr1_pred     = mod_preds$yr1,
      yr2_pred     = mod_preds$yr2,
      yr3_pred     = mod_preds$yr3,
      gaus1_pred   = mod_preds$gaus1,
      gaus2_pred   = mod_preds$gaus2,
      gaus3_pred   = mod_preds$gaus3,
      simpl1_pred  = mod_preds$simpl1,
      simpl2_pred  = mod_preds$simpl2,
      simpl3_pred  = mod_preds$simpl3,
      ridge1_pred  = mod_preds$ridge1,
      ridge2_pred  = mod_preds$ridge2,
      ridge3_pred  = mod_preds$ridge3,
      gaus_pred    = mod_preds$gaus,
      ridge_pred   = mod_preds$ridge,
      simpl_n_pred = mod_preds$simpl_n,      
      
      # Expected Log Predictive Density
      ctrl1_elpd   = mod_elpds$ctrl1,
      yr1_elpd     = mod_elpds$yr1,
      yr2_elpd     = mod_elpds$yr2,
      yr3_elpd     = mod_elpds$yr3,
      gaus1_elpd   = mod_elpds$gaus1,
      gaus2_elpd   = mod_elpds$gaus2,
      gaus3_elpd   = mod_elpds$gaus3,
      simpl1_elpd  = mod_elpds$simpl1,
      simpl2_elpd  = mod_elpds$simpl2,
      simpl3_elpd  = mod_elpds$simpl3,
      ridge1_elpd  = mod_elpds$ridge1,
      ridge2_elpd  = mod_elpds$ridge2,
      ridge3_elpd  = mod_elpds$ridge3,
      gaus_elpd    = mod_elpds$gaus,
      ridge_elpd   = mod_elpds$ridge,
      simpl_n_elpd = mod_elpds$simpl_n )
  
  # df to return
  out         <- left_join(pred_elpd_df, diagnost_df)
  
  # remove stanfit objects (garbage collection)
  rm(fit_ctrl1_crossval) 
  rm(fit_yr1_crossval) 
  rm(fit_yr2_crossval) 
  rm(fit_yr3_crossval) 
  rm(fit_gaus1_crossval) 
  rm(fit_gaus2_crossval) 
  rm(fit_gaus3_crossval)
  rm(fit_simpl1_crossval) 
  rm(fit_simpl2_crossval)
  rm(fit_simpl3_crossval)
  rm(fit_ridge1_crossval)
  rm(fit_ridge2_crossval)
  rm(fit_ridge3_crossval)
  rm(fit_gaus_all_crossval)
  rm(fit_ridge_all_crossval)
  rm(fit_36_nest_crossval)
  
  return( list('means' = out, 'post' = post_df) )
  
}

# spp-specific cross validation
year_inds   <- seq_along(unique(mod_data$resp$year))
cxval_res_ll<- lapply( year_inds, CrossVal, mod_data, response)

# extra formatting
pluck_res   <- function (i, what_i) cxval_res_ll[[i]] %>% pluck(what_i)

# put means/diagnostics, and posterior in separate lists
mean_l      <- lapply( year_inds, pluck_res, 1 )
post_l      <- lapply( year_inds, pluck_res, 2 )

# means/diagnostics as data frames
cxval_pred  <- do.call(rbind, mean_l)
cxval_post  <- do.call(rbind, post_l) %>% arrange(test_i, mod)


# measures of fit -------------------------------------------------------------------------- 

# calculate either mse or deviance
pred_perform <- function(x, mod_data, response, type){
  
  if( type == "mse"){
    res   <- (x - mod_data$resp[,response])^2 %>% mean
  }
  if(type == "deviance"){
    res   <-calc.deviance(x, mod_data$resp[,response],
                          weights = rep(1, length(x) ),
                          family="gaussian", calc.mean = TRUE)
  }
  
  return(res)
  
}

# format results into a data frame
perform_format <- function(x, var){
  
  x %>%
    unlist %>%
    t %>%
    t %>%
    as.data.frame %>%
    tibble::rownames_to_column(var = "model") %>%
    mutate( model = gsub("mod_preds.", "", model))  %>%
    setNames( c("model", var) )
  
}

# order of 'mod_data$resp' and 'cxval_pred' need be the same
expect_true( all.equal(dplyr::select(mod_data$resp, year, population),
                       dplyr::select(cxval_pred,    year, population)) )

print('IMPORTANT: Check names')
names(cxval_pred)

# mean squared error
mse <- cxval_pred %>%
  dplyr::select( grep('_pred',names(.),value=T) ) %>%
  lapply(pred_perform, mod_data, response, "mse") %>%
  perform_format("mse") %>%
  mutate( model = gsub("_pred","",model) )

# Expected Log Predictive Density
elpd <- cxval_pred %>%
  dplyr::select( grep('_elpd',names(.),value=T) ) %>%
  apply(2, sum) %>%
  as.matrix %>%
  as.data.frame %>%
  tibble::add_column(model = rownames(.), .before=1) %>%
  mutate( model = gsub("_elpd","",model) ) %>%
  setNames( c("model", "elpd") )

# measures of fit
# mof  <- merge(mse, devi)
mof  <- merge(mse, elpd)

# store results ---------------------------------------------------------------------------
mod_summs <- Reduce(function(...) merge(...), 
                    list(mod_pars_diag, loo_df, waic_df) ) %>% #, mof
                    arrange( mse )

write.csv(mod_summs,  paste0(args[3], "_mod_summaries_",spp_name,".csv"), row.names = F)
write.csv(posteriors, paste0(args[3], "_posterior_",spp_name,".csv"), row.names = F)
write.csv(cxval_pred, paste0(args[3], "_crossval_pred_diag_",spp_name,".csv"), row.names = F)
write.csv(cxval_post, paste0(args[3], "_crossval_pred_post_",spp_name,".csv"), row.names = F)
