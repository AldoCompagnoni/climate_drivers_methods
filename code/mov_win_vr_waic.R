# bjtwd<-"C:/Users/admin_bjt162/Dropbox/A.Current/Ongoing_Collab_Research/sApropos project/"
rm(list=ls())
source("code/format_data.R")
library(tidyverse)
library(dismo)
library(mgcv)
library(testthat)
library(rstan)
library(loo)

# set rstan options
rstan_options( auto_write = TRUE )
options( mc.cores = 4 )

# climate predictor, response, months back, max. number of knots
response  <- "surv"
clim_var  <- "precip"
m_back    <- 36    
st_dev    <- FALSE

for(ii in 23:34){

# read data -----------------------------------------------------------------------------------------
lam       <- read.csv("data/all_demog_6tr.csv", stringsAsFactors = F)
clim      <- data.table::fread(paste0('data/',clim_var,"_chelsa_hays.csv"),  
                               stringsAsFactors = F)
# clim      <- data.table::fread(paste0(clim_var,"_fc_hays.csv"),  stringsAsFactors = F)
spp       <- lam$SpeciesAuthor %>% unique

# format data --------------------------------------------------------------------------------------

# set up model "family" based on response
if( response == "surv" | response == "grow" )             family = "beta" 
if( response == "fec" )                                   family = "gamma"
if( grepl("PreRep", response) | grepl("Rep", response) )  family = "beta"
if( response == "rho" | response == "react_fsa" )         family = "gamma"
if( response == "log_lambda")                             family = "normal"

expp_beta     <- 20

# set species (I pick Sphaeraclea_coccinea)
ii            <- 30
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

# transform survival/growth - ONLY if less than 30% data points are 1/0
if( response == "surv" | response == "grow" | grepl("PreRep", response) | grepl("Rep", response) ){
  
  raw_x <- mod_data$resp[,response]
  pc_1  <- sum( raw_x == 1 ) / length(raw_x)
  pc_0  <- sum( raw_x == 0 ) / length(raw_x)
  
  # for survival
  if( grepl("surv", response, ignore.case = T) & pc_1 < 0.3 ){
    n     <- length(raw_x)
    new_x <- ( raw_x*(n - 1) + 0.5 ) / n
    mod_data$resp[,response] <- new_x
  }
  
  # for growth
  if( grepl("grow", response, ignore.case = T) & pc_0 < 0.3 ){
    n     <- length(raw_x)
    new_x <- ( raw_x*(n - 1) + 0.5 ) / n
    mod_data$resp[,response] <- new_x
  }
  
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
                  rowMeans(mod_data$climate[,25:36]) ) %>% do.call(rbind, .),
  M       = 12,    # number of months in a year
  K       = ncol(mod_data$climate) / 12,
  S       = mod_data$resp$population %>% unique %>% length,
  site_i  = mod_data$resp$population %>% as.factor %>% as.numeric,
  expp_beta = expp_beta
)

# simulation parameters
sim_pars <- list(
  warmup = 1000, 
  iter = 4000, 
  thin = 2, 
  chains = 1
)

# NULL model (model of the mean)
fit_ctrl1 <- stan(
  file = paste0("code/stan/",family,"_null.stan"),
  data = dat_stan,
  pars = c('alpha', 'y_sd', 'log_lik'),
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains
)

# year t
dat_stan$clim_means  <- rowMeans(mod_data$climate[,1:12 ])
fit_yr1 <- stan(
  file = paste0("code/stan/",family,"_yr.stan"),
  data = dat_stan,
  pars = c('alpha', 'beta', 'y_sd', 'yhat', 'log_lik'),
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains,
  control = list(adapt_delta = 0.99) #, stepsize = 0.001, max_treedepth = 20)
)

# year t-1
dat_stan$clim_means  <- rowMeans(mod_data$climate[,13:24])
fit_yr2 <- stan(
  file = paste0("code/stan/",family,"_yr.stan"),
  data = dat_stan,
  pars = c('alpha', 'beta', 'y_sd', 'yhat', 'log_lik'),
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains,
  control = list(adapt_delta = 0.99)
)

# year t-2
dat_stan$clim_means  <- rowMeans(mod_data$climate[,25:36])
fit_yr3 <- stan(
  file = paste0("code/stan/",family,"_yr.stan"),
  data = dat_stan,
  pars = c('alpha', 'beta', 'y_sd', 'yhat', 'log_lik'),
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains,
  control = list(adapt_delta = 0.99)
)
dat_stan$clim_means  <- rowMeans(mod_data$climate)

# gaussian moving window
fit_gaus <- stan(
  file = paste0("code/stan/",family,"_gaus.stan"),
  data = dat_stan,
  pars = c('sens_mu', 'sens_sd', 'alpha', 'beta', 'y_sd',
           'yhat',  'log_lik'), 
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains,
  control = list(adapt_delta = 0.99)
)

# exponential power moving window
fit_expp <- stan(
  file = paste0("code/stan/",family,"_expp.stan"),
  data = dat_stan,
  pars = c('sens_mu', 'sens_sd', 'alpha', 'beta', 
           'y_sd', 'yhat', 'log_lik'), 
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains,
  control = list(adapt_delta = 0.99)
)

# moving beta hierarchical
fit_mb_h <- stan(
  file = paste0("code/stan/",family,"_movbeta_h.stan"),
  data = dat_stan,
  pars = c('alpha', 'beta', 'y_sd', 
           'mu_beta', 'sigma_beta', 
           'yhat', 'log_lik'),
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains,
  control = list(adapt_delta = 0.99)
)

# moving beta correlated
fit_mb <- stan(
  file = paste0("code/stan/",family,"_movbeta.stan"),
  data = dat_stan,
  pars = c('alpha',   'beta',  'y_sd', 'mu_beta', 
           'yhat', 'log_lik'),
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains,
  control = list(adapt_delta = 0.99)
)

# moving beta hierarchical
fit_mb_h_nest <- stan(
  file = paste0("code/stan/",family,"_movbeta_h_nest.stan"),
  data = dat_stan,
  pars = c('alpha', 'beta', 'y_sd', 
           'mu_beta', 'sigma_beta', 'theta_y',
           'yhat', 'log_lik'),
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains,
  control = list(adapt_delta = 0.99)
)

# moving beta correlated
fit_mb_nest <- stan(
  file = paste0("code/stan/",family,"_movbeta_nest.stan"),
  data = dat_stan,
  pars = c('alpha', 'beta', 'y_sd',
           'mu_beta', 'theta_y', 
           'yhat', 'log_lik'),
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains,
  control = list(adapt_delta = 0.99)
)

# Nested models 
# update data list
dat_stan$clim         <- t(mod_data$climate)
dat_stan$clim1        <- t(mod_data$climate)[1:12 ,]
dat_stan$clim2        <- t(mod_data$climate)[13:24,]
dat_stan$clim3        <- t(mod_data$climate)[25:36,]

# Generalized Extreme Value nested
fit_gaus_nest <- stan(
  file = paste0("code/stan/",family,"_gaus_nest.stan"),
  data = dat_stan,
  pars = c('sens_mu', 'sens_sd', 'theta_y', 'alpha', 'beta', 'y_sd', 
           'yhat', 'log_lik'),
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains,
  control = list(adapt_delta = 0.99)
)

# Power exponential nested 
fit_expp_nest <- stan(
  file = paste0("code/stan/",family,"_expp_nest.stan"),
  data = dat_stan,
  pars = c('sens_mu', 'sens_sd', 'theta_y', 'alpha', 'beta', 'y_sd', 
           'yhat', 'log_lik'),
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains,
  control = list(adapt_delta = 0.99)
)

# Simplex nested
fit_24_nest <- stan(
  file = paste0("code/stan/",family,"_dirichlet_nest.stan"),
  data = dat_stan,
  pars = c('theta_y', 'theta_m', 'alpha', 'beta', 'y_sd', 
           'yhat', 'log_lik'),
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains,
  control = list(adapt_delta = 0.99)
)


# parameter values and diagnostics ----------------------------------------------------------------

# list of model fits
mod_fit   <- list(ctrl1   = fit_ctrl1,     
                  yr1     = fit_yr1,     
                  yr2     = fit_yr2,
                  yr3     = fit_yr3,
                  gaus    = fit_gaus,     expp     = fit_expp,
                  movb_h  = fit_mb_h,     movb     = fit_mb,      
                  movb_h_n = fit_mb_h_nest,
                  movb_n   = fit_mb_nest,  
                  simpl_n  = fit_24_nest, 
                  gaus_n   = fit_gaus_nest,
                  expp_n   = fit_expp_nest )

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

# WAIC model comparison --------------------------------------------------------------------

# wAIC model selection using loo approximation (from library 'loo')
log_liks   <- lapply(mod_fit, extract_log_lik)

# leave-one-out estimates
loo_l      <- lapply(log_liks, loo) %>%
                setNames( c("loo_ctrl1",   
                            "loo_yr1",     "loo_yr2",   "loo_yr3", 
                            "loo_gaus",    "loo_expp",   
                            "loo_movb",    "loo_movb_c",   
                            "loo_simpl_n", "loo_gaus_n", "loo_expp_n") )
loo_df     <- loo::compare(loo_l$loo_ctrl1,   
                           loo_l$loo_yr1,     loo_l$loo_yr2,   loo_l$loo_yr3, 
                           loo_l$loo_gaus,    loo_l$loo_expp,  
                           loo_l$loo_movb,    loo_l$loo_movb_c,  
                           loo_l$loo_simpl_n, loo_l$loo_gaus_n, loo_l$loo_expp_n ) %>%
                as.data.frame %>%
                tibble::add_column(model = gsub("loo_l\\$loo_","",row.names(.) ), .before = 1)

# WAIC estimates
waic_l    <- lapply(log_liks, waic) %>%
                setNames(c("waic_ctrl1",   
                           "waic_yr1",     "waic_yr2",   "waic_yr3",  
                           "waic_gaus",    "waic_expp", 
                           "waic_movb",    "waic_movb_c", 
                           "waic_simpl_n", "waic_gaus_n", "waic_expp_n") )
waic_df   <- loo::compare(waic_l$waic_ctrl1,   
                          waic_l$waic_yr1,     waic_l$waic_yr2,     waic_l$waic_yr3,
                          waic_l$waic_gaus,    waic_l$waic_expp, 
                          waic_l$waic_movb,    waic_l$waic_movb_c, 
                          waic_l$waic_simpl_n, waic_l$waic_gaus_n,   waic_l$waic_expp_n) %>%
                as.data.frame %>%
                tibble::add_column(model = gsub("waic_l\\$waic_","",row.names(.) ), 
                                   .before = 1) %>% 
                # this is useless to me, causes a conflict down the line
                dplyr::select(-elpd_diff)

# store results ---------------------------------------------------------------------------
mod_summs <- Reduce(function(...) merge(...), list(mod_pars_diag, loo_df, waic_df) )

write.csv(mod_summs,  paste0('results/',clim_var, "/mod_summaries_",spp_name,".csv"), row.names = F)
write.csv(posteriors, paste0('results/',clim_var, "/posterior_",spp_name,".csv"), row.names = F)
# write.csv(cxval_pred, paste0('results/',clim_var, "/crossval_pred_diag_",spp_name,".csv"), row.names = F)

}

