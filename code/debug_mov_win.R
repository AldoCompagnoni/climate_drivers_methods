#bjtwd<-"C:/Users/admin_bjt162/Dropbox/A.Current/Ongoing_Collab_Research/sApropos project/"
rm(list=ls())
setwd("C:/cloud/Dropbox/sApropos/")
source("C:/CODE/moving_windows/format_data.R")
library(tidyverse)
library(dismo)
library(mgcv)
library(testthat)
library(rstan)
library(loo)

# set rstan options
rstan_options( auto_write = TRUE )
options( mc.cores = parallel::detectCores() )

# "pipeable" Reduce rbind
rbind_l <- function(x) Reduce(function(...) rbind(...), x)

# climate predictor, response, months back, max. number of knots
response  <- "log_lambda"
clim_var  <- "airt"
m_back    <- 36    
st_dev    <- FALSE

# read data -----------------------------------------------------------------------------------------
lam       <- read.csv("all_demog_6tr.csv", stringsAsFactors = F)
m_info    <- read.csv("MatrixEndMonth_information.csv", stringsAsFactors = F)
clim      <- data.table::fread(paste0(clim_var,"_chelsa_hays.csv"),  stringsAsFactors = F)
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
                  rowMeans(mod_data$climate[,25:36]) ) %>% rbind_l,
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
  chains = 3
)

# NULL model (model of the mean)
fit_ctrl1 <- stan(
  file = paste0("stan/",family,"_null.stan"),
  data = dat_stan,
  pars = c('alpha', 'y_sd', 'log_lik'),
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains
  #control = list(adapt_delta = 0.999, stepsize = 0.001, max_treedepth = 20)
)

# average climate model
fit_ctrl2 <- stan(
  file = paste0("stan/",family,"_ctrl2.stan"),
  data = dat_stan,
  pars = c('alpha', 'beta', 'y_sd', 'log_lik'),
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains
  #control = list(adapt_delta = 0.999, stepsize = 0.001, max_treedepth = 20)
)


# year t
dat_stan$clim_means  <- rowMeans(mod_data$climate[,1:12 ])
fit_yr1 <- stan(
  file = paste0("stan/",family,"_yr.stan"),
  data = dat_stan,
  pars = c('alpha', 'beta', 'y_sd', 'log_lik'),
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains#,
  #control = list(adapt_delta = 0.999, stepsize = 0.001, max_treedepth = 20)
)

# year t-1
dat_stan$clim_means  <- rowMeans(mod_data$climate[,13:24])
fit_yr2 <- stan(
  file = paste0("stan/",family,"_yr.stan"),
  data = dat_stan,
  pars = c('alpha', 'beta', 'y_sd', 'log_lik'),
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains#,
  #control = list(adapt_delta = 0.999, stepsize = 0.001, max_treedepth = 20)
)

# year t-2
dat_stan$clim_means  <- rowMeans(mod_data$climate[,25:36])
fit_yr3 <- stan(
  file = paste0("stan/",family,"_yr.stan"),
  data = dat_stan,
  pars = c('alpha', 'beta', 'y_sd', 'log_lik'),
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains#,
  #control = list(adapt_delta = 0.999, stepsize = 0.001, max_treedepth = 20)
)
dat_stan$clim_means  <- rowMeans(mod_data$climate)

# year weights
fit_yr_weight <- stan(
  file = paste0("stan/",family,"_yr_dirichlet.stan"),
  data = dat_stan,
  pars = c('theta', 'alpha', 'beta', 'y_sd', 'log_lik'),
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains#,
  #control = list(adapt_delta = 0.999, stepsize = 0.001, max_treedepth = 20)
)

# year beta (a different beta for each year)
fit_yr_beta <- stan(
  file = paste0("stan/",family,"_yr_beta.stan"),
  data = dat_stan,
  pars = c('alpha', 'beta', 'y_sd', 'log_lik'),
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains#,
  #control = list(adapt_delta = 0.999, stepsize = 0.001, max_treedepth = 20)
)

# gaussian moving window
fit_gaus <- stan(
  file = paste0("stan/",family,"_gaus.stan"),
  data = dat_stan,
  pars = c('sens_mu', 'sens_sd', 'alpha', 'beta', 'y_sd', 'log_lik'),
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains,
  control = list( adapt_delta = 0.99)
  #control = list(adapt_delta = 0.999, stepsize = 0.001, max_treedepth = 20)
)


# exponential power moving window
fit_expp <- stan(
  file = paste0("stan/",family,"_expp.stan"),
  data = dat_stan,
  pars = c('sens_mu', 'sens_sd', 'alpha', 'beta', 'y_sd', 'log_lik'),
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains,
  control = list( adapt_delta = 0.99)
  #control = list(adapt_delta = 0.999, stepsize = 0.001, max_treedepth = 20)
)

# Generalized extreme value 
fit_gev <- stan(
  file = paste0("stan/",family,"_gev.stan"),
  data = dat_stan,
  pars = c('loc', 'scale', "shape", 'alpha', 'beta', 'y_sd', 'log_lik'),
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains#,
  #control = list(adapt_delta = 0.999, stepsize = 0.001, max_treedepth = 20)
)

# Nested models 
# update data list
dat_stan$clim         <- t(mod_data$climate)
dat_stan$clim1        <- t(mod_data$climate)[1:12 ,]
dat_stan$clim2        <- t(mod_data$climate)[13:24,]
dat_stan$clim3        <- t(mod_data$climate)[25:36,]

# Simplex nested
fit_24_nest <- stan(
  file = paste0("stan/",family,"_dirichlet_nest.stan"),
  data = dat_stan,
  pars = c('theta_y', 'theta_m', 'alpha', 'beta', 'y_sd', 'log_lik'),
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains,
  control = list(adapt_delta = 0.99)
)


# try to debug -----------------------------------------------

get_diverg <- function(x){
  
  sampl_df <- get_sampler_params(x, inc_warmup = F) %>% 
                lapply(as.data.frame) %>% 
                bind_rows 
            
  sampl_df %>% subset( divergent__ == 1 ) %>% nrow

}

get_diverg(fit_expp)

params <- fit_gaus %>% 
            rstan::extract(permuted=FALSE) %>% 
            as.data.frame %>% 
            select( -grep('log_lik',names(.)) ) %>% 
            setNames( gsub('chain:','',names(.)) ) %>% 
            # point out iterations
            mutate( iter = 1:nrow(.) )



par(mar = c(4, 4, 0.5, 0.5))
plot(params$iter, params$`2.sens_mu`, 
     pch=16, cex=0.8,
     xlab="Iteration", ylab="mu")

par(mar = c(4, 4, 0.5, 0.5))
plot(params$`iter`, params$`3.sens_sd`, 
     pch=16, cex=0.8,
     xlab="sens_mu", ylab="sd")



params24 <- fit_24_nest %>% 
              rstan::extract() %>%
              # rstan::extract(permuted=FALSE) %>%
              as.data.frame %>% 
              select( -grep('log_lik',names(.)) ) %>% 
              # point out iterations
              mutate( iter = 1:nrow(.) )

plot(`chain:1.beta` ~ iter,data=params24)
plot(theta_y.2 ~ iter,data=params24)
plot(`chain:3.beta` ~ iter,data=params24)          
plot(theta_y.1 ~ theta_y.3,data=params24)

params24 %>% 
  select( grep('theta_m',names(.)) ) %>% 
  boxplot

params24 %>% 
  select( grep('theta_y',names(.)) ) %>% 
  boxplot


# calculate running mean 
running_means <- sapply(params$iter, 
                        function(n) mean(params$`2.sens_mu`[1:n]))
par(mar = c(4, 4, 0.5, 0.5))
plot(params$iter, running_means, 
     pch=16, cex=0.8,
      xlab="Iteration", ylab="MCMC mean of mu")

# Count divergent transitions
divergent <- get_sampler_params(fit_gaus, inc_warmup=T)[[1]][,'divergent__']
sum(divergent)


# 