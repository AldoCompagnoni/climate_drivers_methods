# simulations to 
# retrieve real parameters for a matrix of simulated models based on
# 1. weighted or non-weighted yearly models
# 2. spatially replicated or not
# 3. sd of normal response: 0.3 or 0.5
# 4. beta of effect: 0.45 or 1.2
source("C:/CODE/moving_windows/format_data.R")
library(shinystan)
library(dplyr)
library(tidyr)
library(ggplot2)
library(viridis)
library(testthat)
library(rstan)
library(evd)  
library(purrr)
library(rmutil)
library(parallel)

# set rstan options
rstan_options( auto_write = TRUE )
options( mc.cores = parallel::detectCores() )


# simulations that need be tested
sim_needed <- expand.grid( weight = c(0,    1),
                           spat_n = c(1,    5),
                           sd     = c(0.3,  0.5),
                           beta   = c(0.45, 1.2) )

# calculate realistic y_sd from data ------------------------------------

# get lambda data from Aldo's personal files
lam       <- read.csv("C:/cloud/Dropbox/sApropos/all_demog_6tr.csv", 
                      stringsAsFactors = F) %>% 
                # remove most dodgy datasets
                subset( !(SpeciesAuthor %in% c('Trillium_ovatum',
                                             'Opuntia_macrorhiza_2',
                                             'Cirsium_pitcheri_4',
                                             'Silene_spaldingii') ) )

# SD per study
sd_df <- lam %>% 
            group_by(SpeciesAuthor) %>% 
            summarise( sd_ll  = sd(log_lambda),
                       rep    = n(),
                       rep_yr = MatrixEndYear %>% 
                                unique %>% 
                                length ) %>% 
            ungroup


plot(sd_ll ~ rep, data=sd_df)
plot(sd_ll ~ rep_yr, data=sd_df)

# sd varies from 0.017 to 1.55. Mean 0.36, Median 0.24
sd_df$sd_ll %>% mean
sd_df$sd_ll %>% median
sd_df$sd_ll %>% max
sd_df$sd_ll %>% min

# format raw climate data (messy code) ----------------------------------------

# read 
clim_x    <- read.csv('C:/CODE/climate_drivers_methods/data/demo_airt.csv')

# format: unfortunately this data is replicated daily to accommodate potential daily data
day_one   <- as.Date( paste0("1/1/", first(clim_x$year) ), 
                    format="%d/%m/%Y" ) 

# introduce monthly information
clim_d    <- as.Date(1:nrow(clim_x), day_one-1) %>%
              as.character %>%
              as.data.frame(stringsAsFactors=F) %>%
              separate_(col=".",into=c("year1","month1","day1"),sep="-") %>%
              bind_cols(clim_x) %>%
              dplyr::select(-year,-day) %>%
              setNames( c("year", "month", "day", "species" ,"population", "value") )

# test there is 1 value per month only
clim_d %>% 
  dplyr::select(year,month,value) %>% 
  unique %>% 
  group_by(year,month) %>% 
  summarise(rep=n()) %>% 
  .$rep %>% 
  unique %>% 
  expect_equal( 1 )


# take unique monthly values and calculate anomalies
anom_c <- clim_d %>% 
            dplyr::select(year,month,value) %>% 
            unique %>% 
            # format data in wide form
            spread( month, value ) %>% 
            # calculate anomalies
            select( -year ) %>% 
            apply(2, FUN = scale, center = T, scale = T) %>% 
            # format as data frame
            as.data.frame

# bind three years of data back
yr_1    <- cbind( year = c(1979:2013), anom_c )
yr_2    <- cbind( year = c(1979:2013), anom_c ) %>% 
              setNames( c('year', 13:24) ) %>% 
              mutate( year = year + 1 )
yr_3    <- cbind( year = c(1979:2013), anom_c ) %>% 
              setNames( c('year', 25:36) ) %>% 
              mutate( year = year + 2 )

# Order the anomalies
anom_a  <- Reduce( function(...) inner_join(...), 
                   list(yr_1, yr_2, yr_3) ) %>% 
              .[,c('year', c( paste0(12:10),
                              paste0('0',9:1), 
                              paste0(24:13),
                              paste0(36:25) ) )] %>% 
              select(-year) %>% 
              as.matrix()


# general dataset simulator ----------------------------------


# simulate data changing beta
sim_beta <- function(weight = T,   spat_n = 5,
                     y_sd   = 0.3, beta_x = 1.2){
  
  if( weight ){
  
    # 'true' weight function. Mean in 5th month, sd = 1.
    pdens  <- dnorm(1:12, 5, 1)
    w_v    <- c( (pdens / sum(pdens)) * (1/3),
                 (pdens / sum(pdens)) * (1/3),
                 (pdens / sum(pdens)) * (1/3) )
    
    x1 <- anom_a[,1:12]  %*% w_v[1:12]  %>% as.numeric
    x2 <- anom_a[,13:24] %*% w_v[13:24] %>% as.numeric
    x3 <- anom_a[,25:36] %*% w_v[25:36] %>% as.numeric
    
  }
  
  if( !weight ){
  
    # calculate yearly anomalies 
    x1 <- anom_a[,1:12]  %>% rowMeans
    x2 <- anom_a[,13:24] %>% rowMeans
    x3 <- anom_a[,25:36] %>% rowMeans
    
  }
  
  # function to produce normally distributed response
  prod_y <- function(x){ 
    set.seed(x)
    rnorm(nrow(anom_a), 
          mean= 0 + 
                x1*-beta_x +
                x2*beta_x +
                x3*0, 
           sd = y_sd)
  }
  
  y_vec <- lapply(1:spat_n, prod_y) %>% 
              Reduce(function(...) c(...), .)
    
  # output data
  data.frame( x1 = rep(x1, spat_n),
              x2 = rep(x2, spat_n),
              x3 = rep(x3, spat_n),
              y  = y_vec,
              stringsAsFactors = F )

}

# exploratory plots
plot(y ~ x1, data=sim_beta())
plot(y ~ x2, data=sim_beta())
plot(y ~ x3, data=sim_beta())


# set up model runs 
setup_model_runs <- function( weight,   
                              spat_n,
                              y_sd  , 
                              beta_x ){
  
  # simulate data
  y_sim  <- sim_beta( weight = weight,   
                      spat_n = spat_n,
                      y_sd   = y_sd, 
                      beta_x = beta_x  )$y
  
  # set up climate variables
  clim_mult <- list( anom_a ) %>% 
                 rep( spat_n ) %>% 
                 Reduce( function(...) rbind(...), .)
  
  # simulation parameters
  sim_pars <- list(
    warmup = 1000, 
    iter   = 4000, 
    thin   = 2, 
    chains = 3
  )

  # organize data into list to pass to stan
  dat_stan <- list(
    n_time  = nrow(clim_mult),
    n_lag   = ncol(clim_mult),
    y       = y_sim,
    clim    = t(clim_mult),
    clim1   = t(clim_mult)[1:12 ,],
    clim2   = t(clim_mult)[13:24,],
    clim3   = t(clim_mult)[25:36,],
    M       = 12,    # number of months in a year
    K       = ncol(clim_mult) / 12,
    expp_beta = 20
  )

  dat_stan
  
}



mod_gaus <- stan_model(file = paste0("code/stan/normal_gaus_nest.stan") )
mod_gaus <- stan_model(file = paste0("code/stan/normal_gaus.stan") )


init_t <- Sys.time()
# gaussian moving window
fit_gaus1 <- sampling(
  object =mod_gaus,
  data = dat_stan,
  pars = c('sens_mu', 'sens_sd', 'theta_y', 'alpha', 'beta', 'y_sd'),
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains,
  control = list(adapt_delta = 0.999)#,
  #control = list(adapt_delta = 0.999, stepsize = 0.001, max_treedepth = 20)
)
  
fit_gaus2 <- sampling(
  object =mod_gaus,
  data = dat_stan,
  pars = c('sens_mu', 'sens_sd', 'theta_y', 'alpha', 'beta', 'y_sd'),
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains,
  control = list(adapt_delta = 0.999)#,
  #control = list(adapt_delta = 0.999, stepsize = 0.001, max_treedepth = 20)
)

fit_gaus3 <- sampling(
  object =mod_gaus,
  data = dat_stan,
  pars = c('sens_mu', 'sens_sd', 'theta_y', 'alpha', 'beta', 'y_sd'),
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains,
  control = list(adapt_delta = 0.999)#,
  #control = list(adapt_delta = 0.999, stepsize = 0.001, max_treedepth = 20)
)
Sys.time() - init_t 


tiff('results/simulations/bimodal/beta_spp_rep.tiff',
     width=6.3, height=8, res=300,unit='in',
     compression='lzw')

par(mfcol = c(3,3), 
    mar   = c(3,3,2,0.1), 
    mgp   = c(1.5,0.5,0) )
# plot(y ~ x1, data=sim_beta(1.2), main='True rel.; year 1')
# abline(0,-1.2,lwd=2)
# plot(y ~ x2, data=sim_beta(1.2), main='True rel.; year 2')
# abline(0,1.2,lwd=2)
# plot(y ~ x3, data=sim_beta(1.2), main='True rel.; year 3')
# abline(0,0,lwd=2)

extract(fit_gaus1)$beta    %>% hist(main = 'beta')
extract(fit_gaus1)$sens_mu %>% hist(main = 'sens_mu')
extract(fit_gaus1)$theta_y %>% boxplot(main='yr weights')

extract(fit_gaus2)$beta    %>% hist(main = 'beta')
extract(fit_gaus2)$sens_mu %>% hist(main = 'sens_mu')
extract(fit_gaus2)$theta_y %>% boxplot(main='yr weights')

extract(fit_gaus3)$beta    %>% hist(main = 'beta')
extract(fit_gaus3)$sens_mu %>% hist(main = 'sens_mu')
extract(fit_gaus3)$theta_y %>% boxplot(main='yr weights')

dev.off()


pryr::object_size(fit_gaus3)



# 
# # 1. single spatial replicate --------------------------------
# 
# # simulate data 
# 
# # "true" climate effects
# betas <- seq(0.2,2.2,by=0.5)
# 
# # list to store models 
# mod_l <- list( '0.2'  = NULL,
#                '0.45' = NULL,
#                '0.7'  = NULL,
#                '1.2'  = NULL,
#                '1.7'  = NULL,
#                '2.2'  = NULL )
# 
# # climate "data"
# clim_m <- select(anom_a[1:20,], -year) %>% as.matrix
# 
# # simulate data changing beta
# sim_beta <- function(beta_x){
#   
#   # 'true' weight function. Mean in 5th month, sd = 1.
#   pdens  <- dnorm(1:12, 5, 1)
#   w_v    <- c( (pdens / sum(pdens)) * (1/3),
#                (pdens / sum(pdens)) * (1/3),
#                (pdens / sum(pdens)) * (1/3) )
#   
#   x1 <- clim_m[,1:12]  %*% w_v[1:12]
#   x2 <- clim_m[,13:24] %*% w_v[13:24]
#   x3 <- clim_m[,25:36] %*% w_v[25:36]
#   
#   # function to produce normally distributed response
#   set.seed(1776)
#   prod_y1 <- function(x) rnorm(nrow(clim_m), 
#                                 mean= 0 + 
#                                       x1*-beta_x +
#                                       x2*beta_x +
#                                       x3*0, 
#                                 sd = log_lambda_sd)
#   
#   # output data
#   data.frame( x1 = x1,
#               x2 = x2,
#               x3 = x3,
#               y = prod_y1(x),
#               stringsAsFactors = F )
# 
# }
# 
# # plot it
# plot(y ~ x1, data=sim_beta(0.5))
# plot(y ~ x2, data=sim_beta(0.5))
# plot(y ~ x3, data=sim_beta(0.5))
# 
# # simulate data changing beta
# sim_weight <- function(beta_x){
#   
#   # 'true' weight function. Mean in 5th month, sd = 1.
#   pdens1 <- dnorm(1:12, 5, 1)
#   pdens2 <- dnorm(1:12, 1, 1)
#   pdens3 <- dnorm(1:12, 5, 1)
#   w_v    <- c( (pdens1 / sum(pdens1)) * (1/3),
#                (pdens2 / sum(pdens2)) * (1/3),
#                (pdens3 / sum(pdens3)) * (1/3) )
#   
#   x      <- clim_m  %*% w_v
#   
#   # function to produce normally distributed response
#   set.seed(101)
#   prod_y <- function(x) rnorm(length(x), 
#                               mean = 0 + x*beta_x, 
#                               sd   = log_lambda_sd)
#   
#   # output data
#   data.frame( x = x,
#               y = prod_y(x),
#               stringsAsFactors = F )
# 
# }
# 
# # plot it
# plot(y ~ x, data=sim_weight(1.2))
# 
# 
# # simulate data changing beta
# sim_weight_beta <- function(beta_x){
#   
#   # 'true' weight function. Mean in 5th month, sd = 1.
#   pdens1 <- dnorm(1:12, 5, 1)
#   pdens2 <- dnorm(1:12, 1, 1)
#   pdens3 <- dnorm(1:12, 5, 1)
#   w_v    <- c( (pdens1 / sum(pdens1)) * (1/3),
#                (pdens2 / sum(pdens2)) * (1/3),
#                (pdens3 / sum(pdens3)) * (1/3) )
#   
#   x1 <- clim_m[,1:12]  %*% w_v[1:12]
#   x2 <- clim_m[,13:24] %*% w_v[13:24]
#   x3 <- clim_m[,25:36] %*% w_v[25:36]
#   
#   # function to produce normally distributed response
#   set.seed(101)
#   prod_y <- function(x) rnorm(length(x1), 
#                               mean= 0 + 
#                                     x1*beta_x +
#                                     x2*-beta_x +
#                                     x3*beta_x, 
#                               sd = log_lambda_sd)
#   
#   
#   # output data
#   data.frame( x1 = x1,
#               x2 = x2,
#               x3 = x3,
#               y = prod_y(x),
#               stringsAsFactors = F )
# 
# }
# 
# # plot it
# par(mfrow=c(2,2), mar=c(2,2,0.1,0.1))
# plot(y ~ x1, data=sim_weight_beta(1.2))
# plot(y ~ x1, data=sim_weight_beta(1.2))
# plot(y ~ x3, data=sim_weight_beta(1.2))
# 
# 
# 
# # simulate data changing beta
# sim_beta_m <- function(beta_x){
#   
#   # 'true' weight function. Mean in 5th month, sd = 1.
#   pdens  <- dnorm(1:12, 5, 1)
#   w_v    <- c( (pdens / sum(pdens)) * (1/3),
#                (pdens / sum(pdens)) * (1/3),
#                (pdens / sum(pdens)) * (1/3) )
#   
#   x1 <- clim_m[,1:12]  %*% w_v[1:12]
#   x2 <- clim_m[,13:24] %*% w_v[13:24]
#   x3 <- clim_m[,25:36] %*% w_v[25:36]
#   
#   # function to produce normally distributed response
#   set.seed(1776)
#   prod_y <- function(x) rnorm(nrow(clim_m), 
#                               mean= 0 + 
#                                     x1*-beta_x +
#                                     x2*beta_x +
#                                     x3*beta_x, 
#                               sd = log_lambda_sd)
#   
#   # output data
#   data.frame( x1 = x1,
#               x2 = x2,
#               x3 = x3,
#               y = prod_y(x),
#               stringsAsFactors = F )
# 
# }
# 
# 
# # fit data ---------------------------------------------------
# 
# # set rstan options
# rstan_options( auto_write = TRUE )
# options( mc.cores = parallel::detectCores() )
# 
# # simulation parameters
# sim_pars <- list(
#   warmup = 1000, 
#   iter   = 4000, 
#   thin   = 2, 
#   chains = 3
# )
# 
# # here we focus on a normally distributed response
# family <- 'normal'
# 
# # compile the models
# mod_gaus <- stan_model(file = paste0("code/stan/",family,"_gaus_nest.stan") )
# mod_expp <- stan_model(file = paste0("code/stan/",family,"_expp_nest.stan") )
# mod_gev  <- stan_model(file = paste0("code/stan/",family,"_gev_nest.stan")  )
# mod_sad  <- stan_model(file = paste0("code/stan/",family,"_dirichlet_nest.stan") )
# 
# 
# # fit models based on data
# fit_mods <- function(df_x){
#   
#   # organize data into list to pass to stan
#   dat_stan <- list(
#     n_time  = nrow(clim_m),
#     n_lag   = ncol(clim_m),
#     y       = df_x$y,
#     clim    = t(clim_m),
#     clim1   = t(clim_m)[1:12 ,],
#     clim2   = t(clim_m)[13:24,],
#     clim3   = t(clim_m)[25:36,],
#     M       = 12,    # number of months in a year
#     K       = ncol(clim_m) / 12,
#     expp_beta = 20
#   )
#   
#   # gaussian moving window
#   fit_gaus <- sampling(
#     object =mod_gaus,
#     data = dat_stan,
#     pars = c('sens_mu', 'sens_sd', 'theta_y', 'alpha', 'beta', 'y_sd'),
#     warmup = sim_pars$warmup,
#     iter = sim_pars$iter,
#     thin = sim_pars$thin,
#     chains = sim_pars$chains,
#     control = list(adapt_delta = 0.999)#,
#     #control = list(adapt_delta = 0.999, stepsize = 0.001, max_treedepth = 20)
#   )
#   
#   # exponential power moving window
#   fit_expp <- sampling(
#     object = mod_expp,
#     data = dat_stan,
#     pars = c('sens_mu', 'sens_sd', 'theta_y', 'alpha', 'beta', 'y_sd'),
#     warmup = sim_pars$warmup,
#     iter = sim_pars$iter,
#     thin = sim_pars$thin,
#     chains = sim_pars$chains,
#     control = list(adapt_delta = 0.999)#,
#     #control = list(adapt_delta = 0.999, stepsize = 0.001, max_treedepth = 20)
#   )
#   
#   # Generalized Extreme Value nested
#   fit_gev <- sampling(
#     object = mod_gev,
#     data = dat_stan,
#     pars = c('loc', 'scale', "shape", 'theta_y', 'alpha', 'beta', 'y_sd'),
#     warmup = sim_pars$warmup,
#     iter = sim_pars$iter,
#     thin = sim_pars$thin,
#     chains = sim_pars$chains#,
#     #control = list(adapt_delta = 0.999, stepsize = 0.001, max_treedepth = 20)
#   )
#   
#   # Simplex nested
#   fit_sad <- sampling(
#     object = mod_sad,
#     data = dat_stan,
#     pars = c('theta_y', 'theta_m', 'alpha', 'beta', 'y_sd'),
#     warmup = sim_pars$warmup,
#     iter = sim_pars$iter,
#     thin = sim_pars$thin,
#     chains = sim_pars$chains,
#     control = list(adapt_delta = 0.999)
#   )
# 
#   list( fit_gaus, 
#         fit_expp,
#         fit_gev,
#         fit_sad )
#   
# }
# 
# model_names <- function(x){
#   x %>% setNames( c('gaus','expp','gev','sad') )
# }
# 
# # set up lists containing models
# mod_b_l  <- 
# mod_w_l  <- 
# mod_wb_l <- list( '1.2' = NULL )
# 
# #fit 
# mod_b_l$'1.2'  <- fit_mods( sim_beta(1.2) )        %>% model_names
# mod_w_l$'1.2'  <- fit_mods( sim_weight(1.2) )      %>% model_names
# mod_wb_l$'1.2' <- fit_mods( sim_weight_beta(1.2) ) %>% model_names
# # mod_b_l$'0.2'  <- fit_mods( sim_beta(0.2) ) %>% model_names
# # mod_b_l$'2.2'  <- fit_mods( sim_beta(2.2) ) %>% model_names
# 
# # posteriors: beta
# extract(mod_b_l$'1.2'[[1]])$beta  %>% hist
# extract(mod_w_l$'1.2'[[1]])$beta  %>% hist
# extract(mod_wb_l$'1.2'[[1]])$beta %>% hist
# 
# # sens mu 
# extract(mod_b_l$'1.2'[[2]])$sens_mu %>% hist
# extract(mod_w_l$'1.2'[[2]])$sens_mu %>% hist
# extract(mod_wb_l$'1.2'[[1]])$sens_mu %>% hist
# 
# # thetas 
# extract(mod_b_l$'1.2'[[1]])$theta_y  %>% boxplot
# extract(mod_w_l$'1.2'[[1]])$theta_y  %>% boxplot
# extract(mod_wb_l$'1.2'[[1]])$theta_y %>% boxplot
# 
# 
# tiff('results/simulations/bimodal/beta.tiff',
#      width=6.3, height=6.3, res=300,unit='in',
#      compression='lzw')
# 
# par(mfrow = c(2,2), 
#     mar   = c(3,3,2,0.1), 
#     mgp   = c(1.5,0.5,0) )
# extract(mod_b_l$'1.2'[[1]])$beta    %>% hist(main = 'beta')
# extract(mod_b_l$'1.2'[[1]])$sens_mu %>% hist(main = 'sens_mu')
# extract(mod_b_l$'1.2'[[1]])$theta_y %>% boxplot(main='yr weights')
# 
# dev.off()
# 
# 
# tiff('results/simulations/bimodal/weight.tiff',
#      width=6.3, height=6.3, res=300,unit='in',
#      compression='lzw')
# 
# par(mfrow = c(2,2), 
#     mar   = c(3,3,2,0.1), 
#     mgp   = c(1.5,0.5,0) )
# extract(mod_w_l$'1.2'[[1]])$beta    %>% hist(main = 'beta')
# extract(mod_w_l$'1.2'[[1]])$sens_mu %>% hist(main = 'sens_mu')
# extract(mod_w_l$'1.2'[[1]])$theta_y %>% boxplot(main='yr weights')
# 
# dev.off()
# 
# 
# 
# tiff('results/simulations/bimodal/beta_weight.tiff',
#      width=6.3, height=6.3, res=300,unit='in',
#      compression='lzw')
# 
# par(mfrow = c(2,2), 
#     mar   = c(3,3,2,0.1), 
#     mgp   = c(1.5,0.5,0) )
# extract(mod_wb_l$'1.2'[[1]])$beta    %>% hist(main = 'beta')
# extract(mod_wb_l$'1.2'[[1]])$sens_mu %>% hist(main = 'sens_mu')
# extract(mod_wb_l$'1.2'[[1]])$theta_y %>% boxplot(main='yr weights')
# 
# dev.off()
# 
# 
# 
# # 2. 5 spatial replicates --------------------------
# 
# # "true" climate effects
# betas <- seq(0.2,2.2,by=0.5)
# 
# # list to store models 
# mod_l <- list( '1.2'  = NULL )
# 
# # climate "data"
# clim_m <- select(anom_a[1:20,], -year) %>% as.matrix
# 
# # simulate data changing beta
# sim_beta <- function(beta_x){
#   
#   # 'true' weight function. Mean in 5th month, sd = 1.
#   pdens  <- dnorm(1:12, 5, 1)
#   w_v    <- c( (pdens / sum(pdens)) * (1/3),
#                (pdens / sum(pdens)) * (1/3),
#                (pdens / sum(pdens)) * (1/3) )
#   
#   x1 <- clim_m[,1:12]  %*% w_v[1:12]
#   x2 <- clim_m[,13:24] %*% w_v[13:24]
#   x3 <- clim_m[,25:36] %*% w_v[25:36]
#   
#   # function to produce normally distributed response
#   set.seed(1)
#   prod_y1 <- function(x) rnorm(nrow(clim_m), 
#                                 mean= 0 + 
#                                       x1*-beta_x +
#                                       x2*beta_x +
#                                       x3*0, 
#                                 sd = log_lambda_sd)
#   set.seed(2)
#   prod_y2 <- function(x) rnorm(nrow(clim_m), 
#                                 mean= 0 + 
#                                       x1*-beta_x +
#                                       x2*beta_x +
#                                       x3*0, 
#                                 sd = log_lambda_sd)
#   set.seed(3)
#   prod_y3 <- function(x) rnorm(nrow(clim_m), 
#                                 mean= 0 + 
#                                       x1*-beta_x +
#                                       x2*beta_x +
#                                       x3*0, 
#                                 sd = log_lambda_sd)
#   set.seed(4)
#   prod_y4 <- function(x) rnorm(nrow(clim_m), 
#                                 mean= 0 + 
#                                       x1*-beta_x +
#                                       x2*beta_x +
#                                       x3*0, 
#                                 sd = log_lambda_sd)
#   set.seed(5)
#   prod_y5 <- function(x) rnorm(nrow(clim_m), 
#                                 mean= 0 + 
#                                       x1*-beta_x +
#                                       x2*beta_x +
#                                       x3*0, 
#                                 sd = log_lambda_sd)
#   
#   # output data
#   data.frame( x1 = rep(x1,5),
#               x2 = rep(x2,5),
#               x3 = rep(x3,5),
#               y  = c(prod_y1(), prod_y2(), prod_y3(),
#                      prod_y4(), prod_y5()),
#               stringsAsFactors = F )
# 
# }
# 
# plot(y ~ x1, data=sim_beta(0.5))
# plot(y ~ x2, data=sim_beta(0.5))
# plot(y ~ x3, data=sim_beta(0.5))
# 
# 
# clim_mult = Reduce( function(...) rbind(...), 
#                     list(clim_m, clim_m, clim_m,
#                          clim_m, clim_m) )
# 
# # organize data into list to pass to stan
# dat_stan <- list(
#   n_time  = nrow(clim_mult),
#   n_lag   = ncol(clim_mult),
#   y       = sim_beta(1.2)$y,
#   clim    = t(clim_mult),
#   clim1   = t(clim_mult)[1:12 ,],
#   clim2   = t(clim_mult)[13:24,],
#   clim3   = t(clim_mult)[25:36,],
#   M       = 12,    # number of months in a year
#   K       = ncol(clim_mult) / 12,
#   expp_beta = 20
# )
# 
# # gaussian moving window
# fit_gaus <- sampling(
#   object =mod_gaus,
#   data = dat_stan,
#   pars = c('sens_mu', 'sens_sd', 'theta_y', 'alpha', 'beta', 'y_sd'),
#   warmup = sim_pars$warmup,
#   iter = sim_pars$iter,
#   thin = sim_pars$thin,
#   chains = sim_pars$chains,
#   control = list(adapt_delta = 0.999)#,
#   #control = list(adapt_delta = 0.999, stepsize = 0.001, max_treedepth = 20)
# )
#   
# 
# 
# tiff('results/simulations/bimodal/beta_spp_rep.tiff',
#      width=6.3, height=8, res=300,unit='in',
#      compression='lzw')
# 
# par(mfrow = c(3,2), 
#     mar   = c(3,3,2,0.1), 
#     mgp   = c(1.5,0.5,0) )
# plot(y ~ x1, data=sim_beta(1.2), main='True rel.; year 1')
# abline(0,-1.2,lwd=2)
# plot(y ~ x2, data=sim_beta(1.2), main='True rel.; year 2')
# abline(0,1.2,lwd=2)
# plot(y ~ x3, data=sim_beta(1.2), main='True rel.; year 3')
# abline(0,0,lwd=2)
# 
# extract(fit_gaus)$beta    %>% hist(main = 'beta')
# extract(fit_gaus)$sens_mu %>% hist(main = 'sens_mu')
# extract(fit_gaus)$theta_y %>% boxplot(main='yr weights')
# 
# dev.off()