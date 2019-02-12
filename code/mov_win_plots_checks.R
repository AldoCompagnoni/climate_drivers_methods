# Script to check:
# data ~ x_ante
# posterior predictive checks
source("C:/CODE/moving_windows/format_data.R")
library(tidyverse)
library(dismo)
library(mgcv)
library(testthat)
library(rstan)
library(evd)

# set rstan options
rstan_options( auto_write = TRUE )
options( mc.cores = parallel::detectCores() )

# "pipeable" Reduce rbind
rbind_l <- function(x) Reduce(function(...) rbind(...), x)

# climate predictor, response, months back, max. number of knots
response  <- "log_lambda"
clim_var  <- "precip"
m_back    <- 36    
st_dev    <- FALSE

# read data -----------------------------------------------------------------------------------------
lam       <- read.csv("C:/cloud/Dropbox/sApropos/all_demog_6tr.csv", stringsAsFactors = F)
m_info    <- read.csv("C:/cloud/Dropbox/sApropos/MatrixEndMonth_information.csv", stringsAsFactors = F)
clim      <- data.table::fread(paste0('C:/cloud/Dropbox/sApropos/',clim_var,"_chelsa_hays.csv"),  stringsAsFactors = F)
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

for(ii in 13:35){
  
  # set species (I pick Sphaeraclea_coccinea)
  # ii            <- 34
  spp_name      <- spp[ii]
  
  if( spp_name == "Brassica_insularis"){
    lam <- lam %>% 
              subset( !(SpeciesAuthor == "Brassica_insularis" & 
                        lambda < 0.6) ) 
  
  }
  
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
  
  
  # year t
  dat_stan$clim_means  <- rowMeans(mod_data$climate[,1:12 ])
  fit_yr1 <- stan(
    file = paste0("code/stan/",family,"_yr.stan"),
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
    file = paste0("code/stan/",family,"_yr.stan"),
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
    file = paste0("code/stan/",family,"_yr.stan"),
    data = dat_stan,
    pars = c('alpha', 'beta', 'y_sd', 'log_lik'),
    warmup = sim_pars$warmup,
    iter = sim_pars$iter,
    thin = sim_pars$thin,
    chains = sim_pars$chains#,
    #control = list(adapt_delta = 0.999, stepsize = 0.001, max_treedepth = 20)
  )
  dat_stan$clim_means  <- rowMeans(mod_data$climate)
  
  
  # gaussian moving window
  fit_gaus <- stan(
    file = paste0("code/stan/",family,"_gaus.stan"),
    data = dat_stan,
    pars = c('sens_mu', 'sens_sd', 'alpha', 'beta', 'y_sd'), #, 'log_lik'
    warmup = sim_pars$warmup,
    iter = sim_pars$iter,
    thin = sim_pars$thin,
    chains = sim_pars$chains,
    control = list(adapt_delta = 0.999, stepsize = 0.001, max_treedepth = 20)
  )
  
  # exponential power moving window
  fit_expp <- stan(
    file = paste0("code/stan/",family,"_expp.stan"),
    data = dat_stan,
    pars = c('sens_mu', 'sens_sd', 'alpha', 'beta', 'y_sd'), #'log_lik'
    warmup = sim_pars$warmup,
    iter = sim_pars$iter,
    thin = sim_pars$thin,
    chains = sim_pars$chains,
    control = list(adapt_delta = 0.999, stepsize = 0.001, max_treedepth = 20)
  )
  
  # Generalized extreme value 
  fit_gev <- stan(
    file = paste0("code/stan/",family,"_gev_scraped.stan"),
    data = dat_stan,
    pars = c('mu', 'sigma', 'xi', 'alpha', 'beta', 'y_sd'), #, "shape", 'log_lik'
    # pars = c('loc', 'scale', 'alpha', 'beta', 'y_sd'), #, "shape", 'log_lik'
    warmup = sim_pars$warmup,
    iter = sim_pars$iter,
    thin = sim_pars$thin,
    init_r = 2,
    chains = sim_pars$chains#,
    # control = list(adapt_delta = 0.999, stepsize = 0.001, max_treedepth = 20)
  )
  
  
  # # Simplex - 36 months
  # dat_stan$clim <- t(mod_data$climate)
  # fit_d <- stan(
  #   file = paste0("code/stan/",family,"_dirichlet.stan"),
  #   data = dat_stan,
  #   pars = c('theta', 'alpha', 'beta', 'y_sd', 'log_lik'),
  #   warmup = sim_pars$warmup,
  #   iter = sim_pars$iter,
  #   thin = sim_pars$thin,
  #   chains = sim_pars$chains#,
  #   #control = list(adapt_delta = 0.999, stepsize = 0.001, max_treedepth = 20)
  # )
  # dat_stan$clim <- mod_data$climate
  #  
  
  # Nested models 
  # update data list
  dat_stan$clim         <- t(mod_data$climate)
  dat_stan$clim1        <- t(mod_data$climate)[1:12 ,]
  dat_stan$clim2        <- t(mod_data$climate)[13:24,]
  dat_stan$clim3        <- t(mod_data$climate)[25:36,]
  
  # Simplex nested
  fit_36_nest <- stan(
    file = paste0("code/stan/",family,"_dirichlet_nest.stan"),
    data = dat_stan,
    pars = c('theta_y', 'theta_m', 'alpha', 'beta', 'y_sd'), #, 'log_lik'
    warmup = sim_pars$warmup,
    iter = sim_pars$iter,
    thin = sim_pars$thin,
    chains = sim_pars$chains#,
    #control = list(adapt_delta = 0.999, stepsize = 0.001, max_treedepth = 20)
  )
  
  # Generalized Extreme Value nested
  fit_gev_nest <- stan(
    file = paste0("code/stan/",family,"_gev_nest.stan"),
    data = dat_stan,
    pars = c('loc', 'scale', "shape", 'theta_y', 'alpha', 'beta', 'y_sd'),# 'log_lik'
    warmup = sim_pars$warmup,
    iter = sim_pars$iter,
    thin = sim_pars$thin,
    chains = sim_pars$chains#,
    #control = list(adapt_delta = 0.999, stepsize = 0.001, max_treedepth = 20)
  )
  
  # Power exponential nested 
  fit_expp_nest <- stan(
    file = paste0("code/stan/",family,"_expp_nest.stan"),
    data = dat_stan,
    pars = c('sens_mu', 'sens_sd', 'theta_y', 'alpha', 'beta', 'y_sd'), #,'log_lik'
    warmup = sim_pars$warmup,
    iter = sim_pars$iter,
    thin = sim_pars$thin,
    chains = sim_pars$chains,
    control = list(adapt_delta = 0.999, stepsize = 0.001, max_treedepth = 20)
  )
  
  
  # plots ------------------------------------------------------------
  
  # store all models 
  all_fits <- list( fit_yr1       = fit_yr1,
                    fit_yr2       = fit_yr2,
                    fit_yr3       = fit_yr3,
                    fit_gaus      = fit_gaus,
                    fit_expp      = fit_expp,
                    fit_gev       = fit_gev,
                    fit_36_nest   = fit_36_nest,
                    fit_gev_nest  = fit_gev_nest,
                    fit_expp_nest = fit_expp_nest )
  
  
  # exponential power distribution
  dexppow <- function(x, mu, sigma, beta) {
      return((beta / (2 * sigma * gamma (1.0/beta)) ) *
               exp(-(abs(x - mu)/sigma)^beta));
  }
  
  # plot it!
  plot_spp <- function(fit_obj,mod){
    
    # parameters
    response  <- "log_lambda"
    family    <- "normal"
    expp_beta <- 20
    m_back    <- 36    
    xx        <- 1:12
    xx_a      <- 1:36
    
    # extract means from posterior
    get_means <- function(fit_obj){
      fit_obj %>%
        rstan::extract() %>% 
        as.data.frame %>% 
        as.matrix %>% 
        colMeans %>% 
        # allow subsetting atomic vectors
        as.data.frame %>% t %>% as.data.frame
    }
    
    # produce data frame to plot data
    plotting_df <- function(x_ante,post_mean,mod_data){
      data.frame( y = mod_data$resp$log_lambda,
                  x = x_ante,
                  stringsAsFactors = F) %>% 
        mutate( alpha = post_mean$alpha,
                beta  = post_mean$beta )
    }
    
    
    # plot model over data --------------------------------------------------------------
    if( grepl('fit_yr', mod) ){
      
      # mean params
      post_mean <- get_means(fit_obj)
      
      # gaus
      if(mod=='fit_yr1'){
        clim_sub <- as.matrix(mod_data$climate[,c(1:12)]) 
      }
      if(mod=='fit_yr2'){
        clim_sub <- as.matrix(mod_data$climate[,c(13:24)])
      }
      if(mod=='fit_yr3'){
        clim_sub <- as.matrix(mod_data$climate[,c(25:36)]) 
      }
      
      # plotting material & plot data frame
      x_ante <- rowMeans(clim_sub)
      pl_df  <- plotting_df(x_ante,post_mean,mod_data)
      
    }
    
  
    if(mod == 'fit_gaus'){
      
      # mean params
      post_mean <- get_means(fit_obj)
      
      # gaus
      gaus_w_v  <- dnorm(xx_a, post_mean$sens_mu, post_mean$sens_sd )
      w_v      <- gaus_w_v / sum(gaus_w_v)
      
      # plotting material & plot data frame
      x_ante <- as.matrix(mod_data$climate) %*% w_v
      pl_df  <- plotting_df(x_ante,post_mean,mod_data)
      
    }
    
    if(mod == 'fit_expp'){
      
      # mean params
      post_mean <- get_means(fit_obj)
      
      # gaus
      expp_w_v  <- dexppow(xx_a, post_mean$sens_mu, 
                                 post_mean$sens_sd,
                                 expp_beta)
      w_v       <- expp_w_v / sum(expp_w_v)
      
      # plotting material & plot data frame
      x_ante <- as.matrix(mod_data$climate) %*% w_v
      pl_df  <- plotting_df(x_ante,post_mean,mod_data)
      
    }
    
    # gev
    if(mod == 'fit_gev'){
      
      # let's start from SAD
      post_mean <- get_means(fit_obj)
      
      # gev
      gev_w_v  <- dgev(xx_a, post_mean$mu, post_mean$sigma,
                             post_mean$xi )
      w_v      <- gev_w_v / sum(gev_w_v)
                     
      # plotting material & plot data frame
      x_ante <- as.matrix(mod_data$climate) %*% w_v
      pl_df  <- plotting_df(x_ante,post_mean,mod_data)
      
    }
    
    if(mod == 'fit_expp_nest'){
    
      # let's start from SAD
      post_mean <- get_means(fit_obj)
      
      # expp
      expp_w_v <- dexppow(xx, post_mean$sens_mu, 
                              post_mean$sens_sd, 20)
      w_v <- c( ((expp_w_v / sum(expp_w_v))*post_mean$theta_y.1),
                ((expp_w_v / sum(expp_w_v))*post_mean$theta_y.2),
                ((expp_w_v / sum(expp_w_v))*post_mean$theta_y.3) )
     
      # plotting material & plot data frame
      x_ante <- as.matrix(mod_data$climate) %*% w_v
      pl_df  <- plotting_df(x_ante,post_mean,mod_data)
          
    }
    
    if(mod == 'fit_gev_nest'){
    
      # let's start from SAD
      post_mean <- get_means(fit_obj)
      
      # expp
      expp_w_v <- dgev(xx, post_mean$loc, post_mean$scale, 
                           post_mean$shape)
      w_v <- c( ((expp_w_v / sum(expp_w_v))*post_mean$theta_y.1),
                ((expp_w_v / sum(expp_w_v))*post_mean$theta_y.2),
                ((expp_w_v / sum(expp_w_v))*post_mean$theta_y.3) )
     
      # plotting material & plot data frame
      x_ante <- as.matrix(mod_data$climate) %*% w_v
      pl_df  <- plotting_df(x_ante,post_mean,mod_data)
          
    }
    
    if(mod == 'fit_36_nest'){
      
      post_mean <- get_means(fit_obj)
      
      # x antecedent
      w_v    <- c(as.numeric(post_mean[paste0('theta_m.',1:12)]*post_mean$theta_y.1),
                  as.numeric(post_mean[paste0('theta_m.',1:12)]*post_mean$theta_y.2),
                  as.numeric(post_mean[paste0('theta_m.',1:12)]*post_mean$theta_y.3) )
      
      # plotting material & plot data frame
      x_ante <- as.matrix(mod_data$climate) %*% w_v
      pl_df  <- plotting_df(x_ante,post_mean,mod_data)
      
    }
    
    # plot it out
    clim_var_print <- gsub('ip','',clim_var)
    
    # plot it out
    ggplot(pl_df, aes(x, y) ) +
      geom_point() + 
      geom_abline( intercept = unique(pl_df$alpha), 
                   slope     = unique(pl_df$beta) ) + 
      ylab( expression('log('*lambda*')') ) +
      xlab( 'X antecedent' ) + 
      ggtitle( spp_name ) +
      ggsave( paste0('results/checks/',
                     clim_var_print,'/',
                     spp_name,'_',
                     mod,'.tiff'),
              width = 6.3, height = 5,
              compression="lzw")
  
  }
  
  Map(plot_spp, all_fits, names(all_fits))
  print(ii)
  
  
  # Bayesian checks -----------------------------------------
  
  # sampled posterior predictive checks
  sppc <- function(fit_obj, mod){
    
    # parameters
    response  <- "log_lambda"
    family    <- "normal"
    expp_beta <- 20
    m_back    <- 36    
    xx        <- 1:12
    xx_a      <- 1:36
    
    # extract ONE sample from posterior
    get_sample <- function(fit_obj){
      fit_obj %>%
        rstan::extract() %>% 
        as.data.frame %>% 
        .[100,]
    }
    
    # produce data frame to plot data
    sppc_df <- function(x_ante,post_mean,mod_data){
      data.frame( y = mod_data$resp$log_lambda,
                  x = x_ante,
                  stringsAsFactors = F) %>% 
        mutate( alpha = post_mean$alpha,
                beta  = post_mean$beta,
                y_sd  = post_mean$y_sd )
    }
    
    
    # plot model over data --------------------------------------------------------------
    if( grepl('fit_yr', mod) ){
      
      # mean params
      post_smpl <- get_sample(fit_obj)
      
      # gaus
      if(mod=='fit_yr1'){
        clim_sub <- as.matrix(mod_data$climate[,c(1:12)]) 
      }
      if(mod=='fit_yr2'){
        clim_sub <- as.matrix(mod_data$climate[,c(13:24)])
      }
      if(mod=='fit_yr3'){
        clim_sub <- as.matrix(mod_data$climate[,c(25:36)]) 
      }
      
      # plotting material & plot data frame
      x_ante <- rowMeans(clim_sub)
      pl_df  <- sppc_df(x_ante,post_smpl,mod_data)
      
    }
    
  
    if(mod == 'fit_gaus'){
      
      # mean params
      post_smpl <- get_sample(fit_obj)
      
      # gaus
      gaus_w_v  <- dnorm(xx_a, post_smpl$sens_mu, post_smpl$sens_sd )
      w_v       <- gaus_w_v / sum(gaus_w_v)
      
      # plotting material & plot data frame
      x_ante <- as.matrix(mod_data$climate) %*% w_v
      pl_df  <- sppc_df(x_ante,post_smpl,mod_data)
      
    }
    
    if(mod == 'fit_expp'){
      
      # mean params
      post_smpl <- get_sample(fit_obj)
      
      # gaus
      expp_w_v  <- dexppow(xx_a, post_smpl$sens_mu, 
                                 post_smpl$sens_sd,
                                 expp_beta)
      w_v       <- expp_w_v / sum(expp_w_v)
      
      # plotting material & plot data frame
      x_ante <- as.matrix(mod_data$climate) %*% w_v
      pl_df  <- sppc_df(x_ante,post_smpl,mod_data)
      
    }
    
    if(mod == 'fit_gev'){
      
      # let's start from SAD
      post_smpl <- get_sample(fit_obj)
      
      # gev
      gev_w_v  <- dgev(xx_a, post_smpl$mu, post_smpl$sigma,
                             post_smpl$xi )
      w_v      <- gev_w_v / sum(gev_w_v)
                     
      # plotting material & plot data frame
      x_ante <- as.matrix(mod_data$climate) %*% w_v
      pl_df  <- sppc_df(x_ante,post_smpl,mod_data)
      
    }
    
    if(mod == 'fit_expp_nest'){
    
      # let's start from SAD
      post_smpl <- get_sample(fit_obj)
      
      # expp
      expp_w_v <- dexppow(xx, post_smpl$sens_mu, 
                              post_smpl$sens_sd, 20)
      w_v <- c( ((expp_w_v / sum(expp_w_v))*post_smpl$theta_y.1),
                ((expp_w_v / sum(expp_w_v))*post_smpl$theta_y.2),
                ((expp_w_v / sum(expp_w_v))*post_smpl$theta_y.3) )
     
      # plotting material & plot data frame
      x_ante <- as.matrix(mod_data$climate) %*% w_v
      pl_df  <- sppc_df(x_ante,post_smpl,mod_data)
          
    }
    
    if(mod == 'fit_gev_nest'){
    
      # let's start from SAD
      post_smpl <- get_sample(fit_obj)
      
      # expp
      expp_w_v <- dgev(xx, post_smpl$loc, post_smpl$scale, 
                           post_smpl$shape)
      w_v <- c( ((expp_w_v / sum(expp_w_v))*post_smpl$theta_y.1),
                ((expp_w_v / sum(expp_w_v))*post_smpl$theta_y.2),
                ((expp_w_v / sum(expp_w_v))*post_smpl$theta_y.3) )
     
      # plotting material & plot data frame
      x_ante <- as.matrix(mod_data$climate) %*% w_v
      pl_df  <- sppc_df(x_ante,post_smpl,mod_data)
          
    }
    
    if(mod == 'fit_36_nest'){
      
      post_smpl <- get_sample(fit_obj)
      
      # x antecedent
      w_v    <- c(as.numeric(post_smpl[paste0('theta_m.',1:12)]*post_smpl$theta_y.1),
                  as.numeric(post_smpl[paste0('theta_m.',1:12)]*post_smpl$theta_y.2),
                  as.numeric(post_smpl[paste0('theta_m.',1:12)]*post_smpl$theta_y.3) )
      
      # plotting material & plot data frame
      x_ante <- as.matrix(mod_data$climate) %*% w_v
      pl_df  <- sppc_df(x_ante,post_smpl,mod_data)
      
    }
    
    # plot it out
    clim_var_print <- gsub('ip','',clim_var)
    
      
    # Bayesian p-value calculation: sample posterior
    sim_post <- function(ss){
      dis_df <- pl_df %>% 
                  mutate( y_pred = alpha + beta*x ) %>%
                  mutate( y_rep  = rnorm(y_pred, 
                                         mean=y_pred, 
                                         sd=y_sd) ) %>% 
                  # calculate squared residuals
                  mutate( sr_dat = (y     - y_pred)^2,
                          sr_sim = (y_rep - y_pred)^2 )
      
      # Discrepancy (sum of squared residuals)
      discr <- dis_df %>% 
                  select(sr_dat, sr_sim) %>% 
                  as.matrix %>%
                  colSums
       
      if( discr['sr_sim'] > discr['sr_dat'] ) 1 else 0
              
    }
    
    res   <- sapply(1:1000, sim_post)
    p_val <- sum(res) / 1000
  
    data.frame( clim_var  = clim_var_print,
                mod       = mod,
                spp_name  = spp_name,
                p_val     = p_val,
                stringsAsFactors = F)
  
  }
  
  sppc_l <- Map(sppc, all_fits, names(all_fits))
  
  # store results
  write.csv(sppc_l %>% bind_rows,
            paste0('results/checks/bayes_pvalues/',
                   spp_name,'_',clim_var_print,
                   '.csv'),
            row.names=F)

}
  
# # plot dirichlet
# fit_36_nest %>% 
#   rstan::extract() %>% 
#   as.data.frame %>%
#   select( paste0('theta_y.',1:3) ) %>% 
#   gather(theta,value,theta_y.1:theta_y.3) %>% 
#   mutate( month = theta ) %>% 
#   mutate( month = gsub('theta_y.','',month) ) %>% 
#   mutate( month = factor(month,
#                          levels = paste0(1:3)) ) %>% 
#   ggplot( aes(x=month,y=value)) +
#   geom_boxplot() +
#   theme( axis.text = element_text(angle=90) )
# 
# fit_36_nest %>% 
#   rstan::extract() %>% 
#   as.data.frame %>%
#   select( paste0('theta_m.',1:12) ) %>% 
#   gather(theta,value,theta_m.1:theta_m.12) %>% 
#   mutate( month = theta ) %>% 
#   mutate( month = gsub('theta_m.','',month) ) %>% 
#   mutate( month = factor(month,
#                          levels = paste0(1:12)) ) %>% 
#   ggplot( aes(x=month,y=value)) +
#   geom_boxplot() +
#   theme( axis.text = element_text(angle=90) )
