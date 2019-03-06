o# Read data
# Read model results
# plot spp. specific results
rm(list=ls())
options(stringsAsFactors=F)
source("C:/CODE/moving_windows/format_data.R")
library(tidyverse)
library(dismo)
library(mgcv)
library(testthat)
library(rstan)
library(evd)

# climate predictor, response, months back, max. number of knots
response  <- "log_lambda"
clim_var  <- "airt"
m_back    <- 36    
st_dev    <- FALSE

# read data -----------------------------------------------------------------------------------------

# model data
lam       <- read.csv("data/all_demog_6tr.csv", stringsAsFactors = F)
clim      <- data.table::fread(paste0('data/',clim_var,"_chelsa_hays.csv"),  stringsAsFactors = F)
spp       <- lam$SpeciesAuthor %>% unique

# model data to update species list
mod_dir   <- 'I:/sie/101_data_AC/'
sum_f_l   <- grep('posterior',
                  list.files(paste0(mod_dir,'results/',clim_var)),
                  value=T)
spp_names <- gsub('posterior_|.csv','',sum_f_l)
spp       <- intersect(spp,spp_names)


# format data --------------------------------------------------------------------------------------

# set up model "family" based on response
if( response == "surv" | response == "grow" )             family = "beta" 
if( response == "fec" )                                   family = "gamma"
if( grepl("PreRep", response) | grepl("Rep", response) )  family = "beta"
if( response == "rho" | response == "react_fsa" )         family = "gamma"
if( response == "log_lambda")                             family = "normal"

expp_beta     <- 20

for(ii in 1:33){
  
  # set species (I pick Sphaeraclea_coccinea)
  # ii            <- 33
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
  
  # model posterior
  post_all <- data.table::fread(paste0(mod_dir,'results/',clim_var,'/',
                    grep(spp_name,sum_f_l,value=T))) %>% 
              # get the means of each parameter
              group_by( model ) %>% 
              summarise_all( mean ) %>% 
              ungroup
  
  
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
                    rowMeans(mod_data$climate[,25:36]) ) %>% 
                do.call(rbind, .),
    M       = 12,    # number of months in a year
    K       = ncol(mod_data$climate) / 12,
    S       = mod_data$resp$population %>% unique %>% length,
    site_i  = mod_data$resp$population %>% as.factor %>% as.numeric,
    expp_beta = expp_beta
  )
  
  # exponential power distribution
  dexppow <- function(x, mu, sigma, beta) {
      return((beta / (2 * sigma * gamma (1.0/beta)) ) *
               exp(-(abs(x - mu)/sigma)^beta));
  }
  
  # plot it!
  plot_spp <- function(mod){
    
    # parameters
    response  <- "log_lambda"
    family    <- "normal"
    expp_beta <- 20
    m_back    <- 36    
    xx        <- 1:12
    xx_a      <- 1:36
    
    # produce data frame to plot data
    plotting_df <- function(x_ante,post_mean,mod_data){
      data.frame( y = mod_data$resp$log_lambda,
                  x = x_ante,
                  stringsAsFactors = F) %>% 
        mutate( alpha = post_mean$alpha,
                beta  = post_mean$beta )
    }
    
    
    # plot model over data --------------------------------------------------------------
    if( grepl('yr', mod) ){
      
      # mean params
      post_mean <- post_all %>% subset(model == mod )
      
      # gaus
      if(mod=='yr1'){
        clim_sub <- as.matrix(mod_data$climate[,c(1:12)]) 
      }
      if(mod=='yr2'){
        clim_sub <- as.matrix(mod_data$climate[,c(13:24)])
      }
      if(mod=='yr3'){
        clim_sub <- as.matrix(mod_data$climate[,c(25:36)]) 
      }
      
      # plotting material & plot data frame
      x_ante <- rowMeans(clim_sub)
      pl_df  <- plotting_df(x_ante,post_mean,mod_data)
      
    }
    
    if(mod == 'gaus'){
      
      # mean params
      post_mean <- post_all %>% subset(model == mod )
      
      # gaus
      gaus_w_v  <- dnorm(xx_a, post_mean$sens_mu, 
                               post_mean$sens_sd )
      w_v      <- gaus_w_v / sum(gaus_w_v)
      
      # plotting material & plot data frame
      x_ante <- as.matrix(mod_data$climate) %*% w_v
      pl_df  <- plotting_df(x_ante,post_mean,mod_data)
      
    }
    
    if(mod == 'expp'){
      
      # mean params
      post_mean <- post_all %>% subset(model == mod )
      
      # gaus
      expp_w_v  <- dexppow(xx_a, post_mean$sens_mu, 
                                 post_mean$sens_sd,
                                 expp_beta)
      w_v       <- expp_w_v / sum(expp_w_v)
      
      # plotting material & plot data frame
      x_ante <- as.matrix(mod_data$climate) %*% w_v
      pl_df  <- plotting_df(x_ante,post_mean,mod_data)
      
    }
    
    # expp_n
    if(mod == 'expp_n'){
    
      # post means
      post_mean <- post_all %>% subset(model == mod )
      
      # expp
      expp_w_v <- dexppow(xx, post_mean$sens_mu, 
                              post_mean$sens_sd, 20)
      w_v <- c( ((expp_w_v / sum(expp_w_v))*post_mean$theta_y_1),
                ((expp_w_v / sum(expp_w_v))*post_mean$theta_y_2),
                ((expp_w_v / sum(expp_w_v))*post_mean$theta_y_3) )
     
      # plotting material & plot data frame
      x_ante <- as.matrix(mod_data$climate) %*% w_v
      pl_df  <- plotting_df(x_ante,post_mean,mod_data)
          
    }
    
    if(mod == 'gaus_n'){
    
      # let's start from SAD
      post_mean <- post_all %>% subset(model == mod )
      
      # expp
      expp_w_v <- dnorm(xx, post_mean$sens_mu, 
                            post_mean$sens_sd)
      w_v <- c( ((expp_w_v / sum(expp_w_v))*post_mean$theta_y_1),
                ((expp_w_v / sum(expp_w_v))*post_mean$theta_y_2),
                ((expp_w_v / sum(expp_w_v))*post_mean$theta_y_3) )
     
      # plotting material & plot data frame
      x_ante <- as.matrix(mod_data$climate) %*% w_v
      pl_df  <- plotting_df(x_ante,post_mean,mod_data)
          
    }
    
    if(mod == 'simpl_n'){
      
      # means
      post_mean <- post_all %>% subset(model == mod )
      
      # x antecedent
      w_v    <- c(as.numeric(post_mean[paste0('theta_m_',1:12)]*post_mean$theta_y_1),
                  as.numeric(post_mean[paste0('theta_m_',1:12)]*post_mean$theta_y_2),
                  as.numeric(post_mean[paste0('theta_m_',1:12)]*post_mean$theta_y_3) )
      
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
      ggsave( paste0('results/y_vs_ante/',
                     clim_var_print,'/',
                     spp_name,'_',
                     mod,'.tiff'),
              width = 6.3, height = 5,
              compression="lzw")
  
  }
  
  lapply(c('yr1','yr2','yr3',
           'gaus','expp',
           'gaus_n','expp_n','simpl_n'),plot_spp)
  print(ii)
  
}

  
# 
# # posteriors -----------------------------------------------------------
# par_post <- function(ii){
#   
#   interval    <- NULL
#   pre_chelsa  <- NULL # '_pre_chelsa'
#   clim_var  <- input_df$clim_var[ii]
#   response  <- input_df$response[ii]
#   interval  <- input_df$interval[ii]
#   resp_clim <- paste0("_",response,"_",clim_var)
#   
#   # read lambda/clim data 
#   lam     <- read.csv("C:/cloud/Dropbox/sAPROPOS/all_demog_6tr.csv", 
#                       stringsAsFactors = F) #%>%
#                 #subset( SpeciesAuthor != "Purshia_subintegra" )
# 
#   # summary info 
#   
#   # result folder name
#   res_folder<- paste0("C:/cloud/Dropbox/sAPROPOS/results/moving_windows/",
#                       response,
#                       "/",
#                       clim_var,pre_chelsa,interval) 
#   # summary files names
#   sum_files <- list.files(res_folder)[grep("posterior", list.files(res_folder) )] %>% 
#                   stringr::str_subset(resp_clim)
#   # read files
#   mod_summ  <- lapply(sum_files, function(x) fread(paste0(res_folder,"/",x)) ) %>%
#                   setNames( gsub("posterior", "", sum_files ) ) %>%
#                   setNames( gsub(paste0(resp_clim,".csv"), "", names(.) ) )
#   # all model selection summaries
#   all_sums  <- Map(function(x,y) tibble::add_column(x, species = y, .before = 1), 
#                    mod_summ, names(mod_summ) ) %>% 
#                   # # selec ONLY these model selection variables
#                   # lapply(function(x) x %>% 
#                   #                     dplyr::select(species, 
#                   #                                   model,
#                   #                                   sens_mu, sens_sd, loc, scale, shape,     
#                   #                                   theta_y_1, theta_y_2, theta_y_3, 
#                   #                                   theta_m_1, theta_m_2, theta_m_3, theta_m_4,
#                   #                                   theta_m_5, theta_m_6, theta_m_7, theta_m_8,
#                   #                                   theta_m_9, theta_m_10, theta_m_11, theta_m_12,
#                   #                                   beta, alpha) ) %>% 
#                   bind_rows %>% 
#                   subset( model %in% c("simpl_n", "gev_n", "expp_n") ) %>% 
#                   mutate( clim_var = clim_var )
#         
#   all_sums
# 
# }
# 
#  
# # all models results
# input_df    <- expand.grid( clim_var = c("precip","airt"),
#                             response = "log_lambda", #c("surv","grow","fec",
#                             interval = "",
#                             stringsAsFactors = F)
# 
# # ALL model information
# post_df  <- lapply(1:nrow(input_df), par_post) %>%
#                   bind_rows %>% 
#                   arrange(species) %>% 
#                   mutate( species = gsub('^_','',species) )
# 
# 
# # plot models ----------------------------------------------------------------
# plot_spp <- function(ii){
#   
#   clim_var  <- plot_in$clim_var[ii]
#   mod       <- plot_in$model[ii]
#   spp_name  <- plot_in$species[ii]
#   
#   # store climate info
#   clim      <- clim_l[clim_var][[1]]
#   
#   # parameters
#   response  <- "log_lambda"
#   family    <- "normal"
#   expp_beta <- 20
#   m_back    <- 36    
#   st_dev    <- FALSE
#   xx        <- 1:12
#   
#   # format data 
#   
#   # lambda 
#   spp_resp      <- format_species(spp_name, lam, response)
#   # climate
#   clim_separate <- clim_list(spp_name, clim, spp_resp)
#   clim_detrnded <- lapply(clim_separate, clim_detrend, clim_var)
#   clim_mats     <- Map(clim_long, clim_detrnded, spp_resp, m_back)
#   # model 
#   mod_data          <- lambda_plus_clim(spp_resp, clim_mats, response)
#   mod_data$climate  <- mod_data$climate #/ diff(range(mod_data$climate))
#   
#   
#   # plot model over data --------------------------------------------------------------
#   if(mod == 'simpl_n'){
#     
#     col_n <- c(paste0('theta_y_',1:3),
#                paste0('theta_m_',1:12),
#                'beta', 'alpha')
#     
#     # extract means
#     post_mean <- post_df %>% 
#                     subset( species == spp_name ) %>% 
#                     subset( model == 'simpl_n' ) %>% 
#                     select( col_n ) %>% 
#                     as.matrix %>% 
#                     colMeans
#     
#     # x antecedent
#     w_v    <- c(post_mean[paste0('theta_m_',1:12)]*post_mean['theta_y_1'],
#                 post_mean[paste0('theta_m_',1:12)]*post_mean['theta_y_2'],
#                 post_mean[paste0('theta_m_',1:12)]*post_mean['theta_y_3'])
#     x_ante <- as.matrix(mod_data$climate) %*% w_v
#     b01    <- list( alpha = post_mean['alpha'],
#                     beta  = post_mean['beta'] )
#     pl_df  <- data.frame( y = mod_data$resp$log_lambda,
#                           x = x_ante,
#                           stringsAsFactors = F) 
#   }
#   
#   # gev
#   if(mod == 'gev_n'){
#     
#     # let's start from SAD
#     post_mean <- post_df %>% 
#                     subset( species == spp_name ) %>% 
#                     subset( model == 'gev_n' ) %>% 
#                     select( c('loc','scale','shape',
#                               paste0('theta_y_',1:3),
#                               'alpha', 'beta') ) %>% 
#                     colMeans %>% t %>% 
#                     as.data.frame 
#     
#     # gev
#     gev_w_v  <- dgev(xx, post_mean$loc, post_mean$scale,
#                          post_mean$shape )
#     w_v      <- c( ((gev_w_v / sum(gev_w_v))*post_mean$theta_y_1),
#                    ((gev_w_v / sum(gev_w_v))*post_mean$theta_y_2),
#                    ((gev_w_v / sum(gev_w_v))*post_mean$theta_y_3) )
#     
#     # plotting material  
#     x_ante <- as.matrix(mod_data$climate) %*% w_v
#     b01    <- list( alpha = post_mean['alpha'],
#                     beta  = post_mean['beta'] )
#     pl_df  <- data.frame( y = mod_data$resp$log_lambda,
#                           x = x_ante,
#                           stringsAsFactors = F) 
#   }
#   
#   # expp
#   if(mod == 'expp_n'){
#   
#     # let's start from SAD
#     post_mean <- post_df %>% 
#                     subset( species == spp_name ) %>% 
#                     subset( model == 'expp_n' ) %>% 
#                     select( c('sens_mu','sens_sd',
#                               paste0('theta_y_',1:3),
#                               'alpha', 'beta') ) %>% 
#                     colMeans %>% t %>% 
#                     as.data.frame 
#     
#     # expp
#     expp_w_v <- dexppow(xx, post_mean$sens_mu, 
#                             post_mean$sens_sd, 20)
#     w_v <- c( ((expp_w_v / sum(expp_w_v))*post_mean$theta_y_1),
#               ((expp_w_v / sum(expp_w_v))*post_mean$theta_y_2),
#               ((expp_w_v / sum(expp_w_v))*post_mean$theta_y_3) )
#    
#     # plotting material
#     x_ante <- as.matrix(mod_data$climate) %*% w_v
#     b01    <- list( alpha = post_mean['alpha'],
#                     beta  = post_mean['beta'] )
#     pl_df  <- data.frame( y = mod_data$resp$log_lambda,
#                           x = x_ante,
#                           stringsAsFactors = F) 
#         
#   }
#     
#   # plot it out
#   clim_var_print <- gsub('ip','',clim_var)
#   
#   # plot it out
#   ggplot(pl_df, aes(x, y) ) +
#     geom_point() + 
#     geom_abline( intercept = as.numeric(b01$alpha), 
#                  slope     = as.numeric(b01$beta) ) + 
#     ylab( expression('log('*lambda*')') ) +
#     xlab( 'X antecedent' ) + 
#     ggtitle( gsub('_',' ',spp_name) ) +
#     ggsave( paste0('results/y_vs_ante/',
#                    clim_var_print,'/',
#                    mod,'/',
#                    spp_name,'.tiff'),
#             width = 6.3, height = 5,
#             compression="lzw")
# 
# }
# 
# # 'Eryngium_cuneifolium' 
# 
# # download core data
# r_dir     <- 'C:/cloud/Dropbox/sApropos/'
# lam       <- read.csv(paste0(r_dir,"all_demog_6tr.csv"), 
#                       stringsAsFactors = F)
# clim_l    <- list( airt = data.table::fread(
#                             paste0(r_dir, 'airt_chelsa_hays.csv') ),
#                    precip = data.table::fread(
#                             paste0(r_dir, 
#                                    'precip_chelsa_hays.csv') ) )
# spp       <- lam$SpeciesAuthor %>% 
#                 unique %>% 
#                 Filter(function(x) x !='Daphne_rodriguezii',
#                        .)
#   
# # check all error messages
# plot_in <- expand.grid( model    = c('simpl_n', 'gev_n', 'expp_n'),
#                         clim_var = c('precip', 'airt'),
#                         species  = spp,
#                         stringsAsFactors = F )
# 
# # exponential power distribution
# dexppow <- function(x, mu, sigma, beta) {
#     return((beta / (2 * sigma * gamma (1.0/beta)) ) *
#              exp(-(abs(x - mu)/sigma)^beta));
# }
# 
# # plot it all out
# lapply(1:nrow(plot_in), plot_spp)
