#rm(list=ls())
options(stringsAsFactors=F)
source("C:/CODE/moving_windows/format_data.R")
library(tidyverse)
library(dismo)
library(evd)
library(testthat)
library(data.table)
library(rstan)
library(loo)


# read parameters posterior
par_post <- function(ii){
  
  interval    <- NULL
  pre_chelsa  <- NULL # '_pre_chelsa'
  clim_var  <- input_df$clim_var[ii]
  response  <- input_df$response[ii]
  interval  <- input_df$interval[ii]
  resp_clim <- paste0("_",response,"_",clim_var)
  
  # read lambda/clim data 
  lam     <- read.csv("C:/cloud/Dropbox/sAPROPOS/all_demog_6tr.csv", 
                      stringsAsFactors = F) #%>%
                #subset( SpeciesAuthor != "Purshia_subintegra" )

  # summary info 
  
  # result folder name
  res_folder<- paste0("C:/cloud/Dropbox/sAPROPOS/results/moving_windows/",
                      response,
                      "/",
                      clim_var,pre_chelsa,interval) 
  # summary files names
  sum_files <- list.files(res_folder)[grep("posterior", list.files(res_folder) )] %>% 
                  stringr::str_subset(resp_clim)
  # read files
  mod_summ  <- lapply(sum_files, function(x) fread(paste0(res_folder,"/",x)) ) %>%
                  setNames( gsub("posterior", "", sum_files ) ) %>%
                  setNames( gsub(paste0(resp_clim,".csv"), "", names(.) ) )
  # all model selection summaries
  all_sums  <- Map(function(x,y) tibble::add_column(x, species = y, .before = 1), 
                   mod_summ, names(mod_summ) ) %>% 
                  # # selec ONLY these model selection variables
                  # lapply(function(x) x %>% 
                  #                     dplyr::select(species, 
                  #                                   model,
                  #                                   sens_mu, sens_sd, loc, scale, shape,     
                  #                                   theta_y_1, theta_y_2, theta_y_3, 
                  #                                   theta_m_1, theta_m_2, theta_m_3, theta_m_4,
                  #                                   theta_m_5, theta_m_6, theta_m_7, theta_m_8,
                  #                                   theta_m_9, theta_m_10, theta_m_11, theta_m_12,
                  #                                   beta, alpha) ) %>% 
                  bind_rows %>% 
                  subset( model %in% c("simpl_n", "gev_n", "expp_n") ) %>% 
                  mutate( clim_var = clim_var )
        
  all_sums

}

 
# all models results
input_df    <- expand.grid( clim_var = c("precip","airt"),
                            response = "log_lambda", #c("surv","grow","fec",
                            interval = "",
                            stringsAsFactors = F)

# ALL model information
post_df  <- lapply(1:nrow(input_df), par_post) %>%
                  bind_rows %>% 
                  arrange(species) %>% 
                  mutate( species = gsub('^_','',species) )



# check all error messages
plot_in <- expand.grid( model    = c('simpl_n', 'gev_n', 'expp_n'),
                        clim_var = c('precip', 'airt'),
                        species  = spp,
                        stringsAsFactors = F )

# sampled posterio predictive checks ------------------------------
sppc <- function(ii){
  
  clim_var  <- plot_in$clim_var[ii]
  mod       <- plot_in$model[ii]
  spp_name  <- plot_in$species[ii]
  
  # store climate info
  clim      <- clim_l[clim_var][[1]]
  
  # parameters
  response  <- "log_lambda"
  family    <- "normal"
  expp_beta <- 20
  m_back    <- 36    
  st_dev    <- FALSE
  xx        <- 1:12
  
  # format data 
  
  # lambda 
  spp_resp      <- format_species(spp_name, lam, response)
  # climate
  clim_separate <- clim_list(spp_name, clim, spp_resp)
  clim_detrnded <- lapply(clim_separate, clim_detrend, clim_var)
  clim_mats     <- Map(clim_long, clim_detrnded, spp_resp, m_back)
  # model 
  mod_data          <- lambda_plus_clim(spp_resp, clim_mats, response)
  mod_data$climate  <- mod_data$climate #/ diff(range(mod_data$climate))
  
  
  # plot model over data --------------------------------------------------------------
  if(mod == 'simpl_n'){
    
    col_n <- c(paste0('theta_y_',1:3),
               paste0('theta_m_',1:12),
               'beta', 'alpha', 'y_sd')
    
    # extract means
    post_smpl <- post_df %>% 
                    subset( species == spp_name ) %>% 
                    subset( model == 'simpl_n' ) %>% 
                    select( col_n ) %>% 
                    as.matrix %>% 
                    # random sample from posterior
                    .[100,]
    
    # x antecedent
    w_v    <- c(post_smpl[paste0('theta_m_',1:12)]*post_smpl['theta_y_1'],
                post_smpl[paste0('theta_m_',1:12)]*post_smpl['theta_y_2'],
                post_smpl[paste0('theta_m_',1:12)]*post_smpl['theta_y_3'])
    x_ante <- as.matrix(mod_data$climate) %*% w_v
    b01    <- list( alpha = post_smpl['alpha'],
                    beta  = post_smpl['beta'],
                    y_sd  = post_smpl['y_sd'] )
    pl_df  <- data.frame( y     = mod_data$resp$log_lambda,
                          x     = x_ante,
                          stringsAsFactors = F ) %>% 
                  mutate( alpha = post_smpl['alpha'],
                          beta  = post_smpl['beta'],
                          y_sd  = post_smpl['y_sd'] ) 
  }
  
  # gev
  if(mod == 'gev_n'){
    
    # let's start from SAD
    post_smpl <- post_df %>% 
                    subset( species == spp_name ) %>% 
                    subset( model == 'gev_n' ) %>% 
                    select( c('loc','scale','shape',
                              paste0('theta_y_',1:3),
                              'alpha', 'beta', 'y_sd') ) %>% 
                    as.data.frame %>% 
                    # random sample from posterior
                    .[100,]
    
    # gev
    gev_w_v  <- dgev(xx, post_smpl$loc, post_smpl$scale,
                         post_smpl$shape )
    w_v      <- c( ((gev_w_v / sum(gev_w_v))*post_smpl$theta_y_1),
                   ((gev_w_v / sum(gev_w_v))*post_smpl$theta_y_2),
                   ((gev_w_v / sum(gev_w_v))*post_smpl$theta_y_3) )
    
    # plotting material  
    x_ante <- as.matrix(mod_data$climate) %*% w_v
    pl_df  <- data.frame( y     = mod_data$resp$log_lambda,
                          x     = x_ante,
                          stringsAsFactors = F ) %>% 
                  mutate( alpha = post_smpl$alpha,
                          beta  = post_smpl$beta,
                          y_sd  = post_smpl$y_sd ) 
  }
  
  # expp
  if(mod == 'expp_n'){
  
    # let's start from SAD
    post_smpl <- post_df %>% 
                    subset( species == spp_name ) %>% 
                    subset( model == 'expp_n' ) %>% 
                    select( c('sens_mu','sens_sd',
                              paste0('theta_y_',1:3),
                              'alpha', 'beta', 'y_sd') ) %>% 
                    # random sample from posterior
                    .[100,]
    
    # expp
    expp_w_v <- dexppow(xx, post_smpl$sens_mu, 
                            post_smpl$sens_sd, 20)
    w_v <- c( ((expp_w_v / sum(expp_w_v))*post_smpl$theta_y_1),
              ((expp_w_v / sum(expp_w_v))*post_smpl$theta_y_2),
              ((expp_w_v / sum(expp_w_v))*post_smpl$theta_y_3) )
   
    # plotting material
    x_ante <- as.matrix(mod_data$climate) %*% w_v
    pl_df  <- data.frame( y     = mod_data$resp$log_lambda,
                          x     = x_ante,
                          stringsAsFactors = F ) %>% 
                  mutate( alpha = post_smpl$alpha,
                          beta  = post_smpl$beta,
                          y_sd  = post_smpl$y_sd ) 
  }
    
  # Bayesian p-value calculation: sample posterior
  sim_post <- function(ss){
    pl_df <- pl_df %>% 
              mutate( y_pred  = alpha + beta*x ) %>%
              mutate( y_rep   = rnorm(y_pred,y_pred,y_sd) ) %>% 
              mutate( dis_dat = (y     - y_pred)^2,
                      dis_sim = (y_rep - y_pred)^2 )
    
    # Discrepancy (sum of squared residuals)
    discr <- pl_df %>% 
                select(dis_dat, dis_sim) %>% 
                as.matrix %>%
                colSums
     
    if( discr['dis_sim'] > discr['dis_dat'] ) 1 else 0
            
  }
  
  res   <- sapply(1:1000, sim_post)
  p_val <- sum(res) / 1000

  data.frame( clim_var  = plot_in$clim_var[ii],
              mod       = plot_in$model[ii],
              spp_name  = plot_in$species[ii],
              p_val     = p_val,
              stringsAsFactors = F)

}

cia_l1 <- lapply(101:200, sppc)
cia_l <- lapply(101:100, sppc)
cia_l1 %>% bind_rows
cia_l %>% bind_rows

lapply(201:204, sppc) %>% bind_rows
