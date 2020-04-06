# Perform Bayesian checks ON LOG_LAMBDA only (for now!)
# IMPORTANT: script works one response var. (airt/precip) at a time.
# therefore, need be at least work for both airt/precip
rm(list=ls())
options(stringsAsFactors=F)
source("code/format_data.R")
library(dplyr)
library(tidyr)
library(mgcv)
library(testthat)
library(rstan)
library(evd)
library(gridExtra)
library(grid)

# read data
lam       <- read.csv("data/all_demog_updt.csv", stringsAsFactors = F)
spp       <- lam$SpeciesAuthor %>% unique

# climate predictor, response, months back, max. number of knots
resp_l    <- c('log_lambda')
clim_var  <- "airt"
m_back    <- 36    
st_dev    <- FALSE

# read in the climate data BEFORE doing the analysis
clim      <- data.table::fread(paste0('data/',clim_var,"_chelsa_prism_hays_2014.csv"),  
                               stringsAsFactors = F)

# where are the files?
mod_dir   <- 'E:/work/sApropos/2019.6.11/'

# file list
sum_f_l   <- grep('posterior_',
                  list.files(mod_dir),
                  value=T) %>%
               grep(clim_var, ., value=T) %>% 
               grep(resp_l, ., value=T)

# path list
path_l    <- paste0(mod_dir, sum_f_l)

# extract species names
spp       <- gsub("posterior_", "", sum_f_l ) %>%
                gsub('array_vr-[0-9]{7}-[0-9]{1,2}_', "", . ) %>%
                gsub('.csv', "", . ) %>%
                gsub( paste0('_',resp_l,'_',clim_var),'', .)
 
# analysis cases  
cases_df  <- expand.grid( resp = resp_l[1],
                          spec = spp,
                          stringsAsFactors = F )


# calculate the Bayesia p-values
bayes_p_val <- function( ii ){
  
  # set up characteristics for the computations
  response <- cases_df$resp[ii]
  spp_name <- cases_df$spec[ii]
  
  # format data --------------------------------------------------------------------------------------
  
  # set up model "family" based on response
  if( response == "surv" | response == "grow" )             family = "beta" 
  if( response == "fec" )                                   family = "gamma"
  if( grepl("PreRep", response) | grepl("Rep", response) )  family = "beta"
  if( response == "rho" | response == "react_fsa" )         family = "gamma"
  if( response == "log_lambda")                             family = "normal"
  
  expp_beta     <- 20
  
   
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
  post_all <- data.table::fread( path_l[ii] )
  
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
  
  # across models, set up data frames for simulations -----------------------
  
  # models 
  mod_l    <- post_all %>% .$model %>% unique
  post_l   <- lapply(mod_l, function(x) 
                              post_all %>% 
                              subset( model == x ) %>% 
                              .[100,]
                     ) %>% 
                setNames( mod_l ) 
  
  # produce data frame to plot data
  sppc_df <- function(x_ante,post_smpl,mod_data){
    data.frame( y = mod_data$resp$log_lambda,
                x = x_ante,
                stringsAsFactors = F ) %>% 
      mutate( alpha = post_smpl$alpha,
              beta  = post_smpl$beta,
              y_sd  = post_smpl$y_sd )
  }   
  
  
  # sampled posterior predictive checks across datasets
  sppc <- function(mod, post_l, mod_data){
    
    # needed variables in enclosing environment
    expp_beta <- 20
    m_back    <- 36    
    xx        <- 1:12
    xx_a      <- 1:36
  
    # null model
    if( mod == 'ctrl1'){
      
      x_ante <- 0 # no need for this in null model
      pl_df  <- sppc_df(x_ante, post_l[mod][[1]], mod_data) %>% 
                  mutate( beta = 0 )
    
    }
    
    # function to set up data frames
    if( grepl('yr', mod) ){
        
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
      pl_df  <- sppc_df(x_ante, post_l[mod][[1]], mod_data)
      
    }
    
    # gaussian model
    if(mod == 'gaus'){
      
      # mean params
      post_smpl <- post_l[mod][[1]]
      
      # gaus
      gaus_w_v  <- dnorm(xx_a, post_smpl$sens_mu, post_smpl$sens_sd )
      w_v       <- gaus_w_v / sum(gaus_w_v)
      
      # plotting material & plot data frame
      x_ante <- as.matrix(mod_data$climate) %*% w_v
      pl_df  <- sppc_df(x_ante, post_l[mod][[1]], mod_data)
      
    }
    
    # Generalized Normal model
    if(mod == 'expp'){
      
      # mean params
      post_smpl <- post_l[mod][[1]]
      
      # expp
      expp_w_v  <- dexppow(xx_a, post_smpl$sens_mu, 
                                 post_smpl$sens_sd,
                                 expp_beta)
      w_v       <- expp_w_v / sum(expp_w_v)
      
      # plotting material & plot data frame
      x_ante <- as.matrix(mod_data$climate) %*% w_v
      pl_df  <- sppc_df(x_ante,post_l[mod][[1]],mod_data)
    
    }
    
    # Nested Gaussian
    if(mod == 'gaus_n'){
    
      # let's start from SAD
      post_smpl <- post_l[mod][[1]]
      
      # gaus
      gaus_w_v <- dnorm(xx, post_smpl$sens_mu, post_smpl$sens_sd)
      w_v      <- c( ((gaus_w_v / sum(gaus_w_v))*post_smpl$theta_y_1),
                     ((gaus_w_v / sum(gaus_w_v))*post_smpl$theta_y_2),
                     ((gaus_w_v / sum(gaus_w_v))*post_smpl$theta_y_3) )
     
      # plotting material & plot data frame
      x_ante <- as.matrix(mod_data$climate) %*% w_v
      pl_df  <- sppc_df(x_ante, post_l[mod][[1]], mod_data)
        
    }  
    
    # Nested Generalized Normal
    if(mod == 'expp_n'){
    
      # let's start from SAD
      post_smpl <- post_l[mod][[1]]
      
      # expp
      expp_w_v <- dexppow(xx, post_smpl$sens_mu, 
                              post_smpl$sens_sd, 
                              expp_beta)
      w_v <- c( ((expp_w_v / sum(expp_w_v))*post_smpl$theta_y_1),
                ((expp_w_v / sum(expp_w_v))*post_smpl$theta_y_2),
                ((expp_w_v / sum(expp_w_v))*post_smpl$theta_y_3) )
     
      # plotting material & plot data frame
      x_ante <- as.matrix(mod_data$climate) %*% w_v
      pl_df  <- sppc_df(x_ante, post_smpl, mod_data)
        
    }  
    
    # SAM
    if(mod == 'simpl_n'){
    
      # let's start from SAD
      post_smpl <- post_l[mod][[1]]
      
      # x antecedent
      w_v    <- c(as.numeric(post_smpl[,paste0('theta_m_',1:12)]*post_smpl$theta_y_1),
                  as.numeric(post_smpl[,paste0('theta_m_',1:12)]*post_smpl$theta_y_2),
                  as.numeric(post_smpl[,paste0('theta_m_',1:12)]*post_smpl$theta_y_3) )
      
      # plotting material & plot data frame
      x_ante <- as.matrix(mod_data$climate) %*% w_v
      pl_df  <- sppc_df(x_ante, post_l[mod][[1]], mod_data)
          
    }  
    
    # moving beta hierarchical
    if( mod == 'mb_h' | mod == 'mb' ){
    
      # let's start from SAD
      post_smpl <- post_l[mod][[1]]
      
      # plotting material & plot data frame
      x_ante    <- rep(1, nrow(mod_data$climate) )
      b_mat     <- matrix(select( post_l[mod][[1]], 
                              paste0('beta_',1:36) ) %>% unlist,
                          36, 1)
      c_mat     <- as.matrix(mod_data$climate)
      mock_post <- data.frame( alpha = post_l[mod][[1]]$alpha,
                               beta  = c_mat %*% b_mat,
                               y_sd  = post_l[mod][[1]]$y_sd ) #%>% 
                      # test computations are correct
                      # mutate( yhat   = select(post_l[mod][[1]], 
                      #                         paste0('yhat_',1:31) ) %>% unlist,
                      #         pred   = alpha + beta )
                    
      # put out plot df 
      pl_df  <- sppc_df(x_ante, mock_post, mod_data)
          
    }  
    
    # nested moving beta hierarchical
    if( mod == 'mb_h_n' | mod == 'mb_n' ){
    
      # let's start from SAD
      post_smpl <- post_l[mod][[1]]
      
      # plotting material & plot data frame
      x_ante    <- rep(1, nrow(mod_data$climate) )
      b_m_mat   <- matrix(select( post_l[mod][[1]], 
                                  paste0('beta_',1:12) ) %>% unlist,
                          12, 1)
      
      b_mat     <- do.call( rbind, 
                            list( b_m_mat * post_l[mod][[1]]$theta_y_1,
                                  b_m_mat * post_l[mod][[1]]$theta_y_2,
                                  b_m_mat * post_l[mod][[1]]$theta_y_3 ) )
      
      c_mat     <- as.matrix(mod_data$climate)
      mock_post <- data.frame( alpha = post_l[mod][[1]]$alpha,
                               beta  = c_mat %*% b_mat,
                               y_sd  = post_l[mod][[1]]$y_sd ) #%>% 
                      # test computations are correct
                      # mutate( yhat   = select(post_l[mod][[1]], 
                      #                         paste0('yhat_',1:31) ) %>% unlist,
                      #         pred   = alpha + beta )
                    
      # put out plot df 
      pl_df  <- sppc_df(x_ante, mock_post, mod_data)
          
    }
    
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
  
    data.frame( mod       = mod,
                spp_name  = spp_name,
                p_val     = p_val,
                stringsAsFactors = F)
  
  }

  # calculate Bayesian p-value
  lapply(c('ctrl1', 
           'yr1','yr2','yr3',
           'gaus',     'expp',
           'gaus_n',   'expp_n','simpl_n',
           'mb_h',   'mb',
           'mb_h_n', 'mb_n'), 
          sppc, post_l, mod_data) %>% 
    bind_rows %>% 
    mutate( response = cases_df$resp[ii],
            clim_var = clim_var )

}

# re-read in the climate data BEFORE doing the analysis
clim  <- data.table::fread(paste0('data/',clim_var,"_chelsa_prism_hays_2014.csv"),  
                           stringsAsFactors = F)

# calculate p-values
bayes_p_l <- lapply(1:nrow(cases_df), bayes_p_val)
bayes_p   <- bind_rows(bayes_p_l) 

# save run (takes too long to let it go!)
saveRDS(bayes_p, paste0('results/checks/',
                        clim_var,'_p_val.rds'),
        )


# make graphs -----------------------------------------------------

# p_value by model 
ggplot(bayes_p, 
       aes( x = mod,
            y = p_val ) ) +
  geom_point( position = position_jitter(w=0.06,h=0),
              alpha = 0.5) +
  geom_hline( yintercept = 0.025, lty=2  ) +
  geom_hline( yintercept = 0.975, lty=2  ) +
  theme( axis.text.x = element_text(angle=90) ) +
  ggsave( paste0('results/checks/p_val_',
                 substr(clim_var,1,4),
                 '_llam_mod.tiff'),
          width=6.3, height=6.3, compression='lzw')
  
# p-value by species
bayes_p %>% 
  mutate( spp_name = substr(spp_name, 1, 15) ) %>% 
  ggplot(aes( x = spp_name,
              y = p_val ) ) +
  geom_point( position = position_jitter(w=0.06,h=0),
              alpha = 0.5) +
  geom_hline( yintercept = 0.025, lty=2  ) +
  geom_hline( yintercept = 0.975, lty=2  ) +
  theme( axis.text.x = element_text(angle=90) ) +
  ggsave( paste0('results/checks/p_val_',
                 substr(clim_var,1,4),
                 'llam_spp.tiff'),
          width=6.3, height=6.3, compression='lzw')
  

# make one big graph! -----------------------------------------

# order of models
mod_ordr <- c('ctrl1',  'yr1', 'yr2', 'yr3',
              'gaus',   'expp',
              'gaus_n', 'expp_n', 'simpl_n',
              'mb_h',    'mb',
              'mb_h_n',  'mb_n')


bayes_at  <- readRDS( paste0('results/checks/airt_p_val.rds') ) %>%
                mutate( mod = as.character(mod) ) %>% 
                mutate( mod = factor(mod, levels = mod_ordr) )
bayes_pr  <- readRDS( paste0('results/checks/precip_p_val.rds') ) %>%
                mutate( mod = as.character(mod) ) %>%  
                mutate( mod = factor(mod, levels = mod_ordr) )
bayes_all <- bind_rows( bayes_at, bayes_pr ) %>% 
                mutate( clim_var = replace(clim_var, 
                                           clim_var == 'airt',
                                           'Air temperature') ) %>% 
                mutate( clim_var = replace(clim_var, 
                                           clim_var == 'precip',
                                           'Precipitation') ) 


# p_value by model 
ggplot(bayes_all, 
       aes( x = mod,
            y = p_val,
            color = clim_var) ) +
  geom_point( position = position_jitter(w=0.06,h=0),
              alpha = 0.5) +
  geom_hline( yintercept = 0.025, lty=2  ) +
  geom_hline( yintercept = 0.975, lty=2  ) +
  scale_color_viridis_d() +
  theme( axis.text.x = element_text( angle=90,
                                     hjust=1,
                                     vjust=0.5),
         legend.title = element_text("Clim. var.") ) +
  ylab( 'Bayesian P-value' ) +
  xlab( 'Model' ) +
  labs( color = 'Predictor') +
  ggsave( 'results/checks/p_val_llam_mod.tiff',
          width=6.3, height=4, compression='lzw')

#   
(sum(bayes_all$p_val > 0.975) + sum(bayes_all$p_val < 0.025)) / nrow(bayes_all)
