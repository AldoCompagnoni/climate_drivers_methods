# Read data
# Read model results
# plot spp. specific results
rm(list=ls())
options(stringsAsFactors=F)
source("code/format_data.R")
library(tidyverse)
library(mgcv)
library(testthat)
library(rstan)
library(evd)
library(gridExtra)
library(grid)

# climate predictor, response, months back, max. number of knots
resp_l    <- c('log_lambda', 'surv', 'grow', 'fec')
clim_var  <- "precip"
m_back    <- 36    
st_dev    <- FALSE

# read data -----------------------------------------------------------------------------------------

# model data
lam       <- read.csv("data/all_demog_updt.csv", stringsAsFactors = F)
clim      <- data.table::fread(paste0('data/',clim_var,"_chelsa_prism_hays_2014.csv"),  stringsAsFactors = F)
spp       <- lam$SpeciesAuthor %>% unique

for(rr in 1:1){
  
  response <- resp_l[rr]
    
  # model data to update species list
  sum_f_l   <- grep('posterior',
                    list.files(paste0('results/mod_sel/',clim_var)),
                    value=T) %>% 
                 grep(paste0('_',response,'_'), ., value=T)
  spp_names <- gsub(paste0('posterior_|.csv|',
                      paste0(response,'_'),sep='|'),'',sum_f_l)
  spp       <- intersect(spp,spp_names)
  
  # format data --------------------------------------------------------------------------------------
  
  # set up model "family" based on response
  if( response == "surv" | response == "grow" )             family = "beta" 
  if( response == "fec" )                                   family = "gamma"
  if( grepl("PreRep", response) | grepl("Rep", response) )  family = "beta"
  if( response == "rho" | response == "react_fsa" )         family = "gamma"
  if( response == "log_lambda")                             family = "normal"
  
  expp_beta     <- 20
  
  for(ii in 1:33 ){ #2length(spp)
    
    # set species
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
    post_all <- data.table::fread( paste0('results/mod_sel/',clim_var,'/',
                                   grep(spp_name,sum_f_l,value=T)) ) 
    
    
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
    plot_spp <- function(mod, response){
      
      # parameters
      expp_beta <- 20
      m_back    <- 36    
      xx        <- 1:12
      xx_a      <- 1:36
      slices    <- seq(1,6000,length.out=200) %>% round()
      
      # produce data frame to plot data
      plotting_df <- function(slices,x_ante,post_mean,mod_data,response){
        out_df <- data.frame( y = mod_data$resp[,response],
                    x = x_ante,
                    stringsAsFactors = F) %>% 
                    mutate( alpha = post_mean$alpha[slices],
                            beta  = post_mean$beta[slices] ) 
  
        return(out_df)
           
      }
      
      # produce means and posterior based on model
      mean_post <- function(mod){
        # assign objects to the parent environment
        post_m    <<- post_all %>% subset(model == mod )
        post_mean <<- post_m %>% 
                      group_by( model ) %>% 
                      summarise_all( mean ) %>% 
                      ungroup
      }
      
      # quantiles of year weights
      quant_ind_w <- function(post_in,vec){
          
          lapply(vec,
                 function(x) quantile(as.data.frame(post_in)[,x],
                                      prob=c(0.025,0.5,0.975))
                ) %>% 
            do.call(rbind, .) %>% 
            as.data.frame %>% 
            tibble::add_column(.before=1, index=1:length(vec) ) %>% 
            setNames( c('index','lwr','mid','upr') ) %>% 
            mutate( index = factor( as.character(index), 
                                    levels = as.character(index) ) )
            
      }
    
      # plot model over data --------------------------------------------------------------
      if( grepl('yr', mod) ){
        
        # load relevant params
        mean_post( mod )
        
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
        pl_l   <- lapply(slices,plotting_df,
                         x_ante,post_m,mod_data,response)
        pl_df  <- plotting_df(1,x_ante,post_mean,mod_data,response)
       
        w_v_l  <- NULL
      
      }
      
      if(mod == 'gaus'){
        
        # load relevant params
        mean_post( mod )
        
        # x_antecedent means
        gaus_w_v  <- dnorm(xx_a, post_mean$sens_mu, 
                                 post_mean$sens_sd )
        w_v       <- gaus_w_v / sum(gaus_w_v)
        
        # x_antecedent posterior
        ante_post <- function(ii,post_in){
          gaus_w_v  <- dnorm(xx_a, post_in$sens_mu[ii], 
                                   post_in$sens_sd[ii] )      
          (gaus_w_v / sum(gaus_w_v)) %>% 
            as.data.frame %>% 
            mutate( x = 1:length(gaus_w_v) ) %>% 
            setNames( c('w','x') ) %>% select(x,w)
        }
        
        # MEANS: plotting material
        x_ante <- as.matrix(mod_data$climate) %*% w_v
        pl_df  <- plotting_df(1,x_ante,post_mean,mod_data, response)
        
        # posterior
        pl_l   <- lapply(slices, plotting_df, x_ante, post_m, mod_data, response)
        w_v_l  <- lapply(slices, ante_post, post_m)
        
      }
      
      # if(mod == 'expp'){
      #   
      #   # load relevant params
      #   mean_post( mod )
      #   
      #   # gaus
      #   expp_w_v  <- dexppow(xx_a, post_mean$sens_mu, 
      #                              post_mean$sens_sd,
      #                              expp_beta)
      #   w_v       <- expp_w_v / sum(expp_w_v)
      #   
      #   # x_antecedent posterior
      #   ante_post <- function(ii, post_in){
      #     gaus_w_v  <- dexppow(xx_a, post_in$sens_mu[ii], 
      #                                post_in$sens_sd[ii],
      #                                expp_beta)      
      #     (gaus_w_v / sum(gaus_w_v)) %>% 
      #       as.data.frame %>% 
      #       mutate( x = 1:length(gaus_w_v) ) %>% 
      #       setNames( c('w','x') ) %>% select(x,w)
      #   }
      #   
      #   # MEANS: plotting material
      #   x_ante <- as.matrix(mod_data$climate) %*% w_v
      #   pl_df  <- plotting_df(1,x_ante,post_mean,mod_data,response)
      #   
      #   # posterior
      #   pl_l   <- lapply(slices, plotting_df, x_ante, post_m, mod_data, response)
      #   w_v_l  <- lapply(slices, ante_post, post_m)
      #   
      # }
      
      # # moving betas
      # if(mod %in% c('movb', 'movb_h') ){
      #   
      #   # load relevant params
      #   mean_post( mod )
      #   
      #   # consider "mu beta" as "beta"
      #   post_m    <- select(post_m,-beta) %>% rename(beta=mu_beta)
      #   post_mean <- select(post_mean,-beta) %>% rename(beta=mu_beta)
      #   
      #   # # MEANS: plotting material
      #   b_names <- paste0('beta_',1:36)
      #   # x_ante  <- as.matrix(mod_data$climate) %*% as.numeric(post_mean[,b_names])
      #   # pl_df  <- plotting_df(1,x_ante,post_mean,mod_data)
      #   # 
      #   # # posterior
      #   # pl_l   <- lapply(slices, plotting_df, x_ante, post_m, mod_data)
      #   qnt_m  <- quant_ind_w(post_m,b_names)
      #   
      # }
        
      # # moving betas nested
      # if(mod %in% c('movb_n','movb_h_n') ){
      #   
      #   # load relevant params
      #   mean_post( mod )
      #   
      #   # consider "mu beta" as "beta"
      #   post_m    <- select(post_m,-beta) %>% rename(beta=mu_beta)
      #   post_mean <- select(post_mean,-beta) %>% rename(beta=mu_beta)
      #   
      #   # # MEANS: plotting material
      #   b_names <- paste0('beta_',1:12)
      #   # x_ante  <- as.matrix(mod_data$climate) %*% as.numeric(post_mean[,b_names])
      #   # pl_df  <- plotting_df(1,x_ante,post_mean,mod_data)
      #   # 
      #   # # posterior
      #   # pl_l   <- lapply(slices, plotting_df, x_ante, post_m, mod_data)
      #   qnt_m  <- quant_ind_w(post_m,b_names)
      #   qnt_yr <- quant_ind_w(post_m,paste0('theta_y_',1:3))
      #     
      # }
      
      # # expp_n
      # if(mod == 'expp_n'){
      # 
      #   # load relevant params
      #   mean_post( mod )
      #   
      #   # expp
      #   expp_w_v <- dexppow(xx, post_mean$sens_mu, 
      #                           post_mean$sens_sd, 20)
      #   w_v <- c( ((expp_w_v / sum(expp_w_v))*post_mean$theta_y_1),
      #             ((expp_w_v / sum(expp_w_v))*post_mean$theta_y_2),
      #             ((expp_w_v / sum(expp_w_v))*post_mean$theta_y_3) )
      #  
      #   # x_antecedent posterior
      #   ante_post <- function(ii, post_in){
      #     expp_w_v  <- dexppow(xx, post_in$sens_mu[ii], 
      #                              post_in$sens_sd[ii],
      #                              expp_beta)
      #     w_v <- (expp_w_v / sum(expp_w_v))
      #     w_v %>% 
      #       as.data.frame %>% 
      #       mutate( x = 1:length(expp_w_v) ) %>% 
      #       setNames( c('w','x') ) %>% select(x,w)
      #   }
      #   
      #   
      #   # MEANS: plotting material
      #   x_ante <- as.matrix(mod_data$climate) %*% w_v
      #   pl_df  <- plotting_df(1,x_ante,post_mean,mod_data,response)
      #   
      #   # posteriors 
      #   pl_l   <- lapply(slices,plotting_df,x_ante,post_m,mod_data,response)
      #   w_v_l  <- lapply(slices, ante_post, post_m)
      #   qnt_yr <- quant_ind_w(post_m,paste0('theta_y_',1:3))
      #   
      # }
      
      # if(mod == 'gaus_n'){
      # 
      #   # load parameters and posterior
      #   mean_post( mod )
      #   
      #   # expp
      #   expp_w_v <- dnorm(xx, post_mean$sens_mu, 
      #                         post_mean$sens_sd)
      #   w_v <- c( ((expp_w_v / sum(expp_w_v))*post_mean$theta_y_1),
      #             ((expp_w_v / sum(expp_w_v))*post_mean$theta_y_2),
      #             ((expp_w_v / sum(expp_w_v))*post_mean$theta_y_3) )
      #  
      #   # posterior of antecedent
      #   ante_post <- function(ii, post_in){
      #     gaus_w_v  <- dnorm(xx, post_in$sens_mu[ii], 
      #                            post_in$sens_sd[ii] )
      #     w_v <- (gaus_w_v / sum(gaus_w_v))
      #     w_v %>% 
      #       as.data.frame %>% 
      #       mutate( x = 1:length(xx) ) %>% 
      #       setNames( c('w','x') ) %>% select(x,w)
      #   }
      #   
      #   # plotting material & plot data frame
      #   x_ante <- as.matrix(mod_data$climate) %*% w_v
      #   pl_df  <- plotting_df(1,x_ante,post_mean,mod_data,response)
      #     
      #   # posteriors 
      #   pl_l   <- lapply(slices, plotting_df,x_ante,post_m,mod_data,response)
      #   w_v_l  <- lapply(slices, ante_post, post_m)
      #   qnt_yr <- quant_ind_w(post_m,paste0('theta_y_',1:3))
      #     
      # }
      
      if(mod == 'simpl_n'){
        
        # means
        mean_post( mod )
        
        # x antecedent
        w_v    <- c(as.numeric(post_mean[paste0('theta_m_',1:12)]*post_mean$theta_y_1),
                    as.numeric(post_mean[paste0('theta_m_',1:12)]*post_mean$theta_y_2),
                    as.numeric(post_mean[paste0('theta_m_',1:12)]*post_mean$theta_y_3) )
        
        # plotting material & plot data frame
        x_ante <- as.matrix(mod_data$climate) %*% w_v
        pl_df  <- plotting_df(1,x_ante,post_mean,mod_data,response)
        
        # posteriors 
        pl_l   <- lapply(slices, plotting_df,x_ante,post_m,mod_data,response)
        qnt_m  <- quant_ind_w(post_m, paste0('theta_m_',1:12))
        qnt_yr <- quant_ind_w(post_m,paste0('theta_y_',1:3))
        
      }
      
      # plot it out
      clim_var_print <- gsub('ip','',clim_var)
      
      
      # plot one: y_vs_x ----------------------------------
      if( mod %in% c('movb',  'movb_h',
                     'movb_n','movb_h_n') ){
        
        y_x <-  grid::textGrob('NA')
        
      }else{
        
        # sequence of xs
        xs <- seq(min(pl_df$x),max(pl_df$x),length.out=100)
        
        p1 <- ggplot(pl_df, aes(x, y) ) +
                geom_point() +
                ylab( response ) +
                xlab( 'X antecedent' ) + 
                ggtitle( spp_name )
        
        for(ii in 1:200){
        
          y_raw <- unique(pl_l[[ii]]$alpha) +
                   (unique(pl_l[[ii]]$beta) * xs)
          
          if(response == 'log_lambda')        y_hat <- y_raw 
          if(response %in% c('surv','grow') ) y_hat <- boot::inv.logit(y_raw)
          if(response %in% c('fec') )         y_hat <- exp(y_raw)
          
          plot_ii <- data.frame( x = xs,
                                 y = y_hat,
                                 stringsAsFactors = F)
          
          p1 <- p1 +
                geom_line(data=plot_ii,
                          aes(x=x,y=y),
                          color = 'grey' )
                
        }
          
        # y vs. antecedent
        y_raw <- unique(pl_df$alpha) +
                 (unique(pl_df$beta) * xs)
        
        if(response == 'log_lambda')        y_hat <- y_raw 
        if(response %in% c('surv','grow') ) y_hat <- boot::inv.logit(y_raw)
        if(response %in% c('fec') )         y_hat <- exp(y_raw)
        
        plot_ii <- data.frame( x = xs,
                               y = y_hat,
                               stringsAsFactors = F )
        
        y_x   <- p1 + 
                  geom_point( data = pl_df,
                              aes(x=x, y=y) ) +
                  geom_line( data = plot_ii, aes(x=x, y=y),size=1 )
        
      }
      
      # plot month weights/betas -----------------------------
      
      if( mod %in% c('gaus','expp','gaus_n','expp_n') ){
        
        # initiate weights plot
        p2 <- ggplot(data = w_v_l[[1]],
                     aes(x=x, y=w) ) +
                geom_line( color = 'grey',
                           alpha = 0.7 )
          
        for(ii in 2:200){
          p2 <- p2 +
                geom_line(data=w_v_l[[ii]],
                          aes( x=x, y=w),
                          color = 'grey',
                          alpha = 0.7 )
        }
        
        mw_p <- p2
        
      }
      
      if( mod %in% c('simpl_n') ){
        
        mw_p <- ggplot(qnt_m) +
                geom_pointrange(aes(x=index,
                                    y=mid,
                                    ymin=lwr,
                                    ymax=upr)) +
                ylab('weight') +
                xlab('month')  
        
      }
      
      if( mod %in% c('movb',  'movb_h',
                     'movb_n','movb_h_n') ){
        
        mw_p <- ggplot(qnt_m) +
                geom_pointrange(aes(x=index,
                                    y=mid,
                                    ymin=lwr,
                                    ymax=upr)) +
                geom_hline( yintercept=0,lty=2) +
                ylab('beta') +
                xlab('month') +
                ggtitle( spp_name )
        
      }
      
      if( mod %in% paste0('yr',1:3) ){
        mw_p <- grid::textGrob('NA')
      }
      
      # plot year weights -----------------------------
      if( mod %in% c('gaus_n','expp_n','simpl_n',
                     'movb_n','movb_h_n') ){
      
        yw_p <- ggplot(qnt_yr) +
          coord_cartesian(ylim = c(0,1) ) + 
          geom_pointrange(aes(x=index,
                              y=mid,
                              ymin=lwr,
                              ymax=upr)) +
          ylab('weight') +
          xlab('year')
          
      }else{
        yw_p <- textGrob('NA')
      }
      
      # kernel density of beta ------------------------
      if( mod %in% c( 'yr1','yr2','yr3',
                      'gaus',  'expp',
                      'gaus_n','expp_n','simpl_n') ){
        
        dens_obj <- post_all %>% 
                      subset( model == mod) %>% 
                      .$beta %>% density
        
        dens_df  <- data.frame(x = dens_obj$x,
                               y = dens_obj$y)
          
        
        kd_p <- subset(post_all, model == mod) %>% 
          ggplot() +
          geom_histogram(aes(x=beta, y=..density..)) +
          geom_line(data=dens_df,aes(x,y),lwd=1.5,
                    color='brown') +
          geom_vline( xintercept = 0, lty=2, lwd=1) +
          ylab( expression('Kernel density of '*beta) ) +
          xlab( expression(beta*' value') )
      }else{
        kd_p <- textGrob('NA')
      }
      
      # final plot -------------------------------------
      out_p <- grid.arrange(y_x,mw_p,yw_p,kd_p,ncol=2)
      ggsave(file = 
             paste0('results/y_vs_ante/',
                     clim_var_print,'/',
                     response,'/',
                     spp_name,'_',
                     mod,'.tiff'),
            plot = out_p,
            width = 6.3, height = 6.3,
            compression="lzw")
      
    }
    
    lapply(c('yr1','yr2','yr3',
             'gaus',  'expp',
             'gaus_n','expp_n','simpl_n',
             'movb',  'movb_h',
             'movb_n','movb_h_n'), plot_spp, response)
    print(ii)
    
  }

}
