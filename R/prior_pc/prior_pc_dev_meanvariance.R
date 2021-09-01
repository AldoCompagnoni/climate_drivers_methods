library(rstan)
library(MASS)
library(tidyverse)
library(gtools)
library(ggjoy)
library(ggridges)
library(bayesplot)
library(extraDistr)
library(truncdist)
library(mvtnorm)
library(betareg)
library(modEvA)


# General parameters
it        <- 60000
gam_p     <- 0.01

# n_time = 9 or 30 (bimodal distribution)
set.seed( 1776 + 4)
sig       <- matrix(0,12,12)
diag(sig) <- 1
x_9       <- MASS::mvrnorm(9,  rep(0,12), Sigma=sig)
x_30      <- MASS::mvrnorm(30, rep(0,12), Sigma=sig) 


# NORMAL PROCESS ---------------------------------------------------------------

# YEAR (CLIMATE SUMMARY) MODEL simulation

# replicate yearly anomalies by posterior sample size
clim_m   <- rowMeans( x_9 ) %>% 
              t %>% 
              as.data.frame %>% 
              slice( rep(1:n(), each = it) )

# for merging down below
clim_df <- data.frame( yhat_i = paste0('V',1:9),
                       clim_x = rowMeans( x_9 ) )

# preliminary simulations
yr_prelim <- data.frame( 
  alpha = rnorm(  it, 0,    1),
  beta  = rnorm(  it, 0,    1),
  y_sd  = rgamma( it, gam_p, rate = gam_p)
)

# calculate yhat from predictor variable clim_m
yhat_calc   <- function( a, b, clim_m) a + b * clim_m

# simulation data frame
yr_r_sim <- yhat_calc( yr_prelim$alpha, 
                       yr_prelim$beta, 
                       clim_m ) %>% 
  as.data.frame %>% 
  mutate( sim_i = 1:nrow(.) ) %>% 
  bind_cols( yr_prelim, . ) %>% 
  gather( yhat_i, yhat, V1:V9) %>% 
  mutate( y_sim = rnorm( nrow(.), yhat, y_sd ) ) %>% 
  arrange( sim_i ) %>% 
  left_join( clim_df )

# extract R2
extract_R2 <- function( ii ){
  
  print( ii ) 

  yr_r_sim %>% 
    subset( sim_i == ii ) %>% 
    glm( y_sim ~ clim_x, data = .) %>% 
    Dsquared()
    
}

# get those r2!
d2_lm <- sapply(1:1000, extract_R2)



# GAMMA PROCESS ----------------------------------------------------------------

# coefficients
fec_coef <- read.csv('results/prior_pc/fec_meanvar_coef.csv')$Estimate

# extreme variances for fecundity
fec_var  <- data.frame( x      = c(-3.5, -1.5, 0.5 ) ) %>% 
                mutate( y        = fec_coef[1] + fec_coef[2] * x ) %>% 
                # convert to natural scale
                mutate( mean     = exp(x),
                        variance = exp(y) )

design_df <- expand.grid( var_mean = fec_var$mean,
                          var_var  = fec_var$variance ) 

# order 
design_df <- list( design_df[7:9,], 
                   design_df[4:6,], 
                   design_df[1:3,] ) %>% bind_rows

rgamma( it, ((design_df$var_var[6]^2) / 0.1), (design_df$var_var[6] / 0.1) ) %>% sd

# d2 by mean variance relationship
d2_by_meanvar <- function( ii ){
          
  var_mean  <- design_df$var_mean[ii]
  var_var   <- design_df$var_var[ii]
  
  # preliminary simulations
  yr_prelim <- data.frame( 
    alpha = rnorm(  it, var_mean,     1),
    beta  = rnorm(  it, 0,     1),
    y_sd  = rgamma( it, ((var_var^2) / 1), (var_var / 1) )
  )
  
  # calculate yhat from predictor variable clim_m
  yhat_calc   <- function( a, b, clim_m) exp(a + b * clim_m)
  
  # simulation data frame
  yr_r_sim <- yhat_calc( yr_prelim$alpha, 
                         yr_prelim$beta, 
                         clim_m ) %>% 
    as.data.frame %>% 
    mutate( sim_i = 1:nrow(.) ) %>% 
    bind_cols( yr_prelim, . ) %>% 
    gather( yhat_i, yhat, V1:V9) %>% 
    mutate( y_sim = rgamma( nrow(.), shape = (yhat^2 / y_sd), rate = (yhat / y_sd) ) ) %>% 
    arrange( sim_i ) %>% 
    left_join( clim_df )
  
  # extract R2
  extract_Dsquared <- function( ii ){
    
    print( ii ) 
  
    tmp_df <- yr_r_sim %>% subset( sim_i == ii ) 
    
    tryCatch( mod <- glm(y_sim ~ clim_x, data = tmp_df,
                         family  = Gamma(link = "log"), maxit = 1000 ) %>% Dsquared(),
              error = function(cond){
                mod <- NULL 
              },
              warning = function(cond){
                return(NULL)
              })
    
    mod
    
  }
  
  sapply( 1:1000, extract_Dsquared )

}


d2_l   <- lapply( 1:nrow(design_df), d2_by_meanvar)
d2_df  <- d2_l %>% 
            bind_cols %>% 
            as.data.frame %>% 
            setNames( paste0('V',1:9) ) %>% 
            gather( sim_i, d2, V1:V9) 
  
d2_df %>% 
  group_by( sim_i ) %>% 
  summarise( na_n = is.na(d2) %>% sum )


ggplot(d2_df) +
  geom_histogram( aes(d2) ) +
  facet_wrap( ~ sim_i ) + 
  theme_minimal() +
  labs( x = 'Deviance explained',
        y = 'Count' ) +
  theme( axis.title = element_text( size = 20 ) ) + 
  ggsave( 'results/prior_pc/dev_explained_moments.tiff',
          width = 6.3, height = 6.3, compression = 'lzw')
  
  
# ML Gamma regression ------------------------------------------

design_df <- expand.grid( var_mean = fec_var$mean,
                          alpha    = c(0.01, 0.1, 1, 10) ) %>% 
                mutate( sim_i = paste0('V',1:12) ) %>% 
                mutate( descr = paste0('Mean = ',round(var_mean,2),
                                       '; alpha = ',round(alpha,2) ) )

# d2 by mean variance relationship
d2_by_meanvar <- function( ii ){
  
  var_mean  <- design_df$var_mean[ii]
  alpha     <- design_df$alpha[ii]
  
  # preliminary simulations
  yr_prelim <- data.frame( 
    alpha = rnorm(  it, log(var_mean),     1),
    beta  = rnorm(  it, 0,     1),
    y_sd  = rgamma( it, alpha, alpha )
  )
  
  # calculate yhat from predictor variable clim_m
  yhat_calc   <- function( a, b, clim_m) exp(a + b * clim_m)
  
  # simulation data frame
  yr_r_sim <- yhat_calc( yr_prelim$alpha, 
                         yr_prelim$beta, 
                         clim_m ) %>% 
    as.data.frame %>% 
    mutate( sim_i = 1:nrow(.) ) %>% 
    bind_cols( yr_prelim, . ) %>% 
    gather( yhat_i, yhat, V1:V9) %>% 
    mutate( y_sim = rgamma( nrow(.), shape = (yhat^2 / y_sd), rate = (yhat / y_sd) ) ) %>% 
    arrange( sim_i ) %>% 
    left_join( clim_df )
  
  
  # extract R2
  extract_Dsquared <- function( ii ){
    
    print( ii ) 
    
    tmp_df <- yr_r_sim %>% subset( sim_i == ii ) 
    
    tryCatch( mod <- glm(y_sim ~ clim_x, data = tmp_df,
                         family  = Gamma(link = "log"), maxit = 1000 ) %>% Dsquared(),
              error = function(cond){
                mod <- NULL 
              },
              warning = function(cond){
                return(NULL)
              })
    
    mod
    
  }
  
  sapply( 1:1000, extract_Dsquared )
  
}

# d2 for maximum likelihood models
d2_ml_l   <- lapply( 1:nrow(design_df), d2_by_meanvar)

d2_lm_df  <- d2_ml_l %>% 
  bind_cols %>% 
  as.data.frame %>% 
  setNames( paste0('V',1:12) ) %>% 
  gather( sim_i, d2, V1:V12) %>% 
  left_join(design_df)

d2_lm_df %>% 
  group_by( descr ) %>% 
  summarise( na_n = is.na(d2) %>% sum )

# Example
ggplot(d2_lm_df) +
  geom_histogram( aes(d2) ) +
  facet_wrap( ~ descr ) + 
  theme_minimal() +
  labs( x = 'Deviance explained',
        y = 'Count' ) +
  theme( axis.title = element_text( size = 15 ),
         strip.text = element_text( size = 8 ) ) + 
  ggsave( 'results/prior_pc/dev_explained_fec_ml.tiff',
          width = 6.3, height = 5.5, compression = 'lzw')


# Gamma Prior -------------------------------------------------


data.frame( 
  alpha_0.01 = rgamma( 1000, 0.01, 0.01 ),
  alpha_0.1  = rgamma( 1000, 0.1,  0.1 ),
  alpha_1    = rgamma( 1000, 1,    1 ),
  alpha_10   = rgamma( 1000, 10,   10 ) ) %>% 
  gather( alpha_val, value, alpha_0.01:alpha_10 ) %>% 
  ggplot() + 
  geom_histogram( aes( value ) ) + 
  facet_wrap( ~ alpha_val, scale = 'free' ) +
  theme_minimal() + 
  labs( x = expression(epsilon),
        y = 'Count' ) +
  theme( axis.title = element_text( size = 20 ),
         strip.text = element_text( size = 15 ) ) + 
  ggsave( 'results/prior_pc/epsilon_prior.tiff',
          width = 6.3, height = 5.5, compression = 'lzw')
