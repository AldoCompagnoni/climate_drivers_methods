library(tidyverse)
library(gtools)
library(ggjoy)
library(ggridges)
library(extraDistr)
library(truncdist) # Not sure necessary to simulate half-cauchy...
library(mvtnorm)

# General parameters
it       <- 60000
gam_p    <- 0.01
alpha_sd <- 0.5
beta_sd  <- 1

# Simulate climatic predictors. 9 years, 12 months.
set.seed( 1776 + 4)
sig       <- diag(12)
x_9       <- MASS::mvrnorm(9,  rep(0,12), Sigma=sig) 


# NULL model R simulation ------------------------------------------------------

# Simulate parameters from prior
null_prior <- data.frame( 
  alpha = rnorm(  it, 0,     alpha_sd) %>% exp,
  y_sd  = rgamma( it, gam_p, gam_p)
)

# calculate yhats
ysim_calc <- function(alpha, y_sd){
  
  rgamma( 9, shape = (alpha^2 / y_sd^2), 
             rate  = (alpha / y_sd^2) ) %>% 
    t %>% 
    as.data.frame 
  
}

# simulate predicted (yhat) and observed (y_sim) data
null_r_sim <- Map(ysim_calc, 
                  null_prior$alpha, 
                  null_prior$y_sd) %>% 
  bind_rows %>% 
  bind_cols( null_prior, . ) %>% 
  gather( yhat_i, y_sim, V1:V9) %>% 
  rename( yhat = alpha ) 



# YEAR (CLIMATE SUMMARY) MODEL simulation --------------------------------------

# replicate yearly anomalies by posterior sample size
clim_m   <- rowMeans( x_9 ) %>% 
              t %>% 
              as.data.frame %>% 
              slice(rep(1:n(), each = it))

# preliminary simulations
yr_prior <- data.frame( 
              alpha = rnorm(  it, 0,    alpha_sd),
              beta  = rnorm(  it, 0,    beta_sd),
              y_sd  = rgamma( it, gam_p, gam_p)
              )
              
# calculate yhat from predictor variable clim_m
yhat_calc   <- function( a, b, clim_m) exp(a + b * clim_m)

# simulate predicted (yhat) and observed (y_sim) data
yr_r_sim <- yhat_calc( yr_prior$alpha, 
                       yr_prior$beta, 
                       clim_m ) %>% 
  as.data.frame %>% 
  bind_cols( yr_prior, . ) %>% 
  gather( yhat_i, yhat, V1:V9) %>% 
  mutate( y_sim = rgamma( nrow(.), 
                          shape = (yhat^2 / y_sd^2), 
                          rate  = (yhat   / y_sd^2) ) 
          )



# GAUS/WMM: R simulation -------------------------------------------------------

# Preliminary simulations
gaus_prior <- data.frame( 
                  alpha   = rnorm(  it, 0,     alpha_sd),
                  beta    = rnorm(  it, 0,     beta_sd),
                  y_sd    = rgamma( it, gam_p, gam_p),
                  i       = 1:it,
                  sens_mu = rtrunc( it, 'norm', 0, 12,  mean = 6.5, sd = 12), 
                  sens_sd = rtrunc( it, 'norm', 0, Inf, mean = 0.5, sd = 12)
                )

# make antecedent
simpl_make <- function( ii ){
  
  sens_raw <- dnorm( 1:12, 
                     gaus_prior$sens_mu[ii], 
                     gaus_prior$sens_sd[ii] )
  sens     <- sens_raw / sum(sens_raw)
  
}
 
# make antecedent
simplex_mat <- lapply(1:nrow(gaus_prior), simpl_make) %>% 
  bind_cols %>% 
  t %>% 
  as.data.frame

# Antecedents (I checked this!)
ante_x      <- x_9 %*% t(simplex_mat) %>% t 

# "Antecedent 1"
ante_1      <- gaus_prior$beta * simplex_mat[,1]

# All antecedents
function(i) gaus_prior$beta * simplex_mat[,i]


# calculate yhat from antecedent
yhat_calc   <- function( a, b, i, ante_x) exp(a + b * ante_x[i,])

# simulate predicted (yhat) and observed (y_sim) data
gaus_r_sim <- yhat_calc( gaus_prior$alpha, 
                         gaus_prior$beta, 
                         gaus_prior$i, 
                         ante_x ) %>% 
  as.data.frame %>% 
  bind_cols( gaus_prior, . ) %>% 
  gather( yhat_i, yhat, V1:V9) %>% 
  # some yhat are NAs, because sens_sd too low, so sens = rep(0,12)
  subset( !is.na(yhat) ) %>% 
  mutate( y_sim = rgamma( nrow(.), 
                          shape = (yhat^2 / y_sd^2), 
                          rate  = (yhat   / y_sd^2) ) 
         ) %>% 
  mutate( ante_1 = ante_1[1:nrow(.)] )



# SAM: R Simulation ------------------------------------------------------------

# Simplex matrix 
simplex_mat <- rdirichlet( it, rep(1,12) )

# Antecedents (I checked this!)
ante_x      <- x_9 %*% t(simplex_mat) %>% t 

# calculate yhat from antecedent
yhat_calc   <- function( a, b, i, ante_x) exp( a + b * ante_x[i,] )

# SAM: R simulation
sam_prior   <- data.frame( 
                        alpha = rnorm(  it, 0,    alpha_sd),
                        beta  = rnorm(  it, 0,    beta_sd),
                        y_sd  = rgamma( it, gam_p, gam_p),
                        i     = 1:it
                      ) 
  
# "Antecedent 1"
ante_1      <- sam_prior$beta * simplex_mat[,1]

# simulate predicted (yhat) and observed (y_sim) data
sam_r_sim   <- yhat_calc( sam_prior$alpha, 
                          sam_prior$beta, 
                          sam_prior$i, 
                          ante_x ) %>% 
                  as.data.frame %>% 
                  bind_cols( sam_prior, . ) %>% 
                  gather( yhat_i, yhat, V1:V9) %>% 
                  mutate( y_sim = rgamma( nrow(.), 
                                          shape = (yhat^2 / y_sd^2), 
                                          rate  = (yhat   / y_sd^2) ) 
                  ) %>%
                  mutate( ante_1 = ante_1[1:nrow(.)] )



# Finnish Horseshoe ------------------------------------------------------------

# starting parameters
m0            <- 0.5
M             <- 12 # number of months/predictors
n             <- 9  # Number of data points (N in the formulas, to be sqrt)
slab_df       <- 25 # "v" in the formulas
slab_scale    <- 0.1  # "s" in the formulas

# Simulations
alpha         <- rnorm(  it, 0,     alpha_sd )
y_sd          <- rgamma( it, gam_p, gam_p )
beta_tilde    <- rmvnorm( it, rep(0,12), sigma = diag(12) )
tau_tilde     <- rtrunc( it, 'cauchy', a = 0)
c2_tilde      <- rinvgamma( it, slab_df/2, slab_df/2 )

# Transformed parameters
t0            <- (m0 / (M-m0)) / (y_sd / sqrt(n))  # Eq. 5g
tau           <- t0 * tau_tilde                    # Eq. 5f
c2            <- slab_scale^2 * c2_tilde           # Eq. 5e
lambda        <- lapply( paste0(1:12), function(i) rtrunc( it, 'cauchy', a = 0) ) %>% 
                    setNames( paste0(1:12) ) %>% 
                    bind_cols                      # Eq. 5d
lambda_tilde  <- sqrt( (c2 * lambda^2) / 
                       (c2 + (tau^2 * lambda^2)) ) # Eq. 5c
beta1         <- lambda_tilde * beta_tilde         # Eq. 5b
beta          <- sweep( beta1, 1, tau, FUN = '*') %>% 
                        as.matrix                  # Eq. 5b

# calculate yhats
ycalc <- function( i, beta, alpha ){
  (alpha[i] + (x_9 %*% beta[i,])) %>% 
    as.numeric %>% 
    sapply( exp )
} 

# simulate predicted (yhat) and observed (y_sim) data
fh_r_sim <- lapply(1:it, ycalc, beta, alpha ) %>% 
  setNames( paste0(1:it) ) %>% 
  bind_cols %>% 
  t %>%   
  as.data.frame %>%
  mutate( i    = 1:it,
          y_sd = y_sd ) %>% 
  mutate( ante_1 = beta[1:nrow(.),1],
          beta   = rowMeans(beta),
          alpha  = alpha ) %>% 
  gather( yhat_i, yhat, V1:V9) %>% 
  subset( !is.na(yhat) ) %>% 
  mutate( y_sim = rgamma( nrow(.), 
                          shape = (yhat^2 / y_sd^2), 
                          rate  = (yhat   / y_sd^2) ) 
          )



# Evaluate R simulations -------------------------------------------------------

# Actual data
all_vr_data <- read.csv('results/prior_pc/all_vr_data.csv') %>% 
                 subset( vr == 'fec' ) %>% 
                 rename( model = vr,
                         y_sim = value ) %>% 
                 mutate( model = 'Data' )

# Simulations
sim_df <- list( 
                'NM'    = mutate( null_r_sim, model = 'NM'), 
                'CSM'   = mutate( yr_r_sim,   model = 'CSM'), 
                'WMM'   = mutate( gaus_r_sim, model = 'WMM'), 
                'SAM'   = mutate( sam_r_sim,  model = 'SAM'), 
                'FHM'   = mutate( fh_r_sim,   model = 'FHM')
             ) %>% 
            # lapply( function(x) dplyr::select(x, model,yhat,y_sim, beta) ) %>% 
            bind_rows %>% 
            mutate( yhat_noalpha = yhat - alpha ) %>% 
            bind_rows( all_vr_data ) %>% 
            mutate( model = factor( model,
                                    levels = c('Data', 'NM',  'CSM',
                                               'WMM' , 'SAM', 'FHM') )
            )

# y_sim figures priors
ggplot(sim_df, aes(y_sim) ) + 
  geom_histogram( aes(y = ..density..) ) +
  facet_wrap( ~ model ) +
  xlim(0,10) +
  ylim(0,1) +
  theme_minimal() +
  labs( y = 'Kernel density estimate',
        x = 'Data (Simulated or observed)' ) +
  theme( strip.text = element_text( size = 15),
         axis.title  = element_text( size = 15) ) + 
  ggsave( 'results/prior_pc/fec_y_sim.tiff',
          width = 6.3, height = 4, compression = 'lzw' )
