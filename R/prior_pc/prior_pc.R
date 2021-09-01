library(tidyverse)
library(gtools)
library(ggjoy)
library(ggridges)
library(extraDistr)
library(truncdist) # Not sure necessary to simulate half-cauchy...
library(mvtnorm)

# General parameters
it       <- 60000
gam_p    <- 1
alpha_sd <- 1
beta_sd  <- 1

# Simulate climatic predictors. 9 years, 12 months.
set.seed( 1776 + 4)
sig       <- diag(12)
x_9       <- MASS::mvrnorm(9,  rep(0,12), Sigma=sig) 



# NULL model R simulation ------------------------------------------------------

# Simulate parameters from prior
null_prior <- data.frame( 
  alpha = rnorm(  it, 0,     alpha_sd),
  y_sd  = rgamma( it, gam_p, gam_p)
)

# calculate yhats
yhat_calc <- function(alpha, y_sd){
  
  rnorm(9, alpha, y_sd) %>% t %>% as.data.frame
  
} 

# simulate predicted (yhat) and observed (y_sim) data
null_r_sim <- Map(yhat_calc, 
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
yhat_calc   <- function( a, b, clim_m) a + b * clim_m

# simulate predicted (yhat) and observed (y_sim) data
yr_r_sim <- yhat_calc( yr_prior$alpha, 
                       yr_prior$beta, 
                       clim_m ) %>% 
  as.data.frame %>% 
  bind_cols( yr_prior, . ) %>% 
  gather( yhat_i, yhat, V1:V9) %>% 
  mutate( y_sim = rnorm( nrow(.), yhat, y_sd ) )



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
yhat_calc   <- function( a, b, i, ante_x) a + b * ante_x[i,]

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
  mutate( y_sim  = rnorm( nrow(.), yhat, y_sd ) ) %>% 
  mutate( ante_1 = ante_1[1:nrow(.)] )



# SAM: R Simulation ------------------------------------------------------------

# Simplex matrix 
simplex_mat <- rdirichlet( it, rep(1,12) )

# Antecedents (I checked this!)
ante_x      <- x_9 %*% t(simplex_mat) %>% t 

# calculate yhat from antecedent
yhat_calc   <- function( a, b, i, ante_x) a + b * ante_x[i,]

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
                  mutate( y_sim = rnorm( nrow(.), yhat, y_sd ) ) %>% 
                  mutate( ante_1 = ante_1[1:nrow(.)] )



# Finnish Horseshoe ------------------------------------------------------------

# starting parameters
m0            <- 0.5
M             <- 12 # number of months/predictors
n             <- 9  # Number of data points (N in the formulas, to be sqrt)
slab_df       <- 25 # "v" in the formulas
slab_scale    <- 1  # "s" in the formulas

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
  alpha[i] + (x_9 %*% beta[i,]) %>% 
    as.numeric 
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
  mutate( y_sim = rnorm( nrow(.), yhat, y_sd ) )



# Evaluate R simulations -------------------------------------------------------
sim_df <- list( 
                'Null'  = mutate( null_r_sim, model = 'Null'), 
                'Year'  = mutate( yr_r_sim,   model = 'Year'), 
                'WMM'   = mutate( gaus_r_sim, model = 'Gaus'), 
                'SAM'   = mutate( sam_r_sim,  model = 'SAM'), 
                'Horse' = mutate( fh_r_sim,   model = 'Horse')
             ) %>% 
            # lapply( function(x) dplyr::select(x, model,yhat,y_sim, beta) ) %>% 
            bind_rows %>% 
            mutate( yhat_noalpha = yhat - alpha )

# x11()
# boxplot( yhat ~ model,  data = sim_df)
# boxplot( y_sim ~ model, data = sim_df)
# boxplot( beta ~ model, data = sim_df)


# y_hat figures priors
ggplot(sim_df) + 
  geom_histogram( aes(yhat) ) +
  facet_wrap( ~ model ) +
  xlim(-5,5) +
  ggsave( 'results/prior_pc/log_lambda_y_hat_original.tiff',
          width = 6.3, height = 4, compression = 'lzw' )

# y_sim figures priors
ggplot(sim_df) + 
  geom_histogram( aes(y_sim) ) +
  facet_wrap( ~ model ) +
  xlim(-5,5) +
  ggsave( 'results/prior_pc/log_lambda_y_sim_alpha_sd01.tiff',
          width = 6.3, height = 4, compression = 'lzw' )



min_max_ante_1_df <- sim_df %>% 
  subset( !(model %in% c('Null','Year')) ) %>% 
  group_by( model ) %>% 
  summarise( max_ante = max(ante_1, na.rm=T),
             min_ante = min(ante_1, na.rm=T),
             sd_ante  = sd(ante_1, na.rm=T) ) %>% 
  ungroup

sim_df %>% 
  subset( !is.na(alpha) ) %>% 
  gather( yhat_type, yhat_value, yhat, yhat_noalpha, alpha ) %>% 
  ggplot() + 
  geom_density_ridges( aes(x = yhat_value, y = yhat_type, color = yhat_type),
                       scale = 1) + 
  facet_wrap( ~ model ) + 
  theme_minimal() + 
  ggsave( 'results/prior_pc/y_hat_type_density_joyplot.tiff',
          width = 6.3, height = 4, compression = 'lzw' )

# Single monthly predictions
sim_df %>% 
  subset( !(model %in% c('Null','Year')) ) %>% 
  # subset( !is.na(ante_1) ) %>% 
  ggplot( ) + 
  geom_histogram( aes(ante_1) ) +
  facet_wrap( ~model ) +
  xlim(-4,4) +
  geom_vline( data = min_max_ante_1_df,
              aes( xintercept = min_ante),
              lty = 2) +
  geom_vline( data = min_max_ante_1_df,
              aes( xintercept = max_ante),
              lty = 2) +
  ggsave( 'results/prior_pc/beta_1.tiff',
          width = 6.3, height = 4, compression = 'lzw' )


# Overlapping density plots
sim_df %>% 
  subset( !(model %in% c('Null','Year')) ) %>% 
  # subset( !is.na(ante_1) ) %>% 
  ggplot( ) + 
  geom_density( aes( x = ante_1, color = model),
                alpha = 0.2, lwd= 3 ) +
  theme_minimal() +
  labs( x = 'Antecedent month 1',
        y = 'Density') + 
  ggsave( 'results/prior_pc/ante_1_density_overlay.tiff',
          width = 6.3, height = 4, compression = 'lzw' )

# Overlapping density plots
sim_df %>% 
  subset( !(model %in% c('Null','Year')) ) %>% 
  # subset( !is.na(ante_1) ) %>% 
  ggplot( ) + 
  geom_density_ridges( aes(x = ante_1, y = model, color = model),
                      scale = 1) + 
  # facet_wrap(~ model ) + 
  theme_minimal() +
  # geom_joy( aes( x = ante_1, y = model, color = model),
  #           alpha = 0.5, scale = 2  ) +
  # theme_joy() + 
  labs( x = 'Antecedent month 1',
        y = 'Density' ) + 
  ggsave( 'results/prior_pc/ante_1_density_joyplot.tiff',
          width = 6.3, height = 4, compression = 'lzw' )
  

# compare alpha and beta
sim_df %>% 
  subset( !(model %in% c('Null','Year')) ) %>% 
  gather( param, value, alpha, ante_1 ) %>% 
  ggplot( ) + 
  geom_density_ridges( aes(x = value, y = param, color = model),
                       scale = 1) + 
  facet_wrap( ~ model )  + 
  ggsave( 'results/prior_pc/alpha_vs_beta_density_joyplot.tiff',
          width = 6.3, height = 4, compression = 'lzw' )

aggregate( yhat  ~ model, FUN='sd', data = sim_df)
aggregate( y_sim ~ model, FUN='sd', data = sim_df)
aggregate( beta  ~ model, FUN='sd', data = sim_df)

# # set rstan options
# rstan_options( auto_write = TRUE )
# options( mc.cores = parallel::detectCores() )
# 
# # simulation parameters
# sim_pars <- list(
#   warmup = 10000, 
#   iter = 100000, 
#   thin = 4, 
#   chains = 3
# )
# 
# # extract vector of posterior output
# post_vec <- function( x ){
#   
#   x %>% 
#     rstan::extract() %>% 
#     .$y_sim %>% 
#     as.numeric 
#   
# }
# 
# post_yhat <- function( x ){
#   
#   x %>% 
#     rstan::extract() %>% 
#     .$yhat %>% 
#     as.numeric 
#   
# }
# 
# null_mod    <- stan_model( file = 'R/prior_pc/stan/normal_null.stan' )
# yr_mod      <- stan_model( file = 'R/prior_pc/stan/normal_yr.stan' )
# yr_fix_mod  <- stan_model( file = 'R/prior_pc/stan/normal_yr_fixed.stan' )
# gaus_mod    <- stan_model( file = 'R/prior_pc/stan/normal_gaus.stan' )
# sam_mod     <- stan_model( file = 'R/prior_pc/stan/normal_dirichlet.stan' )
# sam_fix_mod <- stan_model( file = 'R/prior_pc/stan/normal_dirichlet_fixed.stan' )
# fh_mod      <- stan_model( file = 'R/prior_pc/stan/normal_horse.stan' )
# 
# # simulate models
# dat_stan  <- list( n_time          = 9,
#                    n_lag           = 12,
#                    clim            = x_9,
#                    clim_means      = rowMeans( x_9 ),
#                    gamma_shape     = 0.01,
#                    alpha_dir       = 0.01,
#                    hs_df           = 1,   # variance of 
#                    hs_df_global    = 1,   # 
#                    hs_df_slab      = 25,   # slab degrees of freedom
#                    hs_scale_global = (4 / (36-4)) / 36, # global prior scale
#                    hs_scale_slab   = 2    # slab prior scale 
#                    )
# 
# # fit models and store them
# null_fit <- sampling( object = null_mod, 
#                       data = dat_stan,
#                       iter   = sim_pars$iter,
#                       warmup = sim_pars$warmup,
#                       thin   = sim_pars$thin,
#                       chains = sim_pars$chains 
#                        )
# saveRDS( object = null_fit, 
#          file   = 'results/prior_pc/normal_null_9.RDS' )
# rm(null_fit)
# 
# sim_pars$iter
# sim_pars$warmup
# sim_pars$thin
# 
# # Climate summary model
# yr_fit <- sampling( object = yr_mod, 
#                     data   = dat_stan,
#                     iter   = sim_pars$iter,
#                     warmup = sim_pars$warmup,
#                     thin   = sim_pars$thin,
#                     chains = sim_pars$chains 
#                     )
# saveRDS( object = yr_fit, 
#          file   = 'results/prior_pc/normal_yr_9.RDS' )
# rm(yr_fit)
# 
# # Climate summary model FIXED
# yr_fit_fix <- sampling( object = yr_fix_mod, 
#                         data   = dat_stan,
#                         iter   = sim_pars$iter,
#                         warmup = sim_pars$warmup,
#                         thin   = sim_pars$thin,
#                         chains = sim_pars$chains, 
#                         algorithm = 'Fixed_param' 
#                         )
# saveRDS( object = yr_fit_fix, 
#          file   = 'results/prior_pc/normal_yr_fixed_9.RDS' )
# rm(yr_fit_fix)
# 
# # Weighted mean model
# gaus_fit <- sampling( object = gaus_mod, 
#                       data   = dat_stan,
#                       iter   = sim_pars$iter,
#                       warmup = sim_pars$warmup,
#                       thin   = sim_pars$thin,
#                       chains = sim_pars$chains 
#                      )
# saveRDS( object = gaus_fit, 
#          file   = 'results/prior_pc/normal_gaus_9.RDS' )
# rm(gaus_fit)
# 
# # Stochastic antecedent model
# dat_stan$clim <- t(dat_stan$clim)
# sam_fit <- sampling( object = sam_mod, 
#                      data   = dat_stan,
#                      iter   = sim_pars$iter,
#                      warmup = sim_pars$warmup,
#                      thin   = sim_pars$thin,
#                      chains = sim_pars$chains 
#                      )
# saveRDS( object = sam_fit, 
#          file   = 'results/prior_pc/normal_dirichlet_9.RDS' )
# rm(sam_fit)
# 
# 
# # Stochastic antecedent model FIXED
# dat_stan$clim <- t(dat_stan$clim)
# sam_fit_fix <- sampling( object = sam_mod, 
#                      data   = dat_stan,
#                      iter   = sim_pars$iter,
#                      warmup = sim_pars$warmup,
#                      thin   = sim_pars$thin,
#                      chains = sim_pars$chains, 
#                      algorithm = 'Fixed_param'
#                     )
# saveRDS( object = sam_fit_fix, 
#          file   = 'results/prior_pc/normal_dirichlet_fixed_9.RDS' )
# rm(sam_fit_fix)
# 
# # Finnish Horseshoe model
# dat_stan$clim <- t(dat_stan$clim)
# fh_fit <- sampling( object = fh_mod, 
#                     data = dat_stan,
#                     iter   = sim_pars$iter,
#                     warmup = sim_pars$warmup,
#                     thin   = sim_pars$thin,
#                     chains = sim_pars$chains 
#                      )
# saveRDS( object = fh_fit, 
#          file   = 'results/prior_pc/normal_horse_9.RDS' )
# rm(fh_fit)
# 
# 
# 
# # analyze differences --------------------
# 
# 
# fit_null    <- readRDS( file = 'results/prior_pc/normal_null_9.RDS' )
# fit_yr      <- readRDS( file = 'results/prior_pc/normal_yr_9.RDS' )
# fit_yr_fix  <- readRDS( file = 'results/prior_pc/normal_yr_fixed_9.RDS' )
# fit_gaus    <- readRDS( file = 'results/prior_pc/normal_gaus_9.RDS' )
# fit_sam     <- readRDS( file = 'results/prior_pc/normal_dirichlet_9.RDS' )
# fit_sam_fix <- readRDS( file = 'results/prior_pc/normal_dirichlet_9.RDS' )
# fit_fh      <- readRDS( file = 'results/prior_pc/normal_horse_9.RDS' )
# 
# list( 
#       # 'null' = post_vec(fit_null), 
#       'year'    = post_vec(fit_yr),
#       'year_fx' = post_vec(fit_yr_fix),
#       'gaus'    = post_vec(fit_gaus), 
#       'sam'     = post_vec(fit_sam),
#       'sam_fx'  = post_vec(fit_sam),
#       'fh'      = post_vec(fit_fh) ) %>% 
#   bind_cols %>% 
#   boxplot
# 
# 
# 
# 
# 
# list( 'year' = post_yhat(fit_yr),
#       'gaus' = post_yhat(fit_gaus), 
#       'sam'  = post_yhat(fit_sam),
#       'fh'   = post_yhat(fit_fh) ) %>% 
#   bind_cols %>% 
#   boxplot
#   # apply( 2, sd)
# 
# list( 'null' = post_box(fit_null), 
#       'year' = post_box(fit_yr),
#       'gaus' = post_box(fit_gaus), 
#       'sam'  = post_box(fit_sam),
#       'fh'   = post_box(fit_fh) ) %>% 
#   bind_cols %>% 
#   apply( 2, sd )
# 
# 
# post_var <- function( x ){
#   
#   x %>% 
#     rstan::extract() %>% 
#     .$y_sim %>% 
#     apply(2, var, na.rm=T)
#   
# }
# 
# post_var(fit_null) 
# post_var(fit_yr) 
# post_var(fit_gaus)
# post_var(fit_sam) 
# post_var(fit_fh)
# 
# 
# list( 
#   # 'null' = post_vec(fit_null), 
#   'year'    = post_vec(fit_yr),
#   'year_fx' = post_vec(fit_yr_fix),
#   'year_R'  = yr_r_sim$y_sim,
#   'gaus'    = post_vec(fit_gaus), 
#   'sam'     = post_vec(fit_sam),
#   'sam_fx'  = post_vec(fit_sam),
#   'sam_R'   = sam_r_sim$y_sim,
#   'fh'      = post_vec(fit_fh) ) %>% 
#   bind_cols %>% 
#   boxplot
