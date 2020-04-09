library(tidyverse)
# library(brms)
library(MASS)
library(rstan)
library(glmnet)
options(stringsAsFactors = F)

# set rstan options
rstan_options( auto_write = TRUE )
options( mc.cores = parallel::detectCores() )

# arguments from command line
args     <- commandArgs(trailingOnly = TRUE)

# climate variable
job_n    <- as.numeric(args[1])
out_dir  <- args[2]


# Design of simulation --------------------------------------


# design matrix
design_df <- expand.grid( yr      = 30,
                          sp_reps = 2,
                          beta_p  = 3,
                          sl_df   = c(25),
                          glob_sc = (10 / (36-10)) / sqrt(60),
                          rep     = 1:20,
                          mod     = c('gamma','beta'),
                          stringsAsFactors = F )

# beta values
beta_l   <- list( c( rep(0.05, 3), rep(0, 33) ),
                  c( rep(0.05, 6), rep(0, 30) ),
                  c( rep(0.05, 9), rep(0, 27) ) ) %>% 
              setNames( paste0( 1:3 ) )


# Master functions to simulate gamma/beta data --------------------------

# simualate gamma process
sim_gamm<- function(nobs, reps, cfs, variance, cov_mat){
  
  # numer of predictors, their effect size, and their actual sim. vals.
  nvar  <- length(cfs)  
  beta  <- as.matrix(cfs)
  X     <- mvrnorm( nobs, 
                    mu    = rep(0, nvar), 
                    Sigma = cov_mat )
  
  # mean predicted values
  mu    <- exp( X %*% beta )  # add noise if desired + rnorm(N, sd<-.01)
  sig   <- variance^2
  
  # gamma moment matching
  A <- mu^2 / sig
  B <- mu   / sig
  
  # produce response variable for each "spatial" replicate
  sp_reps <- function( ii ){
    
    X %>% 
      as.data.frame %>% 
      mutate( y = rgamma(nobs, A, B) )
    
  }
  
  lapply(1:reps, sp_reps) %>% bind_rows
  
}


# simualate beta process
sim_beta<- function(nobs, reps, cfs, variance, cov_mat){
  
  # numer of predictors, their effect size, and their actual sim. vals.
  nvar  <- length(cfs)  
  beta  <- as.matrix(cfs)
  X     <- mvrnorm( nobs, 
                    mu    = rep(0, nvar), 
                    Sigma = cov_mat )
  
  # mean predicted values
  mu    <- c(  X %*% beta )  # add noise if desired + rnorm(N, sd<-.01)
  
  # beta moment matching (maximum likelihood?)
  mu    <- boot::inv.logit( mu )  # add noise if desired + rnorm(N, sd<-.01)
  phi   <- 10
  A     <- mu*phi
  B     <- (1-mu)*phi
  
  # produce response variable for each "spatial" replicate
  sp_reps <- function( ii ){
    
    X %>% 
      as.data.frame %>% 
      mutate( y = rbeta(nobs, A, B) )
    
  }
  
  lapply(1:reps, sp_reps) %>% bind_rows
  
}



# estimate values using BRMS
est_brms <- function( ii, design_df, beta_l ){
  
  seed_i <- ii + 5454
  set.seed( seed_i )
  
  # simulate dataset
  betas   <- beta_l[[design_df$beta_p[ii]/3]]
  
  if( design_df$mod[ii] == 'gamma' ){
    
    # generate data
    data_df <- sim_gamm( nobs     = design_df$yr[ii], 
                         reps     = design_df$sp_reps[ii],
                         cfs      = betas,
                         variance = 0.05, 
                         cov_mat  = diag(36) )
    
    # load data for horseshoe from brms
    data_in <- list( y = data_df$y,
                     X = dplyr::select(data_df,-y) %>% as.matrix,
                     N = data_df$y %>% length,
                     K = dplyr::select(data_df,-y) %>% ncol,
                     hs_df           = 1,   # variance of 
                     hs_df_global    = 1,   # 
                     hs_df_slab      = design_df$sl_df[ii],   # slab degrees of freedom
                     hs_scale_global = design_df$glob_sc[ii], # global prior scale
                     hs_scale_slab   = 2    # slab prior scale
    )
    
    # fit model
    mod          <- readRDS( 'gamma_horse_brms.RDS' )
    fit_mod      <- sampling(
      object = mod,
      data = data_in,
      pars = c('Intercept','b'),
      warmup = 4000,
      iter   = 8000,
      thin   = 1,
      chains = 2,
      control = list(adapt_delta = 0.99),
      seed    = seed_i
    )
    
  }
  if( design_df$mod[ii] == 'beta' ){
    
    data_df <- sim_beta( nobs     = design_df$yr[ii], 
                         reps     = design_df$sp_reps[ii],
                         cfs      = betas,
                         variance = 0.05, 
                         cov_mat  = diag(36) )
    # fit model
    mod          <- readRDS( 'beta_horse_brms.RDS' )
    fit_mod      <- sampling(
      object = mod,
      data = data_in,
      pars = c('Intercept','b'),
      warmup = 4000,
      iter   = 8000,
      thin   = 1,
      chains = 2,
      control = list(adapt_delta = 0.99),
      seed    = seed_i
    )
    
  }
  
  # put it all out
  data.frame( ds       = ii,
              rep      = design_df$rep[ii],
              scrat    = summary(fit_mod)$summary[,'mean'][1:37],
              scrat_sd = summary(fit_mod)$summary[,'sd'][1:37],
              true     = c(0, beta_l[[design_df$beta_p[ii]/3]]),
              beta_p   = design_df$beta_p[ii],
              sl_df    = design_df$sl_df[ii],
              glob_sc  = design_df$glob_sc[ii] )
  
}

# simulation output
sim_df <- est_brms( job_n, design_df, beta_l )

# store the job!
write.csv( sim_df, 
           paste0(out_dir,'_sim_regular_',job_n,'.csv'), 
           row.names=F )