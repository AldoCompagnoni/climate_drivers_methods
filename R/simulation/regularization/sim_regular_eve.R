library(tidyverse)
# library(brms)
library(MASS)
library(rstan)
library(glmnet)

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
                          sl_df   = c(4,25),
                          glob_sc = c(1,
                                     (3  / (36-3))  / sqrt(60),
                                     (10 / (36-10)) / sqrt(60)),
                          rep     = 1:20 )

# beta values
beta_l   <- list( c( rep(0.05, 3), rep(0, 33) ),
                  c( rep(0.05, 6), rep(0, 30) ),
                  c( rep(0.05, 9), rep(0, 27) ) ) %>% 
              setNames( paste0( 1:3 ) )


# Master function to simulate normal data --------------------------
sim_norm<- function(nobs, reps, cfs, variance, cov_mat){
  
  # numer of predictors, their effect size, and their actual sim. vals.
  nvar  <- length(cfs)  
  beta  <- as.matrix(cfs)
  X     <- mvrnorm( nobs, 
                    mu    = rep(0, nvar), 
                    Sigma = cov_mat )
  
  # mean predicted values
  mu    <- ( X %*% beta )  # add noise if desired + rnorm(N, sd<-.01)
  
  # produce response variable for each "spatial" replicate
  sp_reps <- function( ii ){
    
    X %>% 
      as.data.frame %>% 
      mutate( y = rnorm(nobs, mu, variance) )
    
  }
  
  lapply(1:reps, sp_reps) %>% bind_rows
  
}

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
  norm_df <- sim_norm( nobs     = design_df$yr[ii], 
                       reps     = design_df$sp_reps[ii],
                       cfs      = betas,
                       variance = 0.05, 
                       cov_mat  = diag(36) )
  
  gamm_df <- sim_gamm( nobs     = design_df$yr[ii], 
                       reps     = design_df$sp_reps[ii],
                       cfs      = betas,
                       variance = 0.05, 
                       cov_mat  = diag(36) )
  
  beta_df <- sim_beta( nobs     = design_df$yr[ii], 
                       reps     = design_df$sp_reps[ii],
                       cfs      = betas,
                       variance = 0.05, 
                       cov_mat  = diag(36) )
  
  
  # horseshoe "from scratch", brms parameterization --------------------------
  
  # load data for horseshoe from brms
  data_in <- list( y = norm_df$y,
                   X = dplyr::select(norm_df,-y) %>% as.matrix,
                   N = norm_df$y %>% length,
                   K = dplyr::select(norm_df,-y) %>% ncol,
                   hs_df           = 1,   # variance of 
                   hs_df_global    = 1,   # 
                   hs_df_slab      = design_df$sl_df[ii],   # slab degrees of freedom
                   hs_scale_global = design_df$glob_sc[ii], # global prior scale
                   hs_scale_slab   = 2    # slab prior scale
  )
  
  mod          <- readRDS( 'horseshoe_brms_gauss.RDS' )
  
  # init_t <- Sys.time()
  beta_scratch <- sampling(
    object = mod,
    # file = "R/simulation/regularization/horseshoe_brms_gauss.stan",
    data = data_in,
    pars = c('Intercept','b'),
    warmup = 4000,
    iter   = 8000,
    thin   = 1,
    chains = 2,
    control = list(adapt_delta = 0.99),
    seed    = seed_i
  )
  
  
  # # horseshoe Piironen J. & Vehtari model --------------------------
  # 
  # # load data for horseshoe from brms
  # data_in <- list( y = norm_df$y,
  #                  x = dplyr::select(norm_df,-y) %>% as.matrix,
  #                  n = norm_df$y %>% length,
  #                  d = dplyr::select(norm_df,-y) %>% ncol,
  #                  scale_icept   = 2,  # prior variance for intercept
  #                  nu_local      = 1,  # local degrees of freedom
  #                  nu_global     = 1,  # global degrees of freedom
  #                  slab_df       = 25,  # slab degrees of freedom
  #                  scale_global  = design_df$glob_sc[ii],  # global prior scale
  #                  slab_scale    = 2   # slab prior scale
  # )
  # 
  # 
  # # mod_Piironen <- readRDS( 'R/simulation/regularization/horseshoe_regularized.rds')
  # 
  # # init_t <- Sys.time()
  # beta_Piironen <- stan(
  #   # mod_Piironen,
  #   file    = "R/simulation/regularization/horseshoe_regularized.stan",
  #   data    = data_in,
  #   pars    = c('Intercept','b'),
  #   warmup  = 4000,
  #   iter    = 8000,
  #   thin    = 1,
  #   chains  = 2,
  #   control = list(adapt_delta = 0.99),
  #   seed    = seed_i
  # )
  # 
  # 
  # # horseshoe from Betancourt's blog --------------------------
  # 
  # # load data for horseshoe from brms
  # data_in <- list( y = norm_df$y,
  #                  X = dplyr::select(norm_df,-y) %>% as.matrix %>% t,
  #                  N = norm_df$y %>% length,
  #                  M = dplyr::select(norm_df,-y) %>% ncol
  # )
  # 
  # # mod_betan <- readRDS("R/simulation/regularization/horseshoe_betancourt.rds")
  # 
  # # init_t <- Sys.time()
  # beta_betan <- stan(
  #   # mod_betan,
  #   file = "R/simulation/regularization/horseshoe_betancourt.stan",
  #   data = data_in,
  #   pars = c('alpha','beta'),
  #   warmup = 4000,
  #   iter   = 8000,
  #   thin   = 1,
  #   chains = 2,
  #   control = list(adapt_delta = 0.99),
  #   seed    = seed_i
  # )
  
  # Fit horseshoe model with BRMS --------------------------
  gaus_horse <- brm(y~., data = norm_df,
                    prior = c(prior(normal(0,2),
                                    class="Intercept"),
                              prior(horseshoe( df_slab      = design_df$sl_df[ii],
                                               scale_global = design_df$glob_sc[ii]),
                                    class = "b")),
                    chains  = 2,
                    iter    = 8000,
                    cores   = 2,
                    seed    = seed_i,
                    control = list(adapt_delta = 0.99) )
  
  # fit a gamma process
  gamm_horse <- brm(y~., data = gamm_df,
                    family=Gamma(link="log"),
                    prior = c(prior(normal(0,2),
                                    class="Intercept"),
                              prior(horseshoe( df_slab      = design_df$sl_df[ii],
                                               scale_global = design_df$glob_sc[ii]),
                                    class = "b")),
                    chains  = 2,
                    iter    = 8000,
                    cores   = 2,
                    seed    = seed_i,
                    control = list(adapt_delta = 0.99) )
  
  # fit a beta process
  beta_horse <- brm( y~., data = beta_df,
                     family=Beta(link="logit"),
                     prior = c(prior(normal(0,2),
                                     class="Intercept"),
                               prior(horseshoe( df_slab      = design_df$sl_df[ii],
                                                scale_global = design_df$glob_sc[ii]),
                                     class = "b")),
                     chains  = 2,
                     iter    = 8000,
                     cores   = 2,
                     seed    = seed_i,
                     control = list(adapt_delta = 0.99) )
  
  # fit the glmnet model
  preds <- dplyr::select(norm_df, V1:V36) %>% as.matrix
  y_v   <- norm_df$y
  
  # crossvalidation
  cvfit <- cv.glmnet(x = preds, y = y_v, 
                     family = "gaussian",
                     n=15,alpha=c(1),keep=T)
  
  beta_net <- coef(cvfit, s = cvfit$lambda.1se)[,1]
  
  # print( ii )
  # put it all out
  data.frame( ds       = ii,
              rep      = design_df$rep[ii],
              scrat    = summary(beta_scratch)$summary[,'mean'][1:37],
              # finn     = summary(beta_Piironen)$summary[,'mean'][1:37],
              # betan    = summary(beta_betan)$summary[,'mean'][1:37],
              # brms     = summary(beta_horse)$fixed[,'Estimate'],
              scrat_sd = summary(beta_scratch)$summary[,'sd'][1:37],
              # finn_sd  = summary(beta_Piironen)$summary[,'sd'][1:37],
              # betan_sd = summary(beta_betan)$summary[,'sd'][1:37],
              # brms_sd  = summary(beta_horse$fit)$summary[,'sd'][1:37],
              glmnet   = beta_net,
              true     = c(0, beta_l[[design_df$beta_p[ii]/3]]),
              beta_p   = design_df$beta_p[ii],
              sl_df    = design_df$sl_df[ii],
              glob_sc  = design_df$glob_sc[ii] )
  
}

# simulation output
# sim_l  <- lapply(1:3, est_brms, design_df, beta_l )
sim_df <- est_brms( job_n, design_df, beta_l )

# store the job!
write.csv( sim_df, 
           paste0(out_dir,'_sim_regular_',job_n,'.csv'), 
           row.names=F )