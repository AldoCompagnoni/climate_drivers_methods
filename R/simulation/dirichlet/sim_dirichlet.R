library(tidyverse)
library(MASS)
library(rstan)
library(glmnet)
library(ggthemes)

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
design_df <- expand.grid( yr      = seq(5,60,by=5),
                          beta_p  = c(3,6),
                          ii      = 1:4 )

# beta values
beta_l   <- list( c( rep(0.05, 3), rep(0, 33) ),
                  c( rep(0.05, 6), rep(0, 30) ) ) %>% 
              setNames( paste0( 1:2 ) )


# Master function to simulate normal data --------------------------

# normal process
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


# Estimation ----------------------------------------------------------------

# simulation parameters
sim_pars <- list(
  warmup = 1000, 
  iter = 4000, 
  thin = 2, 
  chains = 3
)

# estimate beta
estimate_beta <- function( ii ){

  print( ii ) 
  
  betas   <- beta_l[[design_df$beta_p[ii]/3]]
  norm_df <- sim_norm( nobs     = design_df$yr[ii], 
                       reps     = 1,
                       cfs      = betas,
                       variance = 0.15, 
                       cov_mat  = diag(36) )
  
  # organize data into list to pass to stan
  dat_stan <- list(
    n_time  = design_df$yr[ii], 
    n_lag   = 36,
    y       = norm_df$y,
    clim    = dplyr::select(norm_df, V1:V36) %>% t %>% as.matrix
  )
  
  # simplex year t
  fit_simpl <- stan(
    file = paste0("R/stan/normal_dirichlet.stan"),
    data = dat_stan,
    pars = c('theta', 'alpha', 'beta', 'y_sd', 
             'yhat',  'log_lik'),
    warmup = sim_pars$warmup,
    iter = sim_pars$iter,
    thin = sim_pars$thin,
    chains = sim_pars$chains#,
    #control = list(adapt_delta = 0.999, stepsize = 0.001, max_treedepth = 20)
  )
  
  summary(fit_simpl)$summary[,'mean'][1:36]

}

simpl_fit_l <- lapply(1:nrow(design_df), estimate_beta)
simpl_fit_l1 <- lapply(49:nrow(design_df), estimate_beta)

cia   <- append(simpl_fit_l, simpl_fit_l1) 

# format output of the models
form <- function( ii ){
  
  data.frame( theta  = as.numeric(cia[ii][[1]]),
              x      = 1:36,
              ii     = design_df$ii[ii],
              yr     = design_df$yr[ii],
              beta_p = design_df$beta_p[ii] )
  
}

# true values
true_df      <- data.frame( x      = c(1:36, 1:36),
                            beta_p = c(rep(3,36), rep(6,36) ),
                            true_p = c( c( rep(1/3,3), rep(0,33) ),
                                        c( rep(1/6,6), rep(0,30) ) ) 
                            )

# format output
simpl_fit_df <- lapply( 1:nrow(design_df), form ) %>% 
                  bind_rows %>% 
                  left_join( true_df )

# root mean squared error
rmse        <- function(y, x) (y - x)^2 %>% mean %>% sqrt


write.csv( simpl_fit_df,
           )

simpl_fit_df %>% 
  subset( true_p == 0 ) %>%
  group_by( yr, beta_p, ii ) %>% 
  summarise( rmse = rmse(true_p, theta) ) %>% 
  ungroup %>% 
  ggplot() +
  geom_point( aes( yr, rmse,
                   color = as.factor(beta_p)) 
              ) +
  scale_color_colorblind()

simpl_fit_df %>% 
  subset( beta_p == 9) %>% 
  ggplot() +
  geom_point( aes(x,theta,
                 color=as.factor(yr),
                 )) +
  scale_color_colorblind()




# Old regularization


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
    control = list( adapt_delta = 0.99 ),
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

# root mean squared error
rmse        <- function(y, x) (y - x)^2 %>% mean %>% sqrt

