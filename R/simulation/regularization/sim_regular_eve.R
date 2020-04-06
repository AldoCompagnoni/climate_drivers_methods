library(tidyverse)
library(brms)
library(MASS)
library(glmnet)

# arguments from command line
args     <- commandArgs(trailingOnly = TRUE)

# climate variable
job_n    <- as.numeric(args[1])
out_dir  <- args[2]


# Design of simulation --------------------------------------

# design matrix
design_df <- expand.grid( yr      = c(6:30),
                          sp_reps = c(1:5),
                          beta_p  = c(3,6,9),
                          rep     = 1:30 )

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
  
  
  # fit the regularized Bayes models--------------------------
  beta_lasso <- brm(y~., data = norm_df, 
                    prior = c(prior(normal(0,2),
                                    class="Intercept"),
                              prior(lasso(), 
                                    class = "b")),
                    chains=2, 
                    iter=8000, 
                    cores=4,
                    control = list(adapt_delta = 0.99) )
  
  beta_horse <- brm(y~., data = norm_df, 
                    prior = c(prior(normal(0,2),
                                    class="Intercept"),
                              prior(horseshoe(), 
                                    class = "b")),
                    chains=2, 
                    iter=8000, 
                    cores=4,
                    control = list(adapt_delta = 0.99) )
  
  # fit the glmnet model
  preds <- dplyr::select(norm_df, V1:V36) %>% as.matrix
  y_v   <- norm_df$y
  
  # crossvalidation
  cvfit <- cv.glmnet(x = preds, y = y_v, 
                     family = "gaussian",
                     n=15,alpha=c(1),keep=T)
  
  beta_net <- coef(cvfit, s = cvfit$lambda.1se)[,1]
  
  data.frame( ds     = ii,
              lasso  = summary(beta_lasso)$fixed[,'Estimate'],
              horse  = summary(beta_horse)$fixed[,'Estimate'],
              glmnet = beta_net,
              true   = c(0, beta_l[[design_df$beta_p[ii]/3]]) )
  
}

# simulation output
sim_df <- est_brms( job_n, design_df, beta_l )

# store the job!
write.csv(sim_df, 
          paste0(out_dir,'_sim_regular_',job_n,'.csv'), 
          row.names=F)