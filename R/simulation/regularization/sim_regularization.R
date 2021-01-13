# Simulation to test performance of horseshoe and lasso regression
# 1. Simulation study for glmnet
# 2. Retrieve betas using glmnet, BRMS-horseshoe, and BRMS-lasso
# 3. Does glmnet is > NULL in cross-validation?
rm(list=ls())
source("code/format_data.R")
library(dplyr)
library(tidyr)
library(mgcv)
library(testthat)
library(rstan)
library(loo)
library(Rfast)
library(brms)
library(MASS)
library(glmnet)
library(ggthemes)
library(parallel)


# 1. Simulation study for glmnet ---------------------------

# MASTER FUNCTION to simulate normal data
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

# design matrix
design_df <- expand.grid( yr      = c(3:30),
                          sp_reps = c(1:5),
                          beta_p  = c(1:3) )

# beta values
beta_l   <- list( c( rep(0.05, 3), rep(0, 33) ),
                  c( rep(0.05, 6), rep(0, 30) ),
                  c( rep(0.05, 9), rep(0, 27) ) ) %>% 
              setNames( paste0( 1:3 ) )

# ridge regression
ridge_rmse <- function( ii ){
  
  # set true betas
  betas   <- beta_l[[design_df$beta_p[ii]]]
  
  # produce data
  norm_df <- sim_norm( nobs     = design_df$yr[ii], 
                       reps     = design_df$sp_reps[ii],
                       cfs      = betas, 
                       variance = 0.05, 
                       cov_mat  = diag(36) )
  
  # prepare data
  preds <- dplyr::select(norm_df, V1:V36) %>% as.matrix
  y_v   <- norm_df$y
  
  # glmnet crossvalidation
  cvfit <- cv.glmnet(x = preds, y = y_v, 
                     family = "gaussian",n=15,alpha=c(1),keep=T)
  
  # estimated beta values
  beta_est <- coef(cvfit, s = cvfit$lambda.1se)[,1][-1]
  
  # output 
  data.frame( yr      = design_df$yr[ii], 
              sp_reps = design_df$sp_reps[ii],
              beta_p  = design_df$beta_p[ii],
              b_t     = betas,
              b_e     = beta_est )
  
}

# replicate even the cross-validation
rep_crval <- function( tt ){
  
  lapply(1:nrow(design_df), ridge_rmse ) %>% 
    bind_rows %>% 
    mutate( iter = tt )
  
}

# detect cores
cores     <- detectCores()
# set up as many clusters as detected by 'detectCores()'
cluster   <- parallel::makePSOCKcluster(cores)

# attach packages that will be needed on each cluster
clusterEvalQ(cluster, list(library(dplyr),    library(tidyr), library(mgcv),
                           library(testthat), library(rstan), library(MASS),
                           library(glmnet)))
# attach objects that will be needed on each cluster
clusterExport(cluster, c('sim_norm', 'design_df',
                         'beta_l',   'ridge_rmse') )

# 100 repetitions, should take 40 min with 3 cores
ridge_l     <- parLapply(cluster, 1:100,  rep_crval)
ridge_df    <- do.call( rbind, ridge_l )

# # store the hard 40 minutes work!
# write.csv(ridge_df,
#           'results/simulations/regularization/ridge_simulation.csv',
#           row.names=F)

# read simulations (if already ran)
ridge_df    <- read.csv('results/simulations/regularization/ridge_simulation.csv') %>% 
                 # adjust beta to reflect ACTUAL number of parameters...
                 mutate( beta_p = replace(beta_p, beta_p == 3, 9),
                         beta_p = replace(beta_p, beta_p == 2, 6) ) %>% 
                 mutate( beta_p = replace(beta_p, beta_p == 1, 3) )

# root mean squared error
rmse        <- function(y, x) (y - x)^2 %>% mean %>% sqrt

# ONLY FOR TRUE BETA == 0: plot R2 versus sample sizes
ridge_df %>% 
  subset( b_t == 0 ) %>% 
  subset( yr > 5 ) %>% 
  group_by( yr, sp_reps, beta_p ) %>% 
  summarise( met = rmse(b_e, b_t) ) %>% 
  ungroup %>% 
  mutate( rep     = yr * sp_reps ) %>% 
  mutate( sp_reps = as.factor(sp_reps),
          beta_p  = as.factor(beta_p) ) %>% 
  ggplot( aes( x = yr,
               y = met ) ) + 
  geom_point( aes( color = sp_reps,
                   shape = beta_p ),
              alpha = 0.9 ) +
  scale_color_colorblind() +
  theme_minimal() + 
  ylim( 0, 0.055) +
  labs( x     = 'Reps (Years)',
        y     = 'RMSE: true zeros vs. beta zeros',
        color = 'Reps (Spatial)' ) + 
  ggsave( 'results/simulations/regularization/RMSE_b_zero_by_rep.tiff',
          width = 6.3, height = 5 )


# ONLY FOR TRUE BETA > 0: plot R2 versus sample sizes
ridge_df %>% 
  subset( yr > 5 ) %>% 
  subset( b_t != 0 ) %>% 
  group_by( yr, sp_reps, beta_p ) %>% 
  summarise( met = rmse(b_e, b_t) ) %>%
  ungroup %>% 
  mutate( rep     = yr * sp_reps ) %>% 
  mutate( sp_reps = as.factor(sp_reps),
          beta_p  = as.factor(beta_p) ) %>% 
  ggplot( aes( x = yr,
               y = met ) ) + 
  geom_point( aes( color = sp_reps,
                   shape = beta_p) ) +
  scale_color_colorblind() +
  theme_minimal() + 
  ylim( 0, 0.055) +
  labs( x     = 'Reps (Years)',
        y     = 'RMSE: true non-zeros vs. beta recovered',
        color = 'Reps (Spatial)' ) +
  ggsave( 'results/simulations/regularization/RMSE_b_nonzero_by_rep.tiff',
          width = 6.3, height = 5 )


# FOR ALL BETAS: plot R2 versus sample sizes
ridge_df %>%
  subset( yr > 5 ) %>% 
  group_by( yr, sp_reps, beta_p ) %>% 
  summarise( met = rmse(b_e, b_t) ) %>%
  ungroup %>% 
  mutate( beta_p = as.factor(beta_p) ) %>% 
  # subset( met > 0.5 )
  # subset( sp_reps == 1 ) %>% 
  mutate( rep     = yr * sp_reps ) %>% 
  mutate( sp_reps = as.factor(sp_reps) ) %>% 
  ggplot( aes( x = yr,
               y = met ) ) + 
  geom_point( aes( color = sp_reps,
                   shape = beta_p ) ) +
  scale_color_colorblind() +
  theme_minimal() + 
  ylim( 0, 0.055) +
  labs( x     = 'Reps (Years)',
        y     = 'RMSE: all true betas vs. beta recover',
        color = 'Reps (Spatial)' ) +
  ggsave( 'results/simulations/regularization/RMSE_b_all_by_rep.tiff',
          width = 6.3, height = 5 )



# 2. Retrieve betas using glmnet, BRMS-horseshoe, and BRMS-lasso -----

# estimate values using BRMS
est_brms <- function( ii ){
  
  seed_i <- ii + 5454
  
  set.seed( seed_i )
  
  # simulate dataset
  betas   <- c( 0.05, 0.05, 0, 0, -0.05, rep(0, 31) )
  norm_df <- sim_norm( nobs     = 30, 
                       reps     = 2,
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
                    cores=3,
                    control = list(adapt_delta = 0.99) )
  
  beta_horse <- brm(y~., data = norm_df, 
                    prior = c(prior(normal(0,2),
                                    class="Intercept"),
                              prior(horseshoe(), 
                                    class = "b")),
                    chains=2, 
                    iter=8000, 
                    cores=3,
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
              true   = c( 0, 0.05, 0.05, 0, 0, -0.05, rep(0, 31) )
  )
  
}

# run estimation
est_l   <- lapply(1:30, est_brms)
est_df  <- bind_rows( est_l ) %>% 
             mutate( order = rep( c(1:37), 30) )
beta_df <- data.frame( order = c(1:37),
                       betas = c( 0, 0.05, 0.05, 0, 0, -0.05, 
                                  rep(0, 31) ) )
true_df <- expand.grid( ds    = c( 1:30 ), 
                        order = c( 1:37 ) ) %>% 
             inner_join( beta_df ) %>% 
             arrange( order )  

# final file for plotting
r2_df <- inner_join( est_df, true_df) %>% 
            # subset( order %in% c(1,2,3,6) ) %>% 
            group_by( ds ) %>% 
            summarise( r2_lasso  = lm( lasso ~ betas) %>% 
                         summary %>% .$r.squared,
                       r2_horse  = lm( horse ~ betas) %>% 
                         summary %>% .$r.squared,
                       r2_glmnet = lm( glmnet ~ betas) %>% 
                         summary %>% .$r.squared 
            ) %>% 
            ungroup %>% 
            gather( type, r2, r2_lasso:r2_glmnet ) %>% 
            mutate( type = gsub('r2_','',type) )

# plot R2 
ggplot(r2_df) +
  geom_boxplot( aes(type, r2 ) ) +
  theme_minimal() +
  labs( x = 'Model type', 
        y = expression('R'^2) ) +
  ggsave('results/simulations/regularization/horse_lasso_glmnet_30yr_2rep.tiff',
         width = 6.3, height = 6.3, compression='lzw')


# # 3. Does glmnet is > NULL in cross-validation? --------------------
# 
# # MASTER FUNCTION to simulate normal data
# sim_norm <- function(nobs, reps, cfs, variance, cov_mat){
#   
#   nvar  <- length(cfs)  
#   beta  <- as.matrix(cfs)
#   X     <- mvrnorm( nobs, 
#                     mu    = rep(0, nvar), 
#                     Sigma = cov_mat )
#   
#   mu    <- ( X %*% beta )  # add noise if desired + rnorm(N, sd<-.01)
#   
#   sp_reps <- function( ii ){
#     
#     X %>% 
#       as.data.frame %>% 
#       mutate( y   = rnorm(nobs, mu, variance),
#               yr  = 1:nobs,
#               rep = ii )
#     
#   }
#   
#   lapply(1:reps, sp_reps) %>% bind_rows
#   
# }
# 
# # design matrix
# design_df <- expand.grid( yr      = c(20:30),
#                           sp_reps = c(1:5) )
# 
# # ridge regression
# cv_ridge_null <- function( ii ){
#   
#   betas   <- c( 0.05, 0.05, 0, 0, -0.05, rep(0, 31) )
#   
#   norm_df <- sim_norm( nobs     = design_df$yr[ii], 
#                        reps     = design_df$sp_reps[ii],
#                        cfs      = betas, 
#                        variance = 0.05, 
#                        cov_mat  = diag(36) )
#   
#   # leave-one-year-out cross-validation
#   cval <- function( yr_i ){
#     
#     train_df <- subset(norm_df, !(yr %in% yr_i) )
#     test_df  <- subset(norm_df, yr == yr_i )  
#     des_df   <- cbind(1, dplyr::select(test_df, V1:V36)) %>% 
#       as.matrix
#     
#     # prepare data
#     preds <- dplyr::select(train_df, V1:V36) %>% as.matrix
#     y_v   <- train_df$y
#     
#     # glmnet prediction
#     cvfit <- cv.glmnet(x = preds, y = y_v, 
#                        family = "gaussian",n=15,alpha=c(1),keep=T)
#     beta_est <- coef(cvfit, s = cvfit$lambda.1se)[,1]
#     pred_reg <- as.numeric(des_df %*% beta_est)
#     
#     # Null model prediction
#     pred_0   <- predict( lm( y_v ~ 1, data = train_df), 
#                          newdata = test_df)
#     
#     data.frame( pred_r = test_df$y - pred_reg,
#                 pred_0 = test_df$y - pred_0 )
#     
#   }
#   
#   pred_v <- lapply(1:design_df$yr[ii], cval) %>% 
#     bind_rows() %>% 
#     apply(2, mean ) %>% 
#     `^` (2) %>% 
#     sqrt
#   
#   data.frame( pred_0   = pred_v[2],
#               pred_reg = pred_v[1],
#               nobs     = design_df$yr[ii], 
#               reps     = design_df$sp_reps[ii] )
#   
# }
# 
# cia <- lapply(1:55, cv_ridge_null)
# bind_rows(cia) %>% 
#   mutate( sel = pred_0 < pred_reg )
# 
# 
# 
# # estimate values using BRMS
# est_glmnet <- function( ii ){
#   
#   seed_i <- ii + 5454
#   
#   set.seed( seed_i )
#   
#   # simulate dataset
#   betas   <- c( 0.05, 0.05, 0, 0, -0.05, rep(0, 31) )
#   norm_df <- sim_norm( nobs     = 30, 
#                        reps     = 2,
#                        cfs      = betas, 
#                        variance = 0.05, 
#                        cov_mat  = diag(36) )
#   
#   # prepare data
#   preds <- dplyr::select(norm_df, V1:V36) %>% as.matrix
#   y_v   <- norm_df$y
#   
#   # crossvalidation
#   cvfit <- cv.glmnet(x = preds, y = y_v, 
#                      family = "gaussian",n=15,alpha=c(1),keep=T)
#   
#   beta_est <- coef(cvfit, s = cvfit$lambda.1se)[,1]
#   
#   data.frame( ds      = ii,
#               order   = 1:37,
#               glmnet  = beta_est )
#   
# }
# 
# 
# glmnet_df <- lapply(1:30, est_glmnet) %>% bind_rows
# 
# 
# 
# # 
# retrieve_norm_beta <- function(ii){
#   
#   # make datasets reproducible
#   set.seed(ii)
#   
#   # simulate dataset
#   norm_df <- sim_norm( nobs     = 10, 
#                        reps     = 2,
#                        cfs      = c( 0.05, 0.05, 0, 0, -0.05, 0, 0, 0, 0, 0, 0, 0 ), 
#                        variance = 0.05, 
#                        cov_mat  = diag(12) )
#   
#   # fit model using horseshoe
#   norm_horse <- brm(y~., data=norm_df,
#                     prior=c(prior(normal(0,2),
#                                   class="Intercept"),
#                             prior( horseshoe (), 
#                                    class = "b")),
#                     chains=2,iter=4000, cores=3,
#                     control = list(adapt_delta = 0.999) )
#   
#   # fit model using lasso
#   norm_lasso <- brm(y~., data=norm_df,
#                     prior=c(prior(normal(0,2),
#                                   class="Intercept"),
#                             prior(lasso(), 
#                                   class = "b")),
#                     chains=2,iter=4000, cores=3,
#                     control = list(adapt_delta = 0.999) )
#   
# }
# 
# 
# 
# # simulate dataset
# norm_df <- sim_norm( nobs     = 3, 
#                      cfs      = c( 0.05, 0.05, 0, 0, -0.05, 0, 0, 0, 0, 0, 0, 0,
#                                    rep(0,24) ) , 
#                      variance = 0.05, 
#                      cov_mat  = diag(36) )
# 
# # replicate simulations ---------------------------------
# 
# rep_sim <- function( ii ){
#   
#   set.seed( ii )
#   
#   sim_norm( nobs     = 5, 
#             cfs      = c( 0.05, 0.05, 0, 0, -0.05, 0, 0, 0, 0, 0, 0, 0,
#                           rep(0,24) ) , 
#             variance = 0.05, 
#             cov_mat  = diag(36) )
#   
# }
# 
# 
# 
# 
# 
# # 
# sim_beta<- function(nobs, cfs, variance, cov_mat){
#   
#   nvar  <- length(cfs)  
#   beta  <- as.matrix(cfs)
#   X     <- mvrnorm( nobs, 
#                     mu    = rep(0, nvar), 
#                     Sigma = cov_mat )
#   
#   mu    <- plogis( X %*% beta )  # add noise if desired + rnorm(N, sd<-.01)
#   phi   <- variance
#   A     <- mu*phi
#   B     <- (1-mu)*phi
#   y     <- rbeta(nobs, A, B)
#   
#   X %>% 
#     as.data.frame %>% 
#     mutate( y = y ) %>% 
#     mutate( y = replace( y, y == 1, 0.99999) ) %>% 
#     mutate( y = replace( y, y == 0, 8.498758e-312) ) 
#   
# }
# 
# # produce the data
# beta_df <- sim_beta( nobs     = 10, 
#                      cfs      = c(0.05, 0.05, 0, 0, -0.05, 0, 0, 0, 0, 0, 0, 0), 
#                      variance = 0.1, 
#                      cov_mat = diag(12) )
# 
# # fit the model
# beta_ridg <- brm(y~., data=beta_df, 
#                  family='beta',
#                  prior=c(prior(normal(0,2),
#                                class="Intercept"),
#                          prior(lasso(), 
#                                class = "b")),
#                  chains=2,iter=4000, cores=3)
# 
# 
# # simulate
# sim_gamma <- function(nobs, cfs, variance, cov_mat){
#   
#   nvar  <- length(cfs)  
#   beta  <- as.matrix(cfs)
#   X     <- mvrnorm( nobs, 
#                     mu    = rep(0, nvar), 
#                     Sigma = cov_mat )
#   
#   mu    <- exp( X %*% beta )  # add noise if desired + rnorm(N, sd<-.01)
#   A     <- variance
#   B     <- variance / mu
#   y     <- rgamma(nobs, A, B)
#   
#   X %>% 
#     as.data.frame %>% 
#     mutate( y = y ) %>% 
#     mutate( y = replace(y, y == 0, 8.498758e-312) )
#   
# }
# 
# # produce the data
# gamma_df <- sim_gamma(nobs = 240, 
#                       cfs = c(3, 1.5, 0, 0, 2, 0, 0, 0), 
#                       variance = 10, cov_mat = diag(8) )
# 
# gamm_ridg <- brm(y~., data=gamma_df, 
#                  family=Gamma(link="log"),
#                  prior=c(prior(normal(0,2),
#                                class="Intercept"),
#                          prior(lasso(), 
#                                class = "b"),
#                          prior(gamma(0.01,0.01),
#                                class="shape")),
#                  chains=2,iter=4000, cores=3)
# 
# 
# 
# 
# fit_yr_old <- stan(
#   file = paste0("code/stan/",family,"_yr_old.stan"),
#   data = dat_stan,
#   pars = c('alpha', 'beta', 'sigma2',# 'y_sd', 
#            'yhat','log_lik'),
#   warmup = sim_pars$warmup,
#   iter = sim_pars$iter,
#   thin = sim_pars$thin,
#   chains = sim_pars$chains#,
#   # control = list(adapt_delta = 0.99)
# )
# 
# 
# fit_yr_new <- stan(
#   file = paste0("code/stan/",family,"_yr_new.stan"),
#   data = dat_stan,
#   pars = c('alpha', 'beta', 'sigma2',# 'y_sd', 
#            'yhat','log_lik'),
#   warmup = sim_pars$warmup,
#   iter = sim_pars$iter,
#   thin = sim_pars$thin,
#   chains = sim_pars$chains#,
#   # control = list(adapt_delta = 0.99)
# )
# 
# 
# fit_yr_dylan <- stan(
#   file = paste0("code/stan/",family,"_yr_dylan.stan"),
#   data = dat_stan,
#   pars = c('alpha', 'beta', 'phi',# 'y_sd', 
#            'yhat'),
#   warmup = sim_pars$warmup,
#   iter = sim_pars$iter,
#   thin = sim_pars$thin,
#   chains = sim_pars$chains#,
#   # control = list(adapt_delta = 0.99)
# )
# 
# 
# # plot results
# plot_df <- data.frame( x = dat_stan$clim_means,
#                        y = dat_stan$y )
# 
# x_seq <- seq( min(dat_stan$clim_means), 
#               max(dat_stan$clim_means), 
#               length.out=100 ) 
# 
# plot( y ~ x, data = plot_df )
# lines(x_seq,
#       exp( summary(fit_yr1)$summary[,'mean'][1] +
#              summary(fit_yr1)$summary[,'mean'][2] * x_seq ),
#       lwd=2, col = 'Black', lty = 2)
# lines(x_seq,
#       exp( summary(fit_yr_old)$summary[,'mean'][1] +
#              summary(fit_yr_old)$summary[,'mean'][2] * x_seq ),
#       lwd=2, col = 'Black', lty = 1)
# 
# summary(fit_yr)$summary[,'mean'][1:3]
# summary(fit_yr_old)$summary[,'mean'][1:3]
# summary(fit_yr_new)$summary[,'mean'][1:3]
# summary(fit_yr_dylan)$summary[,'mean'][1:3]
# 
# # glm model results
# mod <- glm(y ~ x, data = plot_df, family = Gamma(log) )
# lines(x_seq, 
#       exp( coef(mod)[1] + coef(mod)[2] * x_seq),
#       col='red',       lwd=2)
# 
# # 
# brm_mod <- brm(y~x, data=plot_df, family=Gamma(link="log"),
#                prior=c(prior(normal(0,2),class="Intercept"),
#                        prior(normal(0,2),class="b"),
#                        prior(gamma(0.01,0.01),class="shape")),
#                chains=2,iter=1000, cores=4)
# 
# brm_fixed <- brm_mod %>% summary %>% .$fixed %>% .[,1]
# lines(x_seq, 
#       exp( brm_fixed[1] +  brm_fixed[2] * x_seq),
#       col='green',       lwd=1)
# 
# 
# x <- rgamma(1000, 2.5, 1.5 )
# mean(x) 
# mean(x)
# 
# # test Gammareg ----------------------------------
# library(Gammareg)
# x1  <- runif(500, 0, 30)
# mui <- exp( -5 + 0.2*x1 )
# phi <- 0.5
# y   <- rgamma(500, shape=phi, scale= mui/phi)
# 
# 
# data_sim <- list(
#   
#   y = y,
#   clim_means = x1
#   
# )
