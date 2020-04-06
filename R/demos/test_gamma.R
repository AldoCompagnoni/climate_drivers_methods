
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

# set rstan options
rstan_options( auto_write = TRUE )
options( mc.cores = parallel::detectCores() )

# climate predictor, response, months back, max. number of knots
response  <- "surv"
clim_var  <- "precip"
m_back    <- 36    
st_dev    <- FALSE

# read data -----------------------------------------------------------------------------------------
lam       <- read.csv("data/all_demog_updt.csv", stringsAsFactors = F)
# m_info    <- read.csv("C:/cloud/Dropbox/sApropos/MatrixEndMonth_information.csv", stringsAsFactors = F)
clim      <- data.table::fread(paste0('data/',clim_var,"_chelsa_prism_hays_2014.csv"),  stringsAsFactors = F)
spp       <- lam$SpeciesAuthor %>% unique

# format data --------------------------------------------------------------------------------------

# set up model "family" based on response
if( response == "surv" | response == "grow" )             family = "beta" 
if( response == "fec" )                                   family = "gamma"
if( grepl("PreRep", response) | grepl("Rep", response) )  family = "beta"
if( response == "rho" | response == "react_fsa" )         family = "gamma"
if( response == "log_lambda")                             family = "normal"

expp_beta     <- 20

# set species (I pick Sphaeraclea_coccinea)
ii            <- 37
spp_name      <- spp[ii]

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
  clim_yr = list( rowMeans(mod_data$climate[, 1:12]),
                  rowMeans(mod_data$climate[,13:24]),
                  rowMeans(mod_data$climate[,25:36]) ) %>% do.call(rbind,.),
  M       = 12,    # number of months in a year
  K       = ncol(mod_data$climate) / 12,
  S       = mod_data$resp$population %>% unique %>% length,
  site_i  = mod_data$resp$population %>% as.factor %>% as.numeric,
  expp_beta = expp_beta
)

dat_stan$clim_means  <- rowMeans(mod_data$climate[,1:12 ])

# simulation parameters
sim_pars <- list(
  warmup = 1000, 
  iter = 4000, 
  thin = 4, 
  chains = 1
)

mod_df <- mod_data$climate %>% 
             mutate( y = mod_data$resp[,response] ) 
             

# simulate dataset --------------------------------------------------------
library(MASS)
set.seed(31102018)

# cia
sim_beta<- function(nobs, cfs, variance, cov_mat){
  
  nvar  <- length(cfs)  
  beta  <- as.matrix(cfs)
  X     <- mvrnorm( nobs, 
                    mu    = rep(0, nvar), 
                    Sigma = cov_mat )
  
  mu    <- plogis( X %*% beta )  # add noise if desired + rnorm(N, sd<-.01)
  phi   <- variance
  A     <- mu*phi
  B     <- (1-mu)*phi
  y     <- rbeta(nobs, A, B)
  
  X %>% 
    as.data.frame %>% 
    mutate( y = y ) %>% 
    mutate( y = replace( y, y == 1, 0.99999) ) %>% 
    mutate( y = replace( y, y == 0, 8.498758e-312) ) 
  
  # return( list(x=X, y=y, coefficients=cfs, var=variance, covX=cov.mat) )
}

# produce the data
beta_df <- sim_beta(nobs = 240, 
                cfs = c(3, 1.5, 0, 0, 2, 0, 0, 0), 
                variance = 10, cov_mat = diag(8) )

# fit the model
beta_ridg <- brm(y~., data=beta_df, 
                 family='beta',
                 prior=c(prior(normal(0,2),
                               class="Intercept"),
                         prior(lasso(), 
                               class = "b")),
                 chains=2,iter=4000, cores=3)


fit_beta_ridge <- stan(
  file = paste0("code/stan/",'beta',"_ridge.stan"),
  data = dat_stan,
  pars = c('alpha', 'beta', 'y_sd',
           'sigma2','lambda', 'yhat'),
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains#,
  # control = list(adapt_delta = 0.99)
)



# simulate
sim_gamma <- function(nobs, cfs, variance, cov_mat){
  
  nvar  <- length(cfs)  
  beta  <- as.matrix(cfs)
  X     <- mvrnorm( nobs, 
                    mu    = rep(0, nvar), 
                    Sigma = cov_mat )
  
  mu    <- exp( X %*% beta )  # add noise if desired + rnorm(N, sd<-.01)
  A     <- variance
  B     <- variance / mu
  y     <- rgamma(nobs, A, B)
  
  X %>% 
    as.data.frame %>% 
    mutate( y = y ) %>% 
    mutate( y = replace(y, y == 0, 8.498758e-312) )
  
}

# produce the data
gamma_df <- sim_gamma(nobs = 240, 
                     cfs = c(3, 1.5, 0, 0, 2, 0, 0, 0), 
                     variance = 10, cov_mat = diag(8) )

gamm_ridg <- brm(y~., data=gamma_df, 
                 family=Gamma(link="log"),
                prior=c(prior(normal(0,2),
                              class="Intercept"),
                        prior(lasso(), 
                              class = "b"),
                        prior(gamma(0.01,0.01),
                              class="shape")),
                chains=2,iter=4000, cores=3)


# organize data into list to pass to stan
dat_stan <- list(
  n_time  = dplyr::select(gamma_df, -y) %>% nrow,
  n_lag   = dplyr::select(gamma_df, -y) %>% ncol,
  y       = gamma_df$y,
  clim    = dplyr::select(gamma_df, -y),
  sigma   = 1
)

# simulation parameters
sim_pars <- list(
  warmup = 1000, 
  iter = 2000, 
  thin = 2, 
  chains = 3
)

# My gamma ridge regression "from scratch" (from moment matching)
family <- 'gamma'
fit_ridge <- stan(
  file = paste0("code/stan/",family,"_ridge2.stan"),
  data = dat_stan,
  pars = c('alpha', 'beta', #'y_sd'
           'lambda', 'yhat'),
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains#,
  # control = list(adapt_delta = 0.99)
)


# simulate normal data
sim_norm<- function(nobs, cfs, variance, cov_mat){
  
  nvar  <- length(cfs)  
  beta  <- as.matrix(cfs)
  X     <- mvrnorm( nobs, 
                    mu    = rep(0, nvar), 
                    Sigma = cov_mat )
  
  mu    <- ( X %*% beta )  # add noise if desired + rnorm(N, sd<-.01)
  y     <- rnorm(nobs, mu, variance)
  
  X %>% 
    as.data.frame %>% 
    mutate( y = y )
  
}

# produce the data
norm_df <- sim_norm(nobs = 240, 
                    cfs = c(3, 1.5, 0, 0, 2, 0, 0, 0), 
                    variance = 1, cov_mat = diag(8) )

norm_ridg <- brm(y~., data=norm_df,
                 prior=c(prior(normal(0,2),
                               class="Intercept"),
                         prior(lasso(), 
                               class = "b")),
                 chains=2,iter=4000, cores=3)


fit_yr_old <- stan(
  file = paste0("code/stan/",family,"_yr_old.stan"),
  data = dat_stan,
  pars = c('alpha', 'beta', 'sigma2',# 'y_sd', 
           'yhat','log_lik'),
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains#,
  # control = list(adapt_delta = 0.99)
)


fit_yr_new <- stan(
  file = paste0("code/stan/",family,"_yr_new.stan"),
  data = dat_stan,
  pars = c('alpha', 'beta', 'sigma2',# 'y_sd', 
           'yhat','log_lik'),
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains#,
  # control = list(adapt_delta = 0.99)
)


fit_yr_dylan <- stan(
  file = paste0("code/stan/",family,"_yr_dylan.stan"),
  data = dat_stan,
  pars = c('alpha', 'beta', 'phi',# 'y_sd', 
           'yhat'),
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains#,
  # control = list(adapt_delta = 0.99)
)


# plot results
plot_df <- data.frame( x = dat_stan$clim_means,
                       y = dat_stan$y )

x_seq <- seq( min(dat_stan$clim_means), 
              max(dat_stan$clim_means), 
              length.out=100 ) 

plot( y ~ x, data = plot_df )
lines(x_seq,
      exp( summary(fit_yr1)$summary[,'mean'][1] +
           summary(fit_yr1)$summary[,'mean'][2] * x_seq ),
      lwd=2, col = 'Black', lty = 2)
lines(x_seq,
      exp( summary(fit_yr_old)$summary[,'mean'][1] +
           summary(fit_yr_old)$summary[,'mean'][2] * x_seq ),
      lwd=2, col = 'Black', lty = 1)

summary(fit_yr)$summary[,'mean'][1:3]
summary(fit_yr_old)$summary[,'mean'][1:3]
summary(fit_yr_new)$summary[,'mean'][1:3]
summary(fit_yr_dylan)$summary[,'mean'][1:3]

# glm model results
mod <- glm(y ~ x, data = plot_df, family = Gamma(log) )
lines(x_seq, 
      exp( coef(mod)[1] + coef(mod)[2] * x_seq),
      col='red',       lwd=2)

# 
brm_mod <- brm(y~x, data=plot_df, family=Gamma(link="log"),
                prior=c(prior(normal(0,2),class="Intercept"),
                        prior(normal(0,2),class="b"),
                        prior(gamma(0.01,0.01),class="shape")),
                chains=2,iter=1000, cores=4)

brm_fixed <- brm_mod %>% summary %>% .$fixed %>% .[,1]
lines(x_seq, 
      exp( brm_fixed[1] +  brm_fixed[2] * x_seq),
      col='green',       lwd=1)


x <- rgamma(1000, 2.5, 1.5 )
mean(x) 
mean(x)

# test Gammareg ----------------------------------
library(Gammareg)
x1  <- runif(500, 0, 30)
mui <- exp( -5 + 0.2*x1 )
phi <- 0.5
y   <- rgamma(500, shape=phi, scale= mui/phi)


data_sim <- list(
  
  y = y,
  clim_means = x1
  
)