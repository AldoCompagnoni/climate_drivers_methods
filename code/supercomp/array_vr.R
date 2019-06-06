library(dplyr)
library(tidyr)
library(loo)
print('load rstan')
library(rstan)
library(testthat)

# "pipeable" Reduce rbind
rbind_l <- function(x) Reduce(function(...) rbind(...), x)

# arguments from command line
args    <- commandArgs(trailingOnly = TRUE)

print(args)

# climate variable
spp_num  <- args[1]
data_dir <- args[2]
out_dir  <- args[3]
clim_var <- args[4]
response <- args[5]
family   <- args[6]
code_dir <- args[7]
cores    <- args[8]
m_back   <- 36
expp_beta<- 20

# set rstan options
rstan_options( auto_write = TRUE )
options( mc.cores = cores )

# check arguments
args <- args %>% setNames( c('spp_num', 'data_dir','out_dir', 'clim_var',
                             'response','family',  'code_dir','cores'    ) )
print('these are our arguments')
print(args)

print('current working directory')
getwd()

# read format data functions
source(file.path(code_dir,'format_data.R'))

# read data -----------------------------------------------------------------------------------------
lam       <- read.csv( paste0(data_dir, 'all_demog_updt.csv'), stringsAsFactors = F)
m_info    <- read.csv( paste0(data_dir, 'MatrixEndMonth_information.csv'), stringsAsFactors = F)
clim      <- read.csv( paste0(data_dir, clim_var, "_chelsa_prism_hays_2014.csv"),  stringsAsFactors = F)
spp       <- lam$SpeciesAuthor %>% unique

head(lam)
head(m_info)
head(clim)
spp



# format data ---------------------------------------------------------------------------------------
print( "start formatting")

ii            <- as.numeric(args[1])
spp_name      <- spp[ii] # test run w/ spp number 1

print(spp_name)

# lambda data
spp_resp   <- format_species(spp_name, lam, response)

# climate data
clim_separate <- clim_list(spp_name, clim, spp_resp)
clim_detrnded <- lapply(clim_separate, clim_detrend, clim_var)
clim_mats     <- Map(clim_long, clim_detrnded, spp_resp, m_back)

# model data
mod_data          <- lambda_plus_clim(spp_resp, clim_mats, response)
mod_data$climate  <- mod_data$climate #/ diff(range(mod_data$climate))

# throw error if not enough data
if( nrow(mod_data$resp) < 6 ) stop( paste0("not enough temporal replication for '", 
                                              spp_name, "' and response variable '" , response, "'") )

print("done formatting")

# Transform response variables (if needed) ------------------------------------------------------------------

# transform survival/growth - ONLY if less than 30% data points are 1/0
if( grepl("surv", response, ignore.case = T) | 
    grepl("grow", response, ignore.case = T) ){
  
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
if( response == "fec" | response == "RepSSD" ){
  # transform from [0, infinity) to (0, infinity) 
  # I add quantity 2 orders of mag. lower than lowest obs value.
  mod_data$resp[,response] <- mod_data$resp[,response] + 1.54e-12 
} 

if( response == "rho" | response == "react_fsa" ){
  # bound responses to (0,infinity) instead of [1, infinity) 
  mod_data$resp[,response] <- mod_data$resp[,response] - 0.99999
} 


# fit models ----------------------------------------------------------------------------------------
# organize data into list to pass to stan
dat_stan <- list(
  n_time  = nrow(mod_data$climate),
  n_lag   = ncol(mod_data$climate),
  y       = mod_data$resp[,response],
  clim    = mod_data$climate,
  clim_means = rowMeans(mod_data$climate),
  clim_yr = list( rowMeans(mod_data$climate[, 1:12]),
                  rowMeans(mod_data$climate[,13:24]),
                  rowMeans(mod_data$climate[,25:36]) ) %>% do.call(rbind,.),
  M       = 12,    # number of months in a year
  K       = ncol(mod_data$climate) / 12,
  expp_beta = expp_beta
)

# simulation parameters
sim_pars <- list(
  warmup = 1000,
  iter = 4000,
  thin = 2,
  chains = 4
)

print('start first model')
# NULL model (model of the mean)
fit_ctrl1 <- sampling(
  object = readRDS(paste0(family,"_null.RDS")),
  data = dat_stan,
  pars = c('alpha', 'y_sd', 
           'yhat', 'log_lik'),
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains
  #control = list(adapt_delta = 0.999, stepsize = 0.001, max_treedepth = 20)
)
print('end first model')


# year t-1
dat_stan$clim_means  <- rowMeans(mod_data$climate[,13:24])
fit_yr2 <- sampling(
  object = readRDS(paste0(family,"_yr.RDS")),
  data = dat_stan,
  pars = c('alpha', 'beta', 'y_sd',
           'yhat','log_lik'),
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains,
  control = list(adapt_delta = 0.99)
)
print('end yr2')


# year t
dat_stan$clim_means  <- rowMeans(mod_data$climate[,1:12 ])
fit_yr1 <- sampling(
  object = readRDS(paste0(family,"_yr.RDS")),
  data = dat_stan,
  pars = c('alpha', 'beta', 'y_sd',
           'yhat','log_lik'),
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains,
  control = list(adapt_delta = 0.99)
)
print('end yr1')

# year t-2
dat_stan$clim_means  <- rowMeans(mod_data$climate[,25:36])
fit_yr3 <- sampling(
  object = readRDS(paste0(family,"_yr.RDS")),
  data = dat_stan,
  pars = c('alpha', 'beta', 'y_sd',
           'yhat',  'log_lik'),
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains,
  control = list(adapt_delta = 0.99)
)
dat_stan$clim_means  <- rowMeans(mod_data$climate)
print('end yr3')


# gaussian moving window
fit_gaus <- sampling(
  object = readRDS(paste0(family,"_gaus.RDS")),
  data = dat_stan,
  pars = c('sens_mu', 'sens_sd', 'alpha', 'beta', 'y_sd',
           'yhat','log_lik'), #
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains,
  control = list(adapt_delta = 0.99)
)
print('end gaus')


# exponential power moving window
fit_expp <- sampling(
  object = readRDS(paste0(family,"_expp.RDS")),
  data = dat_stan,
  pars = c('sens_mu', 'sens_sd', 'alpha', 'beta', 'y_sd',
           'yhat', 'log_lik'), #'log_lik'
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains,
  control = list(adapt_delta = 0.99)
)
print('end expp')


# moving beta hierarchical
fit_mb_h <- sampling(
  object = readRDS(paste0(family,"_movbeta_h.RDS")),
  data = dat_stan,
  pars = c('alpha', 'beta', 'y_sd', 'mu_beta', 'sigma_beta',
           'yhat', 'log_lik'),
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains,
  control = list(adapt_delta = 0.99)
)
print('end movb_h')


# moving beta correlated
fit_mb <- sampling(
  object = readRDS(paste0(family,"_movbeta.RDS")),
  data = dat_stan,
  pars = c('alpha', 'beta', 'y_sd', 'mu_beta',
           'yhat','log_lik'),
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains,
  control = list(adapt_delta = 0.99)
)
print('end movb')


# moving beta hierarchical
fit_mb_h_n <- sampling(
  object = readRDS(paste0(family,"_movbeta_h_nest.RDS")),
  data = dat_stan,
  pars = c('alpha', 'beta', 'y_sd', 'mu_beta', 'sigma_beta', 'theta_y',
           'yhat',  'log_lik'),
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains,
  control = list(adapt_delta = 0.99)
)
print('end movb_h_n')


# moving beta correlated
fit_mb_n <- sampling(
  object = readRDS(paste0(family,"_movbeta_nest.RDS")),
  data = dat_stan,
  pars = c('alpha', 'beta', 'y_sd', 'mu_beta', 'theta_y',
           'yhat', 'log_lik'),
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains,
  control = list(adapt_delta = 0.99)
)
print('end movb_n')

# Nested models
# update data list
dat_stan$clim         <- t(mod_data$climate)
dat_stan$clim1        <- t(mod_data$climate)[1:12 ,]
dat_stan$clim2        <- t(mod_data$climate)[13:24,]
dat_stan$clim3        <- t(mod_data$climate)[25:36,]

# Simplex nested
fit_24_nest <- sampling(
  object = readRDS(paste0(family,"_dirichlet_nest.RDS")),
  data = dat_stan,
  pars = c('theta_y', 'theta_m', 'alpha', 'beta', 'y_sd',
           'yhat', 'log_lik'),
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains#,
  #control = list(adapt_delta = 0.999, stepsize = 0.001, max_treedepth = 20)
)

# Gaussian nested
fit_gaus_nest <- sampling(
  object = readRDS(paste0(family,"_gaus_nest.RDS")),
  data = dat_stan,
  pars = c('sens_mu', 'sens_sd', 'theta_y', 'alpha', 'beta', 'y_sd',
           'yhat', 'log_lik'),
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains#,
  #control = list(adapt_delta = 0.999, stepsize = 0.001, max_treedepth = 20)
)

# Power exponential nested
fit_expp_nest <- sampling(
  object = readRDS(paste0(family,"_expp_nest.RDS")),
  data = dat_stan,
  pars = c('sens_mu', 'sens_sd', 'theta_y', 'alpha', 'beta', 'y_sd',
           'yhat', 'log_lik'),
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains#,
  #control = list(adapt_delta = 0.999, stepsize = 0.001, max_treedepth = 20)
)


# parameter values and diagnostics ----------------------------------------------------------------

# list of model fits
mod_fit   <- list( ctrl1    = fit_ctrl1,
                   yr1      = fit_yr1,
                   yr2      = fit_yr2,
                   yr3      = fit_yr3,
                   gaus     = fit_gaus,
                   expp     = fit_expp,
                   mb_h     = fit_mb_h,    
                   mb       = fit_mb,
                   mb_h_n   = fit_mb_h_n,  
                   mb_n     = fit_mb_n,
                   simpl_n  = fit_24_nest,
                   expp_n   = fit_expp_nest,
                   gaus_n   = fit_gaus_nest )



# get central tendencies
pars_diag_extract <- function(x){

  # central tendencies
  tmp         <- rstan::extract(x) %>% as.data.frame
  tmp         <- setNames(tmp, gsub("\\.","_",names(tmp)))
  par_means   <- sapply(tmp, function(x) mean(x)) %>%
                    setNames( paste0(names(tmp),"_mean") )
  par_medians <- sapply(tmp, function(x) median(x)) %>%
                    setNames( paste0(names(tmp),"_median") )
  central_tend<- c(par_means, par_medians)

  # diagnostics
  diverg      <- do.call(rbind, args = get_sampler_params(x, inc_warmup = F))[,5]
  n_diverg    <- length(which(diverg == 1))
  df_summ     <- as.data.frame(summary(x)$summary)
  rhat_high   <- length(which(df_summ$Rhat > 1.1))
  n_eff       <- df_summ$n_eff / length(diverg)
  n_eff_low   <- length(which(n_eff < 0.1))
  mcse_high   <- length(which(df_summ$se_mean / df_summ$sd > 0.1))
  diagnostics <- c(n_diverg = n_diverg, rhat_high = rhat_high,
                   n_eff_low = n_eff_low, mcse_high = mcse_high)
  out         <- c( central_tend, diagnostics ) %>% t %>% as.data.frame

  rm(tmp) ; return(out)

}

# extract rhat~n_eff for all parameters
all_diag_extract <- function(x,y){
  
  as.data.frame(summary(x)$summary) %>% 
    select(n_eff, Rhat, se_mean, sd ) %>% 
    tibble::add_column(.before=1, model = y)
  
}

# store posteriors
posterior_extract <- function(model_fit, model_name){

  # central tendencies
  tmp     <- rstan::extract(model_fit) %>% as.data.frame
  post_df <- setNames(tmp, gsub("\\.","_",names(tmp)))
  post_df <- tibble::add_column(post_df,
                                model = model_name, .before=1)

  rm(tmp) ; return(post_df)

}

# calculate central tendencies
pars_diag_l   <- lapply(mod_fit, pars_diag_extract)
mod_pars_diag <- Reduce(function(...) bind_rows(...), pars_diag_l) %>%
                    tibble::add_column(model = names(mod_fit), .before = 1)

# store posteriors
posts_l       <- Map(posterior_extract, mod_fit, names(mod_fit) )
posteriors    <- bind_rows(posts_l)

# extract rhat/neff
diag_df       <- Map( rhat_n_eff_extract, 
                      mod_fit, 
                      names(mod_fit) ) %>% bind_rows

# WAIC model comparison --------------------------------------------------------------------

# wAIC model selection using loo approximation (from library 'loo')
log_liks   <- lapply(mod_fit, extract_log_lik)

# leave-one-out estimates
loo_l      <- lapply(log_liks, loo) %>%
                  setNames( c("loo_ctrl1",
                              "loo_yr1",     "loo_yr2",   "loo_yr3",
                              "loo_gaus",    "loo_expp",
                              "loo_mb_h",    "loo_mb",
                              "loo_mb_h_n",  "loo_mb_n",
                              "loo_simpl_n", "loo_gaus_n", "loo_expp_n") )
loo_df     <- loo::compare(loo_l$loo_ctrl1,
                           loo_l$loo_yr1,     loo_l$loo_yr2,   loo_l$loo_yr3,
                           loo_l$loo_gaus,    loo_l$loo_expp,
                           loo_l$loo_mb_h,    loo_l$loo_mb,
                           loo_l$loo_mb_h_n,  loo_l$loo_mb_n,
                           loo_l$loo_simpl_n, loo_l$loo_gaus_n, loo_l$loo_expp_n ) %>%
                as.data.frame %>%
                tibble::add_column(model = gsub("loo_l\\$loo_","",row.names(.) ), .before = 1) %>% 
                rename( se_diff_loo = se_diff )

# WAIC estimates
waic_l    <- lapply(log_liks, waic) %>%
                setNames(c("waic_ctrl1",
                           "waic_yr1",     "waic_yr2",   "waic_yr3",
                           "waic_gaus",    "waic_expp",
                           "waic_mb_h",    "waic_mb",
                           "waic_mb_h_n",  "waic_mb_n",
                           "waic_simpl_n", "waic_gaus_n", "waic_expp_n") )
waic_df   <- loo::compare(waic_l$waic_ctrl1,
                          waic_l$waic_yr1,     waic_l$waic_yr2,     waic_l$waic_yr3,
                          waic_l$waic_gaus,    waic_l$waic_expp,
                          waic_l$waic_mb_h,    waic_l$waic_mb,
                          waic_l$waic_mb_h_n,  waic_l$waic_mb_n,
                          waic_l$waic_simpl_n, waic_l$waic_gaus_n,   waic_l$waic_expp_n) %>%
                as.data.frame %>%
                tibble::add_column(model = gsub("waic_l\\$waic_","",row.names(.) ),
                                   .before = 1) %>%
                # this is useless to me, causes a conflict down the line
                dplyr::select(-elpd_diff) %>% 
                rename( se_diff_waic = se_diff )


# leave-one-out crossvalidation ------------------------------------------------------------------------

# crossvalidation function
CrossVal <- function(i, mod_data, response){       # i is index for row to leave out

  # identify years
  uniq_yr           <- mod_data$resp$year %>% unique
  test_i            <- which(mod_data$resp$year == uniq_yr[i])

  # put all in matrix form
  x_clim            <- mod_data$climate
  x_clim_means      <- rowMeans(mod_data$climate)   # climate averages over entire window (for control model #2)

  # response variable
  y_train           <- mod_data$resp[-test_i, response]
  y_test            <- mod_data$resp[test_i, response]

  # climate variable
  clim_train        <- x_clim[-test_i,]
  clim_test         <- x_clim[test_i,]

  # climate averages over full 24-month window (for control model #2)
  clim_means_train  <- x_clim_means[-test_i]
  clim_means_test   <- x_clim_means[test_i]

  # organize data into list to pass to stan
  dat_stan_crossval <- list(
    n_train = length(y_train),  # number of data points in train set (length of response var)
    n_test  = length(y_test),   # number of data points in test set
    n_lag   = ncol(clim_train), # maximum lag
    y_train = array(y_train),
    y_test  = array(y_test),
    clim_train = array(clim_train),
    clim_test  = array(clim_test),
    clim_means_train = array(clim_means_train), # climate averages over full 24-month window (for control model #2)
    clim_means_test  = array(clim_means_test),   # climate averages over full 24-month window (for control model #2)
    clim_yr_train    = list( rowMeans(clim_train[, 1:12]),
                             rowMeans(clim_train[,13:24]),
                             rowMeans(clim_train[,25:36]) ) %>% do.call(rbind,.),
    clim_yr_test     = list( rowMeans(clim_test[, 1:12]),
                             rowMeans(clim_test[,13:24]),
                             rowMeans(clim_test[,25:36]) ) %>% do.call(rbind,.),
    expp_beta = expp_beta, # beta paramater for exponential power distribution
    M         = 12,  # number of months in a year
    K         = ncol(clim_train) / 12
  )

  # fit control 1 (intercept only)
  fit_ctrl1_crossval <- sampling(
    object = readRDS(paste0(family,"_null_crossval.RDS")),
    data = dat_stan_crossval,
    pars = c("alpha", "y_sd",
             "pred_y", "log_lik","log_lik_test"),
    warmup = sim_pars$warmup,
    iter = sim_pars$iter,
    thin = sim_pars$thin,
    chains = sim_pars$chains
  )

  # year t
  dat_stan_crossval$clim_means_test   <- rowMeans( mod_data$climate[test_i, 1:12,drop=F] ) %>% array
  dat_stan_crossval$clim_means_train  <- rowMeans( mod_data$climate[-test_i,1:12,drop=F] ) %>% array
  fit_yr1_crossval <- sampling(
    object = readRDS(paste0(family,"_yr_crossval.RDS")),
    data = dat_stan_crossval,
    pars = c("alpha", "beta", "y_sd",
             "pred_y", "log_lik","log_lik_test"),
    warmup = sim_pars$warmup,
    iter = sim_pars$iter,
    thin = sim_pars$thin,
    chains = sim_pars$chains
  )

  # year t-1
  dat_stan_crossval$clim_means_test   <- rowMeans( mod_data$climate[test_i, 13:24,drop=F] ) %>% array
  dat_stan_crossval$clim_means_train  <- rowMeans( mod_data$climate[-test_i,13:24,drop=F] ) %>% array
  fit_yr2_crossval <- sampling(
    object = readRDS(paste0(family,"_yr_crossval.RDS")),
    data = dat_stan_crossval,
    pars = c("alpha", "beta", "y_sd",
             "pred_y", "log_lik","log_lik_test"),
    warmup = sim_pars$warmup,
    iter = sim_pars$iter,
    thin = sim_pars$thin,
    chains = sim_pars$chains
  )

  # year t-2
  dat_stan_crossval$clim_means_test   <- rowMeans( mod_data$climate[test_i, 25:36,drop=F] ) %>% array
  dat_stan_crossval$clim_means_train  <- rowMeans( mod_data$climate[-test_i,25:36,drop=F] ) %>% array
  fit_yr3_crossval <- sampling(
    object = readRDS(paste0(family,"_yr_crossval.RDS")),
    data = dat_stan_crossval,
    pars = c("alpha", "beta", "y_sd",
             "pred_y", "log_lik","log_lik_test"),
    warmup = sim_pars$warmup,
    iter = sim_pars$iter,
    thin = sim_pars$thin,
    chains = sim_pars$chains
  )

  # fit moving window, gaussian
  fit_gaus_crossval <- sampling(
    object = readRDS(paste0(family,"_gaus_crossval.RDS")),
    data = dat_stan_crossval,
    pars = c("sens_mu", "sens_sd", "alpha", "beta", "y_sd",
             "pred_y", "log_lik","log_lik_test"),
    warmup = sim_pars$warmup,
    iter = sim_pars$iter,
    thin = sim_pars$thin,
    chains = sim_pars$chains,
    control = list(adapt_delta = 0.99)
  )

  # fit moving window, exponential power
  fit_expp_crossval <- sampling(
    object = readRDS(paste0(family,"_expp_crossval.RDS")),
    data = dat_stan_crossval,
    pars = c("sens_mu", "sens_sd", "alpha", "beta", "y_sd",
             "pred_y", "log_lik","log_lik_test"),
    warmup = sim_pars$warmup,
    iter = sim_pars$iter,
    thin = sim_pars$thin,
    chains = sim_pars$chains,
    control = list(adapt_delta = 0.99)
  )

  # moving beta model, hierarchical
  fit_mb_h_crossval <- sampling(
    object = readRDS(paste0(family,"_movbeta_h_crossval.RDS")),
    data = dat_stan_crossval,
    pars = c("alpha", "beta", "y_sd", "mu_beta", "sigma_beta",
             "pred_y", "log_lik","log_lik_test"),
    warmup = sim_pars$warmup,
    iter = sim_pars$iter,
    thin = sim_pars$thin,
    chains = sim_pars$chains,
    control = list(adapt_delta = 0.99)
  )

  # moving beta model, NON-hierarchical
  fit_mb_crossval <- sampling(
    object = readRDS(paste0(family,"_movbeta_crossval.RDS")),
    data = dat_stan_crossval,
    pars = c("alpha",  "beta", "y_sd", "mu_beta", "eta", "rho",
             "pred_y", "log_lik","log_lik_test"),
    warmup = sim_pars$warmup,
    iter = sim_pars$iter,
    thin = sim_pars$thin,
    chains = sim_pars$chains,
    control = list(adapt_delta = 0.99)
  )

  # nested moving beta model, hierarchical
  fit_mb_h_n_crossval <- sampling(
    object = readRDS(paste0(family,"_movbeta_h_nest_crossval.RDS")),
    data = dat_stan_crossval,
    pars = c("alpha", "beta", "y_sd", "mu_beta", "sigma_beta",
             "pred_y", "log_lik","log_lik_test"),
    warmup = sim_pars$warmup,
    iter = sim_pars$iter,
    thin = sim_pars$thin,
    chains = sim_pars$chains,
    control = list(adapt_delta = 0.99)
  )

  # nested moving beta model, NON-hierarchical
  fit_mb_n_crossval <- sampling(
    object = readRDS(paste0(family,"_movbeta_nest_crossval.RDS")),
    data = dat_stan_crossval,
    pars = c("alpha",  "beta",   "y_sd", "mu_beta", "eta", "rho",
             "pred_y", "log_lik","log_lik_test"),
    warmup = sim_pars$warmup,
    iter = sim_pars$iter,
    thin = sim_pars$thin,
    chains = sim_pars$chains,
    control = list(adapt_delta = 0.99)
  )

  # update data list
  dat_stan_crossval$clim_train        <- t(dat_stan_crossval$clim_train)
  dat_stan_crossval$clim_test         <- t(dat_stan_crossval$clim_test)

  dat_stan_crossval$clim1_train       <- dat_stan_crossval$clim_train[1:12,]
  dat_stan_crossval$clim1_test        <- dat_stan_crossval$clim_test[ 1:12, ,drop=F]

  dat_stan_crossval$clim2_train       <- dat_stan_crossval$clim_train[13:24,]
  dat_stan_crossval$clim2_test        <- dat_stan_crossval$clim_test[ 13:24, ,drop=F]

  dat_stan_crossval$clim3_train       <- dat_stan_crossval$clim_train[25:36,]
  dat_stan_crossval$clim3_test        <- dat_stan_crossval$clim_test[ 25:36, ,drop=F]

  dat_stan_crossval$clim1_means_train <- colMeans(dat_stan_crossval$clim_train[1:12,])
  dat_stan_crossval$clim1_means_test  <- colMeans(dat_stan_crossval$clim_test[ 1:12, ,drop=F]) %>% array

  dat_stan_crossval$clim2_means_train <- colMeans(dat_stan_crossval$clim_train[13:24,])
  dat_stan_crossval$clim2_means_test  <- colMeans(dat_stan_crossval$clim_test[ 13:24, ,drop=F]) %>% array

  dat_stan_crossval$clim3_means_train <- colMeans(dat_stan_crossval$clim_train[25:36,])
  dat_stan_crossval$clim3_means_test  <- colMeans(dat_stan_crossval$clim_test[ 25:36, ,drop=F]) %>% array

  # fit simplex nested within year
  fit_24_nest_crossval <- sampling(
    object = readRDS(paste0(family,"_dirichlet_nest_crossval.RDS")),
    data = dat_stan_crossval,
    pars = c("theta_m", "theta_y","alpha", "beta", "y_sd", "pred_y", "log_lik","log_lik_test"),
    warmup = sim_pars$warmup,
    iter = sim_pars$iter,
    thin = sim_pars$thin,
    chains = sim_pars$chains,
    control = list(adapt_delta = 0.99)
  )

  
  # fit gaus nested within year
  fit_gaus_nest_crossval <- sampling(
    object = readRDS(paste0(family,"_gaus_nest_crossval.RDS")),
    data = dat_stan_crossval,
    pars = c("sens_mu","sens_sd", "theta_y","alpha", "beta", "y_sd",
             "pred_y", "log_lik","log_lik_test"),
    warmup = sim_pars$warmup,
    iter = sim_pars$iter,
    thin = sim_pars$thin,
    chains = sim_pars$chains,
    control = list(adapt_delta = 0.99)
  )
  
  # fit exponential power nested within year
  fit_expp_nest_crossval <- sampling(
    object = readRDS(paste0(family,"_expp_nest_crossval.RDS")),
    data = dat_stan_crossval,
    pars = c("sens_mu","sens_sd", "theta_y","alpha", "beta", "y_sd", "pred_y", "log_lik","log_lik_test"),
    warmup = sim_pars$warmup,
    iter = sim_pars$iter,
    thin = sim_pars$thin,
    chains = sim_pars$chains,
    control = list(adapt_delta = 0.99)
  )

  # posterior mean prediction for the out-of-sample value
  crossval_mods <- list( ctrl1   = fit_ctrl1_crossval,
                         yr1     = fit_yr1_crossval,
                         yr2     = fit_yr2_crossval,
                         yr3     = fit_yr3_crossval,
                         gaus    = fit_gaus_crossval,
                         expp    = fit_expp_crossval,
                         mb_h    = fit_mb_h_crossval,
                         mb      = fit_mb_crossval,
                         mb_h_n  = fit_mb_h_n_crossval,
                         mb_n    = fit_mb_n_crossval,
                         simpl_n = fit_24_nest_crossval,
                         gaus_n  = fit_gaus_nest_crossval,
                         expp_n  = fit_expp_nest_crossval )

  # predictions
  mod_preds <- lapply(crossval_mods, function(x) rstan::extract(x, "pred_y")$pred_y %>% apply(2,mean) )

  # Expected Log Predictive Density
  mod_elpds <- lapply(crossval_mods, function(x){
                                        rstan::extract(x, "log_lik_test")$log_lik_test %>%
                                          exp %>%
                                          apply(2,mean) %>%
                                          log
                                      } )


  # diagnostics
  diagnostics <- function(fit_obj, name_mod){

    diverg      <- do.call(rbind, args = get_sampler_params(fit_obj, inc_warmup = F))[,5]
    n_diverg    <- length(which(diverg == 1))
    df_summ     <- as.data.frame(summary(fit_obj)$summary)
    rhat_high   <- length(which(df_summ$Rhat > 1.1))
    n_eff       <- df_summ$n_eff / length(diverg)
    n_eff_low   <- length(which(n_eff < 0.1))
    mcse_high   <- length(which(df_summ$se_mean / df_summ$sd > 0.1))
    out         <- data.frame(n_diverg,  rhat_high,
                              n_eff_low, mcse_high)
    # out         <- setNames(out, paste0(names(out),"_",name_mod) )
    return(out)

  }

  # store diagnostics
  # gaus_expp   <- crossval_mods[c("gaus","expp","gev","simpl")]
  diagnost_l  <- Map(diagnostics, crossval_mods, names(crossval_mods))
  diagnost_df <- do.call(cbind, diagnost_l) %>%
                    bind_cols( unique( dplyr::select(mod_data$resp[test_i,],year) ) ) %>%
                    setNames( gsub("\\.","_",names(.)) )

  # function
  pred_elpd_df<- mod_data$resp[test_i,] %>%
                    mutate( # predictions
                            ctrl1_pred   = mod_preds$ctrl1,
                            yr1_pred     = mod_preds$yr1,
                            yr2_pred     = mod_preds$yr2,
                            yr3_pred     = mod_preds$yr3,
                            gaus_pred    = mod_preds$gaus,
                            expp_pred    = mod_preds$expp,
                            mb_pred      = mod_preds$mb,
                            mb_h_pred    = mod_preds$mb_h,
                            mb_n_pred    = mod_preds$mb_n,
                            mb_h_n_pred  = mod_preds$mb_h_n,
                            simpl_n_pred = mod_preds$simpl_n,
                            gaus_n_pred  = mod_preds$gaus_n,
                            expp_n_pred  = mod_preds$expp_n,

                            # Expected Log Predictive Density
                            ctrl1_elpd   = mod_elpds$ctrl1,
                            yr1_elpd     = mod_elpds$yr1,
                            yr2_elpd     = mod_elpds$yr2,
                            yr3_elpd     = mod_elpds$yr3,
                            gaus_elpd    = mod_elpds$gaus,
                            expp_elpd    = mod_elpds$expp,
                            mb_elpd      = mod_elpds$mb,
                            mb_h_elpd    = mod_elpds$mb_h,
                            mb_n_elpd    = mod_elpds$mb_n,
                            mb_h_n_elpd  = mod_elpds$mb_h_n,
                            simpl_n_elpd = mod_elpds$simpl_n,
                            gaus_n_elpd  = mod_elpds$gaus_n,
                            expp_n_elpd  = mod_elpds$expp_n )

  # df to return
  out         <- left_join(pred_elpd_df, diagnost_df)

  # remove stanfit objects (garbage collection)
  rm(fit_ctrl1_crossval)
  rm(fit_yr1_crossval)
  rm(fit_yr2_crossval)
  rm(fit_yr3_crossval)
  rm(fit_gaus_crossval)
  rm(fit_expp_crossval)
  rm(fit_mb_h_crossval)
  rm(fit_mb_crossval)
  rm(fit_mb_h_n_crossval)
  rm(fit_mb_n_crossval)
  rm(fit_24_nest_crossval)
  rm(fit_gaus_nest_crossval)
  rm(fit_expp_nest_crossval)

  return(out)

}

print( "start crossvalidation" )

# spp-specific cross validation
year_inds   <- seq_along(unique(mod_data$resp$year))
cxval_res   <- lapply( year_inds, CrossVal, mod_data, response)
cxval_pred  <- do.call(rbind, cxval_res)

print( "end crossvalidation" )

# measures of fit --------------------------------------------------------------------------

# calculate either mse or deviance
pred_perform <- function(x, mod_data, response, type){

  if( type == "mse"){
    res   <- (x - mod_data$resp[,response])^2 %>% mean
  }
  if(type == "deviance"){
    res   <-calc.deviance(x, mod_data$resp[,response],
                          weights = rep(1, length(x) ),
                          family="gaussian", calc.mean = TRUE)
  }

  return(res)

}

# format results into a data frame
perform_format <- function(x, var){

  x %>%
    unlist %>%
    t %>%
    t %>%
    as.data.frame %>%
    tibble::rownames_to_column(var = "model") %>%
    mutate( model = gsub("mod_preds.", "", model))  %>%
    setNames( c("model", var) )

}

# order of 'mod_data$resp' and 'cxval_pred' need be the same
expect_true( all.equal(dplyr::select(mod_data$resp, year, population),
                       dplyr::select(cxval_pred,    year, population)) )

print('IMPORTANT: Check names')
names(cxval_pred)

# mean squared error
mse <- cxval_pred %>%
          dplyr::select(ctrl1_pred:expp_n_pred) %>%
          lapply(pred_perform, mod_data, response, "mse") %>%
          perform_format("mse") %>%
          mutate( model = gsub("_pred","",model) )

# Expected Log Predictive Density
elpd <- cxval_pred %>%
          dplyr::select(ctrl1_pred:expp_n_elpd) %>%
          apply(2, sum) %>%
          as.matrix %>%
          as.data.frame %>%
          tibble::add_column(model = rownames(.), .before=1) %>%
          mutate( model = gsub("_elpd","",model) ) %>%
          setNames( c("model", "elpd") )

# measures of fit
# mof  <- merge(mse, devi)
mof  <- merge(mse, elpd)

# store results ---------------------------------------------------------------------------
mod_summs <- Reduce(function(...) full_join(...),
                    list(mod_pars_diag, loo_df, waic_df, mof) ) %>%
                    arrange( mse )

write.csv(mod_summs,  paste0(out_dir,"-",spp_num,"_mod_summaries_",     spp_name,"_",response,"_",clim_var,".csv"), row.names = F)
write.csv(posteriors, paste0(out_dir,"-",spp_num,"_posterior_",         spp_name,"_",response,"_",clim_var,".csv"), row.names = F)
write.csv(cxval_pred, paste0(out_dir,"-",spp_num,"_crossval_pred_diag_",spp_name,"_",response,"_",clim_var,".csv"), row.names = F)
write.csv(diag_df,    paste0(out_dir,"-",spp_num,"_diagnostics_",spp_name,"_",response,"_",clim_var,".csv"), row.names = F)
save.image( paste0(out_file,'_workshace_',spp_num,'.RData') )
