# simulations to 
# retrieve real parameters for a matrix of simulated models based on
# 1. weighted or non-weighted yearly models
# 2. spatially replicated or not
# 3. sd of normal response: 0.3 or 0.5
# 4. beta of effect: 0.45 or 1.2
source("C:/CODE/moving_windows/format_data.R")
library(shinystan)
library(dplyr)
library(tidyr)
library(ggplot2)
library(viridis)
library(testthat)
library(rstan)
library(evd)  
library(purrr)
library(rmutil)
library(parallel)

# set rstan options
rstan_options( auto_write = TRUE )
options( mc.cores = parallel::detectCores() )

# simulations that need be tested
sim_needed <- expand.grid( weight = c(0,    1),
                           spat_n = c(1,    5),
                           sd     = c(0.3,  0.5),
                           beta   = c(0.45, 1.2) )

# calculate realistic y_sd from data ------------------------------------

# get lambda data from Aldo's personal files
lam       <- read.csv("C:/cloud/Dropbox/sApropos/all_demog_6tr.csv", 
                      stringsAsFactors = F) %>% 
                # remove most dodgy datasets
                subset( !(SpeciesAuthor %in% c('Trillium_ovatum',
                                             'Opuntia_macrorhiza_2',
                                             'Cirsium_pitcheri_4',
                                             'Silene_spaldingii') ) )

# SD per study
sd_df <- lam %>% 
            group_by(SpeciesAuthor) %>% 
            summarise( sd_ll  = sd(log_lambda),
                       rep    = n(),
                       rep_yr = MatrixEndYear %>% 
                                unique %>% 
                                length ) %>% 
            ungroup


plot(sd_ll ~ rep, data=sd_df)
plot(sd_ll ~ rep_yr, data=sd_df)

# sd varies from 0.017 to 1.55. Mean 0.36, Median 0.24
sd_df$sd_ll %>% mean
sd_df$sd_ll %>% median
sd_df$sd_ll %>% max
sd_df$sd_ll %>% min

# format raw climate data (messy code) ----------------------------------------

# read 
clim_x    <- read.csv('C:/CODE/climate_drivers_methods/data/demo_airt.csv')

# format: unfortunately this data is replicated daily to accommodate potential daily data
day_one   <- as.Date( paste0("1/1/", first(clim_x$year) ), 
                    format="%d/%m/%Y" ) 

# introduce monthly information
clim_d    <- as.Date(1:nrow(clim_x), day_one-1) %>%
              as.character %>%
              as.data.frame(stringsAsFactors=F) %>%
              separate_(col=".",into=c("year1","month1","day1"),sep="-") %>%
              bind_cols(clim_x) %>%
              dplyr::select(-year,-day) %>%
              setNames( c("year", "month", "day", "species" ,"population", "value") )

# test there is 1 value per month only
clim_d %>% 
  dplyr::select(year,month,value) %>% 
  unique %>% 
  group_by(year,month) %>% 
  summarise(rep=n()) %>% 
  .$rep %>% 
  unique %>% 
  expect_equal( 1 )


# take unique monthly values and calculate anomalies
anom_c <- clim_d %>% 
            dplyr::select(year,month,value) %>% 
            unique %>% 
            # format data in wide form
            spread( month, value ) %>% 
            # calculate anomalies
            select( -year ) %>% 
            apply(2, FUN = scale, center = T, scale = T) %>% 
            # format as data frame
            as.data.frame

# bind three years of data back
yr_1    <- cbind( year = c(1979:2013), anom_c )
yr_2    <- cbind( year = c(1979:2013), anom_c ) %>% 
              setNames( c('year', 13:24) ) %>% 
              mutate( year = year + 1 )
yr_3    <- cbind( year = c(1979:2013), anom_c ) %>% 
              setNames( c('year', 25:36) ) %>% 
              mutate( year = year + 2 )

# Order the anomalies
anom_a  <- Reduce( function(...) inner_join(...), 
                   list(yr_1, yr_2, yr_3) ) %>% 
              .[,c('year', c( paste0(12:10),
                              paste0('0',9:1), 
                              paste0(24:13),
                              paste0(36:25) ) )] %>% 
              select(-year) %>% 
              as.matrix()



# general dataset simulator ----------------------------------


# simulate data changing beta
sim_beta <- function(weight = T,   spat_n = 5,
                     y_sd   = 0.3, beta_x = 1.2){
  
  if( weight ){
  
    # 'true' weight function. Mean in 5th month, sd = 1.
    pdens  <- dnorm(1:12, 5, 1)
    w_v    <- c( (pdens / sum(pdens)) * (1/3),
                 (pdens / sum(pdens)) * (1/3),
                 (pdens / sum(pdens)) * (1/3) )
    
    x1 <- anom_a[,1:12]  %*% w_v[1:12]  %>% as.numeric
    x2 <- anom_a[,13:24] %*% w_v[13:24] %>% as.numeric
    x3 <- anom_a[,25:36] %*% w_v[25:36] %>% as.numeric
    
  }
  
  if( !weight ){
  
    # calculate yearly anomalies 
    x1 <- anom_a[,1:12]  %>% rowMeans
    x2 <- anom_a[,13:24] %>% rowMeans
    x3 <- anom_a[,25:36] %>% rowMeans
    
  }
  
  # function to produce normally distributed response
  prod_y <- function(x){ 
    set.seed(x)
    rnorm(nrow(anom_a), 
          mean= 0 + 
                x1*-beta_x +
                x2*beta_x +
                x3*0, 
           sd = y_sd)
  }
  
  y_vec <- lapply(1:spat_n, prod_y) %>% 
              Reduce(function(...) c(...), .)
    
  # output data
  data.frame( x1 = rep(x1, spat_n),
              x2 = rep(x2, spat_n),
              x3 = rep(x3, spat_n),
              y  = y_vec,
              stringsAsFactors = F )

}

# exploratory plots
par(mfcol = c(2,2), mar   = c(3,3,2,0.1), 
    mgp   = c(1.5,0.5,0) )
plot(y ~ x1, data=sim_beta())
plot(y ~ x2, data=sim_beta())
plot(y ~ x3, data=sim_beta())


# set up model runs 
setup_model_runs <- function( weight,   
                              spat_n,
                              y_sd  , 
                              beta_x ){
  
  # simulate data
  y_sim  <- sim_beta( weight = weight,   
                      spat_n = spat_n,
                      y_sd   = y_sd, 
                      beta_x = beta_x  )$y
  
  # set up climate variables
  clim_mult <- list( anom_a ) %>% 
                 rep( spat_n ) %>% 
                 Reduce( function(...) rbind(...), .)
  
  # organize data into list to pass to stan
  dat_stan <- list(
    n_time  = nrow(clim_mult),
    n_lag   = ncol(clim_mult),
    y       = y_sim,
    clim    = t(clim_mult),
    clim1   = t(clim_mult)[1:12 ,],
    clim2   = t(clim_mult)[13:24,],
    clim3   = t(clim_mult)[25:36,],
    M       = 12,    # number of months in a year
    K       = ncol(clim_mult) / 12,
    expp_beta = 20
  )

  dat_stan
  
}

mod_gaus <- stan_model(file = paste0("code/stan/normal_gaus_nest.stan") )
mod_null <- stan_model(file = paste0("code/stan/normal_null.stan") )

# simulation parameters
sim_pars <- list(
    warmup = 1000, 
    iter   = 4000, 
    thin   = 2, 
    chains = 4
  )


dat_stan <- setup_model_runs( weight = T,   
                              spat_n = 4,
                              y_sd   = 0.3, 
                              beta_x = 1.1 )


fit_null <- sampling(
  object =mod_null,
  data = dat_stan,
  pars = c('alpha', 'y_sd'),
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = 1, #sim_pars$chains,
  init  = list( alpha = 1,
                y_sd = 0.3 )
  #control = list(adapt_delta = 0.999, stepsize = 0.001, max_treedepth = 20)
)

# gaussian moving window
fit_gaus1 <- sampling(
  object =mod_gaus,
  data = dat_stan,
  pars = c('sens_mu', 'sens_sd', 'theta_y', 
           'alpha', 'beta', 'y_sd'),
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = 1, #sim_pars$chains,
  init  = list( beta = 0.5,
                sens_mu = 5,
                sens_sd = 1,
                alpha = 1,
                y_sd = 0.3,
                theta_y = c(1/3),
  # init  = list( beta = c(0.5,0.5,-0.5,-0.5),
  #               sens_mu = rep(5,4),
  #               sens_sd = rep(1,4),
  #               alpha = rep(1,4),
  #               y_sd = rep(0.3,4),
  #               theta_y = rep(c(1/3,1/3,1/3),4) 
  #               ),
  control = list(adapt_delta = 0.999)#,
  #control = list(adapt_delta = 0.999, stepsize = 0.001, max_treedepth = 20)
)
  
fit_gaus2 <- sampling(
  object =mod_gaus,
  data = dat_stan,
  pars = c('sens_mu', 'sens_sd', 'theta_y', 'alpha', 'beta', 'y_sd'),
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains,
  control = list(adapt_delta = 0.999)#,
  #control = list(adapt_delta = 0.999, stepsize = 0.001, max_treedepth = 20)
)

fit_gaus3 <- sampling(
  object =mod_gaus,
  data = dat_stan,
  pars = c('sens_mu', 'sens_sd', 'theta_y', 'alpha', 'beta', 'y_sd'),
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains,
  control = list(adapt_delta = 0.999)#,
  #control = list(adapt_delta = 0.999, stepsize = 0.001, max_treedepth = 20)
)

fit_gaus4 <- sampling(
  object =mod_gaus,
  data = dat_stan,
  pars = c('sens_mu', 'sens_sd', 'theta_y', 'alpha', 'beta', 'y_sd'),
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains,
  control = list(adapt_delta = 0.999)#,
  #control = list(adapt_delta = 0.999, stepsize = 0.001, max_treedepth = 20)
)

fit_gaus5 <- sampling(
  object =mod_gaus,
  data = dat_stan,
  pars = c('sens_mu', 'sens_sd', 'theta_y', 'alpha', 'beta', 'y_sd'),
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains,
  control = list(adapt_delta = 0.999)#,
  #control = list(adapt_delta = 0.999, stepsize = 0.001, max_treedepth = 20)
)

fit_gaus6 <- sampling(
  object =mod_gaus,
  data = dat_stan,
  pars = c('sens_mu', 'sens_sd', 'theta_y', 'alpha', 'beta', 'y_sd'),
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains,
  control = list(adapt_delta = 0.999)#,
  #control = list(adapt_delta = 0.999, stepsize = 0.001, max_treedepth = 20)
)

fit_gaus7 <- sampling(
  object =mod_gaus,
  data = dat_stan,
  pars = c('sens_mu', 'sens_sd', 'theta_y', 'alpha', 'beta', 'y_sd'),
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains,
  control = list(adapt_delta = 0.999)#,
  #control = list(adapt_delta = 0.999, stepsize = 0.001, max_treedepth = 20)
)

fit_gaus8 <- sampling(
  object =mod_gaus,
  data = dat_stan,
  pars = c('sens_mu', 'sens_sd', 'theta_y', 'alpha', 'beta', 'y_sd'),
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains,
  control = list(adapt_delta = 0.999)#,
  #control = list(adapt_delta = 0.999, stepsize = 0.001, max_treedepth = 20)
)

fit_gaus9 <- sampling(
  object =mod_gaus,
  data = dat_stan,
  pars = c('sens_mu', 'sens_sd', 'theta_y', 'alpha', 'beta', 'y_sd'),
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains,
  control = list(adapt_delta = 0.999)#,
  #control = list(adapt_delta = 0.999, stepsize = 0.001, max_treedepth = 20)
)

Sys.time() - init_t 


tiff('results/simulations/bimodal/betas_from_same_data.tiff',
     width=6.3, height=6.3, res=300,unit='in',
     compression='lzw')

par(mfcol = c(3,3), 
    mar   = c(3,3,2,0.1), 
    mgp   = c(1.5,0.5,0) )

extract(fit_gaus1)$beta %>% hist(main = 'beta')
extract(fit_gaus2)$beta %>% hist(main = 'beta')
extract(fit_gaus3)$beta %>% hist(main = 'beta')
extract(fit_gaus4)$beta %>% hist(main = 'beta')
extract(fit_gaus5)$beta %>% hist(main = 'beta')
extract(fit_gaus6)$beta %>% hist(main = 'beta')
extract(fit_gaus7)$beta %>% hist(main = 'beta')
extract(fit_gaus8)$beta %>% hist(main = 'beta')
extract(fit_gaus9)$beta %>% hist(main = 'beta')

dev.off()

# examine chaning
extract(fit_gaus1)$beta %>% plot

shinystan::launch_shinystan(fit_gaus1)

# extract(fit_gaus1)$beta    %>% hist(main = 'beta')
# extract(fit_gaus1)$sens_mu %>% hist(main = 'sens_mu')
# extract(fit_gaus1)$theta_y %>% boxplot(main='yr weights')
# 
# extract(fit_gaus2)$beta    %>% hist(main = 'beta')
# extract(fit_gaus2)$sens_mu %>% hist(main = 'sens_mu')
# extract(fit_gaus2)$theta_y %>% boxplot(main='yr weights')
# 
# extract(fit_gaus3)$beta    %>% hist(main = 'beta')
# extract(fit_gaus3)$sens_mu %>% hist(main = 'sens_mu')
# extract(fit_gaus3)$theta_y %>% boxplot(main='yr weights')
# 
# extract(fit_gaus4)$beta    %>% hist(main = 'beta')
# extract(fit_gaus4)$sens_mu %>% hist(main = 'sens_mu')
# extract(fit_gaus4)$theta_y %>% boxplot(main='yr weights')
# 
# extract(fit_gaus5)$beta    %>% hist(main = 'beta')
# extract(fit_gaus5)$sens_mu %>% hist(main = 'sens_mu')
# extract(fit_gaus5)$theta_y %>% boxplot(main='yr weights')
# 
# dev.off()


library(rstan)

scode <- "
parameters {
  real y[2]; 
} 
model {
  y[1] ~ normal(0, 1);
  y[2] ~ double_exponential(0, 2);
} 
"
