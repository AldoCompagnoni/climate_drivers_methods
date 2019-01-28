source("C:/CODE/moving_windows/format_data.R")
library(tidyverse)
library(testthat)
library(rstan)
library(evd)
library(rmutil)

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
month_rep <- clim_d %>% 
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
yr_1    <- cbind(year=c(1979:2013), anom_c )
yr_2    <- cbind(year=c(1979:2013), anom_c ) %>% 
              setNames( c('year', 13:24) ) %>% 
              mutate( year = year + 1 )
yr_3    <- cbind(year=c(1979:2013), anom_c ) %>% 
              setNames( c('year', 25:36) ) %>% 
              mutate( year = year + 2 )

# Order the anomalies
anom_a  <- Reduce( function(...) inner_join(...), 
                   list(yr_1, yr_2, yr_3) ) %>% 
              .[,c('year', c( paste0(12:10),
                              paste0('0',9:1), 
                              paste0(24:13),
                              paste0(36:25) ) )]


# simulate data ------------------------------------------------

# "true" climate effects
betas <- seq(0.2,2.2,by=0.5)

# list to store models 
mod_l <- list( '0.2' = NULL,
               '0.7' = NULL,
               '1.2' = NULL,
               '1.7' = NULL,
               '2.2' = NULL )

# 'true' weight function. Mean in 5th month, sd = 1.
pdens  <- dnorm(1:36, 5, 1)
w_v    <- (pdens / sum(pdens))
# climate "data"
clim_m <- select(anom_a, -year) %>% as.matrix

# simulate data
sim_data <- function(beta_x){
  
  x      <- clim_m %*% w_v

  # function to produce normally distributed response
  prod_y <- function(x) rnorm(length(x), mean=0+x*beta_x,sd=1)
  
  # output data
  data.frame( x = x,
              y = prod_y(x),
              stringsAsFactors = F )

}

# plot it
plot(y ~ x, data=sim_data(0.7))


# fit data ---------------------------------------------------

# set rstan options
rstan_options( auto_write = TRUE )
options( mc.cores = parallel::detectCores() )

# simulation parameters
sim_pars <- list(
  warmup = 1000, 
  iter   = 4000, 
  thin   = 2, 
  chains = 3
)

# here we focus on a normally distributed response
family <- 'normal'

# fit models based on data
fit_mods <- function(df_x){
  
  # organize data into list to pass to stan
  dat_stan <- list(
    n_time  = nrow(clim_m),
    n_lag   = ncol(clim_m),
    y       = df_x$y,
    clim    = clim_m,
    M       = 12,    # number of months in a year
    K       = ncol(clim_m) / 12,
    expp_beta = 20
  )
  
  # gaussian moving window
  fit_gaus <- stan(
    file = paste0("code/stan/",family,"_gaus.stan"),
    data = dat_stan,
    pars = c('sens_mu', 'sens_sd', 'alpha', 'beta', 'y_sd'),
    warmup = sim_pars$warmup,
    iter = sim_pars$iter,
    thin = sim_pars$thin,
    chains = sim_pars$chains,
    control = list(adapt_delta = 0.999)#,
    #control = list(adapt_delta = 0.999, stepsize = 0.001, max_treedepth = 20)
  )
  
  # exponential power moving window
  fit_expp <- stan(
    file = paste0("code/stan/",family,"_expp.stan"),
    data = dat_stan,
    pars = c('sens_mu', 'sens_sd', 'alpha', 'beta', 'y_sd'),
    warmup = sim_pars$warmup,
    iter = sim_pars$iter,
    thin = sim_pars$thin,
    chains = sim_pars$chains,
    control = list(adapt_delta = 0.999)#,
    #control = list(adapt_delta = 0.999, stepsize = 0.001, max_treedepth = 20)
  )
  
  # Generalized extreme value 
  fit_gev <- stan(
    file = paste0("code/stan/",family,"_gev.stan"),
    data = dat_stan,
    pars = c('loc', 'scale', "shape", 'alpha', 'beta', 'y_sd'),
    warmup = sim_pars$warmup,
    iter = sim_pars$iter,
    thin = sim_pars$thin,
    chains = sim_pars$chains,
    control = list(adapt_delta = 0.999)#,
    #control = list(adapt_delta = 0.999, stepsize = 0.001, max_treedepth = 20)
  )
  
  # Nested models: update data list
  dat_stan$clim         <- t(clim_m)
  dat_stan$clim1        <- t(clim_m)[1:12 ,]
  dat_stan$clim2        <- t(clim_m)[13:24,]
  dat_stan$clim3        <- t(clim_m)[25:36,]
  
  # Simplex nested
  fit_sad_nest <- stan(
    file = paste0("code/stan/",family,"_dirichlet_nest.stan"),
    data = dat_stan,
    pars = c('theta_y', 'theta_m', 'alpha', 'beta', 'y_sd'),
    warmup = sim_pars$warmup,
    iter = sim_pars$iter,
    thin = sim_pars$thin,
    chains = sim_pars$chains,
    control = list(adapt_delta = 0.999)
  )

  list( fit_gaus, 
        fit_expp,
        fit_gev,
        fit_sad_nest )
  
}

mod_l$'0.2' <- fit_mods( sim_data(0.2) )
mod_l$'0.7' <- fit_mods( sim_data(0.7) )
mod_l$'1.2' <- fit_mods( sim_data(1.2) )
mod_l$'1.7' <- fit_mods( sim_data(1.7) )
mod_l$'2.2' <- fit_mods( sim_data(2.2) )


# model evaluation ----------------------------------------

# get divergent transitions
get_div  <- function(x){
  get_sampler_params(x, inc_warmup=T)[[1]][,'divergent__'] %>% 
    sum
}
# number of params with Rhat issues 
get_rhat_n <- function(x){
  x %>% 
    summary %>% 
    .$summary %>% 
    .[,'Rhat'] %>% 
    `>` (1.01) %>% 
    # do parameters have Rhat>1.01?
    sum %>% 
    # how many parameters have Rhat>1.01?
    sum
}
# get a parameter
get_par <- function(x,param){
  x %>% 
   summary %>% 
  .$summary %>% 
  .[param,'mean']
}

# get

# get divergent transitions
sapply(mod_l$'0.2',get_div)
sapply(mod_l$'0.7',get_div)
sapply(mod_l$'1.2',get_div)
sapply(mod_l$'1.7',get_div)
sapply(mod_l$'2.2',get_div)

# get Rhat
sapply(mod_l$'0.2',get_rhat_n)
sapply(mod_l$'0.7',get_rhat_n)
sapply(mod_l$'1.2',get_rhat_n)
sapply(mod_l$'1.7',get_rhat_n)
sapply(mod_l$'2.2',get_rhat_n)

# estimated betas
beta_est <- c( sapply(mod_l$'0.2',get_par, 'beta'),
               sapply(mod_l$'0.7',get_par, 'beta'),
               sapply(mod_l$'1.2',get_par, 'beta'),
               sapply(mod_l$'1.7',get_par, 'beta'),
               sapply(mod_l$'2.2',get_par, 'beta') )

# data frame of betas
b_df <- data.frame( beta_true = c(rep(0.2,4),rep(0.7,4),
                                  rep(1.2,4),rep(1.7,4),rep(2.2,4)),
                    mod       = rep(c('gaus','expp','gev','sad'),5),
                    beta_est  = beta_est )

# store figure
tiff('results/simulations/beta_true_est_n_lag_sd.tiff',
     unit="in", width=6.3, height=6.3, res=600,compression="lzw" )
par(mar=c(3.5,3.5,0.2,0.2), 
    mgp=c(1.7,0.7,0))
plot(beta_est ~ beta_true, 
     data=b_df, pch = 16,
     xlab = 'True beta value', ylab = 'Estimated beta value',
     col = as.numeric(mod) )
legend('topleft',
       c('gaus','expp','gev','sad'),
       pch = 16, lwd=3,
       col = b_df$mod %>% unique %>% as.numeric,
       bty='n')
abline(0,1,lwd=3)
dev.off()



# get weight params 
get_weight <- function(x){
  gaus <- c(x[[1]] %>% get_par('sens_mu'),
            x[[1]] %>% get_par('sens_sd') )
  expp <- c(x[[2]] %>% get_par('sens_mu'),
            x[[2]] %>% get_par('sens_sd') )
  gev  <- c(x[[3]] %>% get_par('loc'),
            x[[3]] %>% get_par('scale'),
            x[[3]] %>% get_par('shape') )
  sad  <- c(x[[4]] %>% get_par('theta_m[1]'),
            x[[4]] %>% get_par('theta_m[2]'),
            x[[4]] %>% get_par('theta_m[3]'),
            x[[4]] %>% get_par('theta_m[4]'),
            x[[4]] %>% get_par('theta_m[5]'),
            x[[4]] %>% get_par('theta_m[6]'),
            x[[4]] %>% get_par('theta_m[7]'),
            x[[4]] %>% get_par('theta_m[8]'),
            x[[4]] %>% get_par('theta_m[9]'),
            x[[4]] %>% get_par('theta_m[10]'),
            x[[4]] %>% get_par('theta_m[11]'),
            x[[4]] %>% get_par('theta_m[12]'),
            x[[4]] %>% get_par('theta_y[1]'),
            x[[4]] %>% get_par('theta_y[2]'),
            x[[4]] %>% get_par('theta_y[3]') )
  
  list(gaus=gaus,expp=expp,gev=gev,sad=sad)
  
}

# exponential power distribution
dexppow <- function(x, mu, sigma, beta) {
    return((beta / (2 * sigma * gamma (1.0/beta)) ) *
             exp(-(abs(x - mu)/sigma)^beta));
}

# plot weights
plot_w <- function(weight_l){
  xx <- 1:36
  # gaus
  gaus_w_v <- dnorm(xx, weight_l$gaus[1], weight_l$gaus[2])
  gaus_w_v <- gaus_w_v / sum(gaus_w_v)
  # expp
  expp_w_v <- dexppow(xx, weight_l$gaus[1], weight_l$gaus[2], 20)
  expp_w_v <- expp_w_v / sum(expp_w_v)
  # gev
  gev_w_v  <- dgev(xx,weight_l$gev[1], 
                   weight_l$gev[2], 
                   weight_l$gev[3])
  gev_w_v  <- gev_w_v / sum(gev_w_v)
  #sad
  sad_w_v  <- c(weight_l$sad[1:12]*weight_l$sad[13],
                weight_l$sad[1:12]*weight_l$sad[14],
                weight_l$sad[1:12]*weight_l$sad[15])
  
  plot(1:36, gaus_w_v,type='l',col='red',lwd=2,
       ylim=c(0,max(w_v)+0.01),ylab='Weight',xlab="Month")
  
  lines(1:36,expp_w_v,lwd=2,col='black')
  lines(1:36,gev_w_v,lwd=2,col='green')
  
  lines(1:36,sad_w_v,lwd=2,col='blue')
  lines(1:36,w_v,lty=2,lwd=2,col='grey')
  
}

tiff('results/simulations/month_w_n_lag_sd.tiff',
     unit="in", width=6.3, height=8, res=600,compression="lzw" )
par(mfrow=c(3,2),
    mar=c(3,3,0.2,0.2), 
    mgp=c(1.3,0.3,0) )
plot_w(get_weight(mod_l$'0.2'))
text(25,0.4,pos=1,'beta=0.2')
plot_w(get_weight(mod_l$'0.7'))
text(25,0.4,pos=1,'beta=0.7')
plot_w(get_weight(mod_l$'1.2'))
text(15,0.4,pos=1,'beta=1.2')
legend('topright',
       c('gaus','expp','gev','sad','true'),
       col=c('red','black','green','blue','grey'),
       bty='n', lty=1, cex=0.8)
plot_w(get_weight(mod_l$'1.7'))
text(25,0.4,pos=1,'beta=1.7')
plot_w(get_weight(mod_l$'2.2'))
text(25,0.4,pos=1,'beta=2.2')
dev.off()


vecc <- dgev(1:36,5,1,0.4)
lines(1:36,(vecc / sum(vecc)),col='magenta')

# par_sad <- fit_sad_nest %>% 
#                 rstan::extract() %>%
#                 # rstan::extract(permuted=FALSE) %>%
#                 as.data.frame %>% 
#                 # point out iterations
#                 mutate( iter = 1:nrow(.) )
# 
# par_sad %>% 
#   select( grep('theta_m',names(.)) ) %>% 
#   boxplot( ylim = c(0,1) )
# 
# par_sad %>% 
#   select( grep('theta_y',names(.)) ) %>% 
#   boxplot
