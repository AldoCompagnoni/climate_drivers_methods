source("C:/CODE/moving_windows/format_data.R")
library(tidyverse)
library(testthat)
library(rstan)
library(evd)
library(rmutil)

# format raw climate data (messy code) ----------------------------------------

# set SD of simulated log_lambda
log_lambda_sd <- 0.3

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
anom_a  <- yr_1
              


# simulate data ------------------------------------------------

# "true" climate effects
betas <- seq(0.2,2.2,by=0.5)

# list to store models 
mod_l <- list( '0.2'  = NULL,
               '0.45' = NULL,
               '0.7'  = NULL,
               '1.2'  = NULL,
               '1.7'  = NULL,
               '2.2'  = NULL )

# 'true' weight function. Mean in 5th month, sd = 1.
pdens  <- dnorm(1:12, 5, 1)
w_v    <- (pdens / sum(pdens))
# chose year 1 as climate "data"
clim_m <- dplyr::select(anom_a,-year) %>% 
            as.matrix %>% 
            # remove 2 years: 33 total, 
            # like simulations 36 months back
            .[-c(1:2),]


# simulate data
sim_data <- function(beta_x){
  
  x      <- clim_m %*% w_v

  # function to produce normally distributed response
  prod_y <- function(x) rnorm(length(x), 
                              mean=0+x*beta_x,
                              sd=log_lambda_sd)
  
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

# compile the models
mod_gaus <- stan_model(file = paste0("code/stan/",family,"_gaus.stan") )
mod_expp <- stan_model(file = paste0("code/stan/",family,"_expp.stan") )
mod_gev  <- stan_model(file = paste0("code/stan/",family,"_gev.stan")  )
mod_sad  <- stan_model(file = paste0("code/stan/",family,"_dirichlet.stan") )


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
  fit_gaus <- sampling(
    object =mod_gaus,
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
  fit_expp <- sampling(
    object = mod_expp,
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
  fit_gev <- sampling(
    object = mod_gev,
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
  
  # Simplex nested
  fit_sad <- sampling(
    object = mod_sad,
    data = dat_stan,
    pars = c('theta', 'alpha', 'beta', 'y_sd'),
    warmup = sim_pars$warmup,
    iter = sim_pars$iter,
    thin = sim_pars$thin,
    chains = sim_pars$chains,
    control = list(adapt_delta = 0.999)
  )

  list( fit_gaus, 
        fit_expp,
        fit_gev,
        fit_sad )
  
}

mod_l$'0.2'  <- fit_mods( sim_data(0.2) )
mod_l$'0.45' <- fit_mods( sim_data(0.45) )
mod_l$'0.7'  <- fit_mods( sim_data(0.7) )
mod_l$'1.2'  <- fit_mods( sim_data(1.2) )
mod_l$'1.7'  <- fit_mods( sim_data(1.7) )
mod_l$'2.2'  <- fit_mods( sim_data(2.2) )

# save this image
save.image( paste0('results/simulations/',
                   'yr1_sd_',
                   log_lambda_sd,'.Rdata') )


# model evaluation ----------------------------------------

# get divergent transitions
get_div  <- function(x){
  get_sampler_params(x,inc_warmup=F) %>% 
    # get divergent transitions for each chain
    sapply( function(x) x[,'divergent__'] ) %>% 
    as.vector %>% 
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
# format the list of diagnostics
format_list <- function(diag_l){
  do.call(rbind, diag_l) %>% 
    as.data.frame %>% 
    setNames( c('gaus', 'expp', 'gev', 'sad') ) %>% 
    mutate( beta = c(0.2,0.45,0.7,1.2,1.7,2.2) ) %>% 
    .[,c('beta', 'gaus', 'expp', 'gev', 'sad')]
}
  
# get divergent transitions
div_df <- list( sapply(mod_l$'0.2',get_div), 
                sapply(mod_l$'0.45',get_div), 
                sapply(mod_l$'0.7',get_div), 
                sapply(mod_l$'1.2',get_div), 
                sapply(mod_l$'1.7',get_div), 
                sapply(mod_l$'2.2',get_div) ) %>% 
            format_list 

# get Rhat
rhat_df<- list(sapply(mod_l$'0.2',get_rhat_n),
               sapply(mod_l$'0.45',get_rhat_n), 
               sapply(mod_l$'0.7',get_rhat_n),
               sapply(mod_l$'1.2',get_rhat_n),
               sapply(mod_l$'1.7',get_rhat_n),
               sapply(mod_l$'2.2',get_rhat_n) ) %>% 
            format_list

# store tables 
write.csv(div_df, 
          paste0('results/simulations/divergent_1yr_sd_',
                 log_lambda_sd,'.csv'), row.names=F)
write.csv(rhat_df, 
          paste0('results/simulations/rhat_1yr_sd_',
                 log_lambda_sd,'.csv'), row.names=F)

# betas ---------------------------------------------------------

# estimated betas
beta_est <- c( sapply(mod_l$'0.2', get_par, 'beta'),
               sapply(mod_l$'0.45',get_par, 'beta'),
               sapply(mod_l$'0.7', get_par, 'beta'),
               sapply(mod_l$'1.2', get_par, 'beta'),
               sapply(mod_l$'1.7', get_par, 'beta'),
               sapply(mod_l$'2.2', get_par, 'beta') )


# data frame of betas
b_df <- data.frame( beta_true = c(rep(0.2,4),
                                  rep(0.45,4),
                                  rep(0.7,4),
                                  rep(1.2,4),
                                  rep(1.7,4),rep(2.2,4)),
                    mod       = rep(c('gaus','expp','gev','sad'),6),
                    beta_est  = beta_est ) %>% 
          mutate( mod = as.factor(mod))

# store figure
tiff(paste0('results/simulations/betas/beta_true_est_yr1_sd_',
            log_lambda_sd,'.tiff'),
     unit="in", width=6.3, height=6.3, res=600,compression="lzw" )
par(mfrow=c(1,1), mar=c(3.5,3.5,0.2,0.2), 
    mgp=c(1.7,0.7,0))
plot(beta_est ~ jitter(beta_true), 
     data=b_df, pch = 16,
     xlab = 'True beta value', 
     ylab = 'Estimated beta value',
     col = as.numeric(as.factor(mod)) )
legend('topleft',
       c('gaus','expp','gev','sad'),
       pch = 16, lwd=3,
       col = b_df$mod %>% unique %>% as.numeric,
       bty='n')
abline(0,1,lwd=3)
dev.off()


# stack post
stack_post <- function(x, beta_val, param){
  
  # get the posterior
  get_post <- function(x,param){
    x %>% 
      rstan::extract() %>% 
      as.data.frame %>% 
      .[,param]  
  }

  # get posterior and "stack it "gather" it
  x %>% 
    purrr::pluck( as.character(beta_val) ) %>% 
    sapply(get_post, param) %>% 
    as.data.frame %>% 
    setNames( c('gaus','expp','gev','sad') ) %>% 
    stack %>% 
    rename( model = ind,
            beta  = values ) %>% 
    mutate( true_beta = beta_val )
  
}

# extract posteriors
beta_post <- list( stack_post(mod_l, 0.2, 'beta'),
                   stack_post(mod_l, 0.45, 'beta'),
                   stack_post(mod_l, 0.7, 'beta'),
                   stack_post(mod_l, 1.2, 'beta'),
                   stack_post(mod_l, 1.7, 'beta'),
                   stack_post(mod_l, 2.2, 'beta') ) %>% 
                bind_rows %>% 
                group_by( true_beta, model ) %>% 
                summarise( beta_m  = mean(beta),
                           beta_05 = quantile(beta, probs = 0.05),
                           beta_95 = quantile(beta, probs = 0.95) ) %>% 
                ungroup %>% 
                # offset true betas for representation
                mutate( true_beta = replace(true_beta,
                                            model == 'gaus',
                                            true_beta[model == 'gaus']-0.04)) %>% 
                mutate( true_beta = replace(true_beta,
                                            model == 'expp',
                                            true_beta[model == 'expp']-0.02)) %>% 
                mutate( true_beta = replace(true_beta,
                                            model == 'gev',
                                            true_beta[model == 'gev']+0.02)) %>% 
                mutate( true_beta = replace(true_beta,
                                            model == 'sad',
                                            true_beta[model == 'sad']+0.04))

# Create and save 
ggplot( mutate(beta_post, 
               true_beta = true_beta ), 
  aes(true_beta, beta_m) ) +
  viridis::scale_color_viridis(discrete=TRUE) + 
  geom_pointrange( aes(ymin = beta_05,
                       ymax = beta_95,
                       color = model ) ) + 
  geom_abline() + 
  ylab('Estimated beta') +
  xlab('True beta') + 
  ggsave(paste0('results/simulations/betas/',
                'beta_true_est_yr1_post_sd_',
                 log_lambda_sd,'.tiff'),
         width = 6.3, height = 6.3)


# weights ------------------------------------------------------

# get weight params 
get_weight <- function(x){
  gaus <- c(x[[1]] %>% get_par('sens_mu'),
            x[[1]] %>% get_par('sens_sd') )
  expp <- c(x[[2]] %>% get_par('sens_mu'),
            x[[2]] %>% get_par('sens_sd') )
  gev  <- c(x[[3]] %>% get_par('loc'),
            x[[3]] %>% get_par('scale'),
            x[[3]] %>% get_par('shape') )
  sad  <- c(x[[4]] %>% get_par('theta[1]'),
            x[[4]] %>% get_par('theta[2]'),
            x[[4]] %>% get_par('theta[3]'),
            x[[4]] %>% get_par('theta[4]'),
            x[[4]] %>% get_par('theta[5]'),
            x[[4]] %>% get_par('theta[6]'),
            x[[4]] %>% get_par('theta[7]'),
            x[[4]] %>% get_par('theta[8]'),
            x[[4]] %>% get_par('theta[9]'),
            x[[4]] %>% get_par('theta[10]'),
            x[[4]] %>% get_par('theta[11]'),
            x[[4]] %>% get_par('theta[12]'))
  
  list(gaus=gaus,expp=expp,gev=gev,sad=sad)
  
}

# exponential power distribution
dexppow <- function(x, mu, sigma, beta) {
    return((beta / (2 * sigma * gamma (1.0/beta)) ) *
             exp(-(abs(x - mu)/sigma)^beta));
}

# plot weights
plot_w <- function(weight_l){
  
  xx       <- 1:12
  # gaus
  gaus_w_v <- dnorm(xx, weight_l$gaus[1], weight_l$gaus[2])
  gaus_w_v <- gaus_w_v / sum(gaus_w_v)
  # expp
  expp_w_v <- dexppow(xx, weight_l$expp[1], weight_l$expp[2], 20)
  expp_w_v <- expp_w_v / sum(expp_w_v)
  # gev
  gev_w_v  <- dgev(xx,weight_l$gev[1], 
                      weight_l$gev[2], 
                      weight_l$gev[3])
  gev_w_v  <- gev_w_v / sum(gev_w_v)
  #sad
  sad_w_v  <- c( weight_l$sad[1:12] )
  
  plot(xx, gaus_w_v,type='l',col='red',lwd=2,
       ylim=c(0,max(w_v)+0.01),ylab='Weight',xlab="Month")
  
  lines(xx,expp_w_v,lwd=2,col='magenta')
  lines(xx,gev_w_v,lwd=2,col='green')
  
  lines(xx,sad_w_v,lwd=2,col='blue')
  lines(xx,w_v,lty=2,lwd=3,col='black')
  
}


tiff(paste0('results/simulations/weights/month_w_yr1_sd_',
            log_lambda_sd,'.tiff'),
     unit="in", width=6.3, height=8, res=600,compression="lzw" )
par(mfrow=c(3,2),
    mar=c(3,3,0.2,0.2), 
    mgp=c(1.3,0.3,0) )
plot_w(get_weight(mod_l$'0.2')) 
co <- par('usr')
text(10, co[4],pos=1,'beta=0.2', cex = 2)
plot_w(get_weight(mod_l$'0.45'))
text(10,co[4],pos=1,'beta=0.45', cex = 2)
plot_w(get_weight(mod_l$'0.7'))
text(10, co[4],pos=1,'beta=0.7', cex = 2)
plot_w(get_weight(mod_l$'1.2'))
text(2.5,co[4],pos=1,'beta=1.2', cex = 2)
legend('topright',
       c('gaus','expp','gev','sad','true'),
       col=c('red','magenta','green','blue','black'),
       bty='n', lwd=c(2,2,2,2,3),
                lty=c(1,1,1,1,2), cex=2)
plot_w(get_weight(mod_l$'1.7'))
text(10,co[4],pos=1,'beta=1.7', cex = 2)
plot_w(get_weight(mod_l$'2.2'))
text(10,co[4],pos=1,'beta=2.2', cex = 2)
dev.off()


# examine chains ------------------------------------

# looks at divergent transitions and problem Rhats
div_df
rhat_df

mod_bet <- c('0.2','0.45','0.7','1.2','1.7','2.2')

# choose model (gaus, expp, ...)
for(mod_ii in 1:4){
  
  # chose beta
  for(ii in 1:length(mod_bet)){
    
    par_ch <- mod_l[mod_bet[ii]][[1]][[mod_ii]] %>% 
                  rstan::extract(permuted=F) %>% 
                  as.data.frame %>% 
                  mutate( iter = 1:nrow(.) )
    
    if( mod_ii == 1 ){
      tiff(paste0('results/simulations/chains/',
                  'gaus/beta_',
                  mod_bet[ii],'_yr1_sd_',log_lambda_sd,'.tiff'),
           unit="in", width=6.3, height=8, res=600,compression="lzw" )
      par(mfcol=c(3,2))
      plot(par_ch$iter,par_ch$'chain:1.sens_mu',
           type='l',ylab="sens_mu",ylim =c(0,12),
           xlab="Iteraction")
      plot(par_ch$iter,par_ch$'chain:2.sens_mu',
           type='l',ylab="sens_mu",ylim =c(0,12),
           xlab="Iteraction")
      plot(par_ch$iter,par_ch$'chain:3.sens_mu',
           type='l',ylab="sens_mu",ylim =c(0,12),
           xlab="Iteraction")
      
      plot(par_ch$iter,par_ch$'chain:1.sens_sd',
           type='l',ylab="sens_sd",ylim =c(0,12),
           xlab="Iteraction")
      plot(par_ch$iter,par_ch$'chain:2.sens_sd',
           type='l',ylab="sens_sd",ylim =c(0,12),
           xlab="Iteraction")
      plot(par_ch$iter,par_ch$'chain:3.sens_sd',
           type='l',ylab="sens_sd",ylim =c(0,12),
           xlab="Iteraction")
      dev.off()
    }
    
    if( mod_ii == 2 ){
      tiff(paste0('results/simulations/chains/',
                  'expp/beta_',
                  mod_bet[ii],'_yr1_sd_',log_lambda_sd,'.tiff'),
           unit="in", width=6.3, height=8, res=600,compression="lzw" )
      par(mfcol=c(3,2))
      plot(par_ch$iter,par_ch$'chain:1.sens_mu',
           type='l',ylab="sens_mu",ylim =c(0,12),
           xlab="Iteraction")
      plot(par_ch$iter,par_ch$'chain:2.sens_mu',
           type='l',ylab="sens_mu",ylim =c(0,12),
           xlab="Iteraction")
      plot(par_ch$iter,par_ch$'chain:3.sens_mu',
           type='l',ylab="sens_mu",ylim =c(0,12),
           xlab="Iteraction")
      
      plot(par_ch$iter,par_ch$'chain:1.sens_sd',
           type='l',ylab="sens_sd",ylim =c(0,12),
           xlab="Iteraction")
      plot(par_ch$iter,par_ch$'chain:2.sens_sd',
           type='l',ylab="sens_sd",ylim =c(0,12),
           xlab="Iteraction")
      plot(par_ch$iter,par_ch$'chain:3.sens_sd',
           type='l',ylab="sens_sd",ylim =c(0,12),
           xlab="Iteraction")
      dev.off()
    }
    
    if( mod_ii == 3){
      tiff(paste0('results/simulations/chains/',
                  'gev/beta_',
                  mod_bet[ii],'_yr1_sd_',log_lambda_sd,'.tiff'),
           unit="in", width=6.3, height=8, res=600,compression="lzw" )
      par(mfcol=c(3,2))
      plot(par_ch$iter,par_ch$'chain:1.loc',
           type='l',ylab="loc",ylim =c(0,12),
           xlab="Iteraction")
      plot(par_ch$iter,par_ch$'chain:2.loc',
           type='l',ylab="loc",ylim =c(0,12),
           xlab="Iteraction")
      plot(par_ch$iter,par_ch$'chain:3.loc',
           type='l',ylab="loc",ylim =c(0,12),
           xlab="Iteraction")
      
      plot(par_ch$iter,par_ch$'chain:1.scale',
           type='l',ylab="scale",ylim =c(0,12),
           xlab="Iteraction")
      plot(par_ch$iter,par_ch$'chain:2.scale',
           type='l',ylab="scale",ylim =c(0,12),
           xlab="Iteraction")
      plot(par_ch$iter,par_ch$'chain:3.scale',
           type='l',ylab="scale",ylim =c(0,12),
           xlab="Iteraction")
      dev.off()
    }
    
    if( mod_ii == 4){
      tiff(paste0('results/simulations/chains/',
                  'sad/beta_',
                  mod_bet[ii],'_yr1_sd_',log_lambda_sd,'.tiff'),
           unit="in", width=6.3, height=8, res=600,compression="lzw" )
      par(mfcol=c(3,2))
      plot(par_ch$iter,par_ch$'chain:1.theta_m[5]',
           type='l',ylab="theta_m[5]",ylim =c(0,1),
           xlab="Iteraction")
      plot(par_ch$iter,par_ch$'chain:2.theta_m[5]',
           type='l',ylab="theta_m[5]",ylim =c(0,1),
           xlab="Iteraction")
      plot(par_ch$iter,par_ch$'chain:3.theta_m[5]',
           type='l',ylab="theta_m[5]",ylim =c(0,1),
           xlab="Iteraction")
      
      plot(par_ch$iter,par_ch$'chain:1.theta_m[12]',
           type='l',ylab="theta_m[12]",ylim =c(0,1),
           xlab="Iteraction")
      plot(par_ch$iter,par_ch$'chain:2.theta_m[12]',
           type='l',ylab="theta_m[12]",ylim =c(0,1),
           xlab="Iteraction")
      plot(par_ch$iter,par_ch$'chain:3.theta_m[12]',
           type='l',ylab="theta_m[12]",ylim =c(0,1),
           xlab="Iteraction")
      dev.off()
    }
    
  }
}
