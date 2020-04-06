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

# set log lambda sd
log_lambda_sd <- 0.3


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
mod_l <- list( '0.2'  = NULL,
               '0.45' = NULL,
               '0.7'  = NULL,
               '1.2'  = NULL,
               '1.7'  = NULL,
               '2.2'  = NULL )

# 'true' weight function. Mean in 5th month, sd = 1.
pdens  <- dnorm(1:12, 5, 1)
w_v    <- c( (pdens / sum(pdens)) * 0.2,
             (pdens / sum(pdens)) * 0.7,
             (pdens / sum(pdens)) * 0.1 )
# climate "data"
clim_m <- select(anom_a, -year) %>% as.matrix

# simulate data
sim_data <- function(beta_x){
  
  x      <- clim_m %*% w_v

  # function to produce normally distributed response
  set.seed(101)
  prod_y <- function(x) rnorm(length(x), mean=0+x*beta_x, 
                              sd = log_lambda_sd)
  
  # output data
  data.frame( x = x,
              y = prod_y(x),
              stringsAsFactors = F )

}

# plot it
plot(y ~ x, data=sim_data(1.2))


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
mod_gaus <- stan_model(file = paste0("code/stan/",family,"_gaus_nest.stan") )
mod_expp <- stan_model(file = paste0("code/stan/",family,"_expp_nest.stan") )
mod_gev  <- stan_model(file = paste0("code/stan/",family,"_gev_nest.stan")  )
mod_sad  <- stan_model(file = paste0("code/stan/",family,"_dirichlet_nest.stan") )


# fit models based on data
fit_mods <- function(df_x){
  
  # organize data into list to pass to stan
  dat_stan <- list(
    n_time  = nrow(clim_m),
    n_lag   = ncol(clim_m),
    y       = df_x$y,
    clim    = t(clim_m),
    clim1   = t(clim_m)[1:12 ,],
    clim2   = t(clim_m)[13:24,],
    clim3   = t(clim_m)[25:36,],
    M       = 12,    # number of months in a year
    K       = ncol(clim_m) / 12,
    expp_beta = 20
  )
  
  # gaussian moving window
  fit_gaus <- sampling(
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
  
  # exponential power moving window
  fit_expp <- sampling(
    object = mod_expp,
    data = dat_stan,
    pars = c('sens_mu', 'sens_sd', 'theta_y', 'alpha', 'beta', 'y_sd'),
    warmup = sim_pars$warmup,
    iter = sim_pars$iter,
    thin = sim_pars$thin,
    chains = sim_pars$chains,
    control = list(adapt_delta = 0.999)#,
    #control = list(adapt_delta = 0.999, stepsize = 0.001, max_treedepth = 20)
  )
  
  # Generalized Extreme Value nested
  fit_gev <- sampling(
    object = mod_gev,
    data = dat_stan,
    pars = c('loc', 'scale', "shape", 'theta_y', 'alpha', 'beta', 'y_sd'),
    warmup = sim_pars$warmup,
    iter = sim_pars$iter,
    thin = sim_pars$thin,
    chains = sim_pars$chains#,
    #control = list(adapt_delta = 0.999, stepsize = 0.001, max_treedepth = 20)
  )
  
  # Simplex nested
  fit_sad <- sampling(
    object = mod_sad,
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
        fit_sad )
  
}

model_names <- function(x){
  x %>% setNames( c('gaus','expp','gev','sad') )
}

mod_l$'0.2'  <- fit_mods( sim_data(0.2) ) %>% model_names
mod_l$'0.45' <- fit_mods( sim_data(0.45) ) %>% model_names
mod_l$'0.7'  <- fit_mods( sim_data(0.7) ) %>% model_names
mod_l$'1.2'  <- fit_mods( sim_data(1.2) ) %>% model_names
mod_l$'1.7'  <- fit_mods( sim_data(1.7) ) %>% model_names
mod_l$'2.2'  <- fit_mods( sim_data(2.2) ) %>% model_names

# # save this image
# save.image( paste0('results/simulations/',
#                    'nested_sd_',
#                    log_lambda_sd,'.Rdata') )

         
# model evaluation 

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
          paste0('results/simulations/divergent_sd_',
                 log_lambda_sd,'.csv'), row.names=F)
write.csv(rhat_df, 
          paste0('results/simulations/rhat_sd_',
                 log_lambda_sd,'.csv'), row.names=F)

# Bayesian checks -----------------------------------------------



# standardized betas --------------------------------------------


# choose which function to calculate x_antecedent with
x_ante_post_choose <- function(mod,post_m){
  
  # calculate antecedents ----------------------------------
  if(mod == 'sad'){
  
    x_ante_post <- function(pp,post_m){
      
      # x antecedent
      w_v    <- c(post_m[pp,paste0('theta_m.',1:12)]*post_m[pp,'theta_y.1'],
                  post_m[pp,paste0('theta_m.',1:12)]*post_m[pp,'theta_y.2'],
                  post_m[pp,paste0('theta_m.',1:12)]*post_m[pp,'theta_y.3']) %>% 
                  unlist
      # x_ante
      (clim_m %*% w_v) %>% as.vector
      
    }
    
  }
  
  if(mod == 'expp'){
  
    x_ante_post <- function(pp,post_m){
      
      # x antecedent
      expp_w_v <- dexppow(xx, post_m[pp,'sens_mu'], 
                              post_m[pp,'sens_sd'], 20)
      w_v <- c( ((expp_w_v / sum(expp_w_v))*post_m[pp,'theta_y.1']),
                ((expp_w_v / sum(expp_w_v))*post_m[pp,'theta_y.2']),
                ((expp_w_v / sum(expp_w_v))*post_m[pp,'theta_y.3']) )
     
      # x_ante
      (clim_m %*% w_v) %>% as.vector
      
    }
    
  }
  
  if(mod == 'gev'){
  
    x_ante_post <- function(pp,post_m){
      
      # gev
      gev_w_v  <- dgev(xx, post_m[pp,'loc'], 
                           post_m[pp,'scale'], 
                           post_m[pp,'shape'] )
      w_v      <- c( ((gev_w_v / sum(gev_w_v))*post_m[pp,'theta_y.1']),
                     ((gev_w_v / sum(gev_w_v))*post_m[pp,'theta_y.2']),
                     ((gev_w_v / sum(gev_w_v))*post_m[pp,'theta_y.3']) )
      
      (clim_m %*% w_v) %>% as.vector
      
    }
    
  }
  
  if(mod == 'gaus'){
  
    x_ante_post <- function(pp,post_m){
      
      # x antecedent
      gaus_w_v <- dnorm(xx, post_m[pp,'sens_mu'], 
                            post_m[pp,'sens_sd'])
      w_v <- c( ((gaus_w_v / sum(gaus_w_v))*post_m[pp,'theta_y.1']),
                ((gaus_w_v / sum(gaus_w_v))*post_m[pp,'theta_y.2']),
                ((gaus_w_v / sum(gaus_w_v))*post_m[pp,'theta_y.3']) )
     
      # x_ante
      (clim_m %*% w_v) %>% as.vector
      
    }
    
  }

  x_ante_post
  
}


# calculate standardized betas
st_beta_df <- expand.grid( b_input = c('0.2',0.45,0.7,1.2,1.7,2.2),
                           mod      = c("gaus","expp","gev","sad"),
                           stringsAsFactors = F )

beta_st <- function(ii,st_beta_df, mod_l){

  b_input <- st_beta_df$b_input[ii]
  mod      <- st_beta_df$mod[ii]

  # extract posterior
  post_m <- purrr::pluck(mod_l, b_input) %>% 
                purrr::pluck(mod) %>% 
                rstan::extract() %>% 
                as.data.frame 

  # produce x_antecedent
  x_ante_post <- x_ante_post_choose(mod, post_m)
  x_ante      <- lapply(1:nrow(post_m), x_ante_post, post_m)
  
  # calculate rest of the posterior
  post_pred <- function(pp,post_m){
    
    b01    <- post_m[pp,c('alpha','beta')] %>% as.numeric
    pl_df  <- data.frame( y      = sim_data(as.numeric(b_input))$y,
                          x      = x_ante[[pp]],
                          stringsAsFactors = F) %>% 
                mutate( y_pred = b01[1] + b01[2]*x ) %>% 
                mutate( perf   = y - y_pred )
      
    # beta standardized
    data.frame( b_st = b01[2] * sd(pl_df$x) / sd(pl_df$y),
                beta = b01[2],
                sd_x = sd(pl_df$x),
                sd_y = sd(pl_df$y) )
      
  }
  
  # get the posterior
  beta_post <- lapply(1:nrow(post_m), post_pred, post_m)
  # beta_post <- lapply(1:2, post_pred) 
  
  # spit the row out 
  beta_post %>% 
    bind_rows %>% 
    # remove NAs
    subset( !(is.na(b_st) | is.na(beta) | is.na(sd_x) | is.na(sd_y)) ) %>% 
    colMeans() %>% 
    as.data.frame %>% t %>% 
    as.data.frame %>% 
    mutate( model   = mod,
            b_sim   = b_input)
   
}

# exponential power distribution
dexppow <- function(x, mu, sigma, beta) {
    return((beta / (2 * sigma * gamma (1.0/beta)) ) *
             exp(-(abs(x - mu)/sigma)^beta));
}

xx <- 1:12

# use maximum amount of cores suggested (max is 20, suggested is 10)
cluster   <- parallel::makePSOCKcluster(4)

# attach packages that will be needed on each cluster
clusterEvalQ(cluster, list(library(lme4), library(dplyr), library(tidyr),
                           library(rstan), library(purrr), library(evd)) )

# attach objects that will be needed on each cluster
clusterExport(cluster, c('clim_m', 'mod_l', 'clim_m', 'st_beta_df',
                         'dexppow','x_ante_post_choose', 'sim_data',
                         'xx', 'w_v', 'log_lambda_sd') )

# try parallel run
beta_st_l  <- parLapply(cluster, 1:nrow(st_beta_df), beta_st, st_beta_df, mod_l)
beta_st_df <- beta_st_l %>% bind_rows

# store all this bounty
write.csv(beta_st_df,
          'results/simulations/beta_st_0.3_sim.csv',
          row.names=F)


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
                    beta_est  = beta_est )

# store figure
tiff(paste0('results/simulations/betas/beta_true_est_nest_sd_',
            log_lambda_sd,'.tiff'),
     unit="in", width=6.3, height=6.3, res=600,compression="lzw" )
par(mfrow=c(1,1), mar=c(3.5,3.5,0.2,0.2), 
    mgp=c(1.7,0.7,0))
plot(beta_est ~ jitter(beta_true), 
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
  ggsave(paste0('results/simulations/betas/beta_true_est_post_nest_sd_',
                 log_lambda_sd,'.tiff'),
         width = 6.3, height = 6.3)


# weights ------------------------------------------------------

# get weight params 
get_weight <- function(x){
  
  # get a parameter
  get_par <- function(x,param){
    x %>% 
     summary %>% 
    .$summary %>% 
    .[param,'mean']
  }

  gaus <- c(x[[1]] %>% get_par('sens_mu'),
            x[[1]] %>% get_par('sens_sd'),
            x[[1]] %>% get_par('theta_y[1]'),
            x[[1]] %>% get_par('theta_y[2]'),
            x[[1]] %>% get_par('theta_y[3]') )
  expp <- c(x[[2]] %>% get_par('sens_mu'),
            x[[2]] %>% get_par('sens_sd'),
            x[[1]] %>% get_par('theta_y[1]'),
            x[[1]] %>% get_par('theta_y[2]'),
            x[[1]] %>% get_par('theta_y[3]') )
  gev  <- c(x[[3]] %>% get_par('loc'),
            x[[3]] %>% get_par('scale'),
            x[[3]] %>% get_par('shape'),
            x[[3]] %>% get_par('theta_y[1]'),
            x[[3]] %>% get_par('theta_y[2]'),
            x[[3]] %>% get_par('theta_y[3]'))
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

# plot weights
plot_w <- function(weight_l){
  
  xx <- 1:12
  # gaus
  gaus_w_v <- dnorm(xx, weight_l$gaus[1], weight_l$gaus[2])
  gaus_w_v <- c( (gaus_w_v / sum(gaus_w_v)*weight_l$gaus[3]),
                 (gaus_w_v / sum(gaus_w_v)*weight_l$gaus[4]),
                 (gaus_w_v / sum(gaus_w_v)*weight_l$gaus[5]) )
  # expp
  expp_w_v <- dexppow(xx, weight_l$expp[1], weight_l$expp[2], 20)
  expp_w_v <- c( (expp_w_v / sum(expp_w_v)*weight_l$expp[3]),
                 (expp_w_v / sum(expp_w_v)*weight_l$expp[4]),
                 (expp_w_v / sum(expp_w_v)*weight_l$expp[5]) )
  # gev
  gev_w_v  <- dgev(xx,weight_l$gev[1], weight_l$gev[2], 
                      weight_l$gev[3])
  gev_w_v  <- c( (gev_w_v / sum(gev_w_v)*weight_l$gev[4]),
                 (gev_w_v / sum(gev_w_v)*weight_l$gev[5]),
                 (gev_w_v / sum(gev_w_v)*weight_l$gev[6]) )
  #sad
  sad_w_v  <- c(weight_l$sad[1:12]*weight_l$sad[13],
                weight_l$sad[1:12]*weight_l$sad[14],
                weight_l$sad[1:12]*weight_l$sad[15])
  
  plot(1:36, gaus_w_v,type='l',col='red',lwd=2,
       ylim=c(0,max(w_v)+0.01),ylab='Weight',xlab="Month")
  
  lines(1:36,expp_w_v,lwd=2,col='magenta')
  lines(1:36,gev_w_v,lwd=2,col='green')
  
  lines(1:36,sad_w_v,lwd=2,col='blue')
  lines(1:36,w_v,lty=2,lwd=3,col='black')
  
}


tiff(paste0('results/simulations/weights/month_w_nested_sd_',
            log_lambda_sd,'.tiff'),
     unit="in", width=6.3, height=8, res=600,compression="lzw" )
par(mfrow=c(3,2),
    mar=c(3,3.2,0.2,0.2), 
    mgp=c(1.7,0.3,0),
    cex.lab = 2)
plot_w(get_weight(mod_l$'0.2')) 
co <- par('usr')
text(28, co[4],pos=1,'beta=0.2', cex = 2)
plot_w(get_weight(mod_l$'0.45'))
text(28,co[4],pos=1,'beta=0.45', cex = 2)
plot_w(get_weight(mod_l$'0.7'))
text(28, co[4],pos=1,'beta=0.7', cex = 2)
plot_w(get_weight(mod_l$'1.2'))
text(8,co[4],pos=1,'beta=1.2', cex = 2)
legend('topright',
       c('gaus','expp','gev','sad','true'),
       col=c('red','magenta','green','blue','black'),
       bty='n', lwd=c(2,2,2,2,3),
                lty=c(1,1,1,1,2), cex=2)
plot_w(get_weight(mod_l$'1.7'))
text(28,co[4],pos=1,'beta=1.7', cex = 2)
plot_w(get_weight(mod_l$'2.2'))
text(28,co[4],pos=1,'beta=2.2', cex = 2)
dev.off()


# vecc <- dgev(1:36,5,1,0.4)
# lines(1:36,(vecc / sum(vecc)),col='magenta')


# examine chains ---------------------------------------------

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
                  mod_bet[ii],'_sd_',log_lambda_sd,'.tiff'),
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
                  mod_bet[ii],'_sd_',log_lambda_sd,'.tiff'),
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
                  mod_bet[ii],'_sd_',log_lambda_sd,'.tiff'),
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
                  mod_bet[ii],'_sd_',log_lambda_sd,'.tiff'),
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

# plot beta st -------------------------------------------------

beta_st_03 <- read.csv('results/simulations/beta_st_sim.csv') %>% mutate(sd_sim=0.3)
beta_st_05 <- read.csv('results/simulations/beta_st_0.5_sim.csv') %>% mutate(sd_sim=0.5)
beta_st_df <- bind_rows(beta_st_03, beta_st_05) %>% 
                mutate( sd_sim = as.character(sd_sim))


# beta_st and betas
ggplot(beta_st_df, aes(x=b_sim, y=b_st)) +
  geom_point( aes(shape=sd_sim, color=model),
              size = 3) +
  geom_abline( intercept=0, slope=1 ) +
  scale_colour_viridis( discrete = T ) +
  xlab( expression('Simulated '*beta) ) + 
  ylab( expression('Retrieved standardized '*beta) ) +
  geom_point( aes(x=b_sim, y=beta), shape=15 ) +
  ggsave('results/simulations/betas/beta_st.tiff',
          width = 6.3, height = 6.3, compression='lzw')

# residual sd
ggplot(beta_st_df) + 
  geom_point( aes(x=b_sim,y=sd_y,
                  color=model,
                  shape=sd_sim),
              size=3) +
  scale_colour_viridis( discrete = T) +
  ylab( 'sd of y (simulations)' ) +
  xlab( expression('Simulated '*beta) ) +
  ggsave('results/simulations/betas/sd_y.tiff',
         width = 6.3, height = 6.3, compression='lzw')

# variance of x_antecedent
ggplot(beta_st_df) + 
  geom_point( aes(x=b_sim,y=sd_x,
                  color=model,
                  shape=sd_sim),
              size=3) +
  scale_colour_viridis( discrete = T) +
  ylab( 'sd of x_antecedent (simulations)' ) +
  xlab( expression('Simulated '*beta) ) +
  ggsave(paste0('results/simulations/betas/sd_x.tiff'),
         width = 6.3, height = 6.3, compression='lzw') 
  

# examine chains using ShinyStan ---------------------------------------
launch_shinystan(mod_l$'0.2'[[3]])
