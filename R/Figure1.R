# Fig1.
# A all 36 months, yr2, gev, simplex, spline
# B year nested, gev_n, simplex_n
# C Spline
rm(list=ls())
setwd("C:/cloud/Dropbox/sApropos/")
source("C:/CODE/moving_windows/format_data.R")
library(tidyverse)
# library(dismo)
# library(mgcv)
library(testthat)
# library(rstan)
# library(loo)
# library(evd)
library(DirichletReg)
# library(gtools)
library(ggthemes)
library(gridExtra)


# FIGURE for first submission --------------------------------------------------

# normal
prod_norm <- function( m, s ){
  dnorm( 1:12, mean = m, sd = s ) / 
    sum( dnorm( 1:12, mean = m, sd = s) )
}

set.seed(1400)
offset  <- 0.25
dalph <- replace( rep(1,12), 3, 6)
w_csm <- data.frame( x     = 1:12,
                     w     = 1/12,
                     shape = 1,
                     model = 'CSM')
w_wmm <- data.frame( x     = 1:12 - offset,
                     w     = prod_norm(2, 2),
                     shape = 2,
                     model = 'WMM')
w_sam <- data.frame( x     = 1:12 + offset,
                     shape = 3,
                     w     = rdirichlet(1, dalph )[1,],
                     model = 'SAM')
w_mok <- data.frame( x     = 1,
                     w     = -5,
                     shape = 18,
                     model = 'FHM')
eff <- data.frame( x       = 1:12,
                   eff = rnorm(12, 0, 0.001) ) %>% 
  mutate( eff = replace( eff, 3, 0.25) ) %>% 
  mutate( eff = replace( eff, 10, -0.1) ) %>% 
  mutate( ymin = eff - runif( nrow(.), 0.1, 0.2),
          ymax = eff + runif( nrow(.), 0.1, 0.2) )


p1 <- list( w_csm,
            w_wmm,
            w_sam,
            w_mok) %>% 
  bind_rows %>% 
  mutate( shape = as.factor(shape) ) %>% 
  mutate( model = factor(model, levels = c('CSM','WMM','SAM','FHM')) ) %>% 
  mutate( ymin = w - runif( nrow(.), 0.01, 0.02),
          ymax = w + runif( nrow(.), 0.01, 0.02) ) %>% 
  mutate( ymin = replace(ymin, model == 'CSM', w[model == 'CSM']),
          ymax = replace(ymax, model == 'CSM', w[model == 'CSM']) ) %>% 
  ggplot() +
  geom_pointrange( aes(x, w,
                       ymin  = ymin,
                       ymax  = ymax,
                       group = model,
                       color = model,
                       shape = model) ) +
  scale_color_colorblind() +
  theme_minimal() +
  labs( x     = 'Month',
        y     = 'Weight',
        color = 'Model',
        pch   = 'Model') + 
  scale_x_continuous( breaks = seq(1,12,by=1) ) +
  scale_shape_manual( values = c(16,17,15,18) ) + 
  annotate("text", x=-1.7, y=0.3875,label="A)",size=4)+
  coord_cartesian( ylim=c(0, 0.375),
                   xlim=c(0, 13 ),
                   clip="off") +
  theme( axis.text.x = element_text( angle = 90,
                                     vjust = 0.5 ) )


p2 <- ggplot( eff ) +
  geom_pointrange( aes(x, eff,
                       ymin  = ymin,
                       ymax  = ymax),
                   shape = 18,
                   color = '#009E73') + 
  theme_minimal() + 
  labs( x     = 'Month',
        y     = 'Effect size' ) + 
  scale_x_continuous( breaks = seq(1,12,by=1) ) +
  annotate("text", x=-1.7, y=0.465,label="B)",size=4)+
  coord_cartesian( ylim=c(-0.3, 0.43),
                   xlim=c(0, 13 ),
                   clip="off" ) +
  theme( plot.margin = margin(t = 0, r = 72.5, b = 0, l = 0),
         axis.text.x = element_text( angle = 90,
                                     vjust = 0.5 ) )


ggsave( 'results/fig1_horiz.png', 
        plot = grid.arrange( p1, p2, ncol=2, widths=c(2,1.5)), 
        height = 2.5, width = 6.3)

ggsave( 'results/fig1_vertical.png', 
        plot   = grid.arrange( p1, p2, nrow=2, ncol = 1 ),
        height = 4.3, width = 3.15 )

ggsave( 'results/fig1A.png', 
        plot = p1, height = 2.5, width = 3.15)
ggsave( 'results/fig1B.png', 
        plot = p2, height = 2.5, width = 2.5)


# SIMPLIFIED FIGURE FIRST! ----------------------------

tiff( 'results/figure1.tiff',
      # width = 480, height = 480, units = "px", pointsize = 12,
      width = 6.3, height = 6.3, unit = 'in', res = 600,
      compression = 'lzw')

par(mfrow = c(1,1), mar = c(3.5,3.5,0.2,0.2), 
    mgp = c(2,0.6,0),
    oma = c(0,0,0,0) )

# A
plot(1:12,rep(1,12), type = 'n',
     ylim = c(0,0.2), 
     ylab = 'Weight', 
     xlab = 'Month', cex.lab = 2 )
abline( h = 0, lty=2)

# year 2 mean
segments(1,0,1,1/12, lwd = 4,     col = "#E69F00")
segments(1,1/12,12,1/12, lwd = 4, col = "#E69F00")
segments(12,1/12,12,0, lwd = 4,   col = "#E69F00")

# Dirichlet
set.seed(107)
lines(1:12, rdirichlet( 1, rep(3,12)),
      col = '#59B4E9', lwd=4)

# normal
prod_norm <- function( m, s ){
  dnorm( 1:12, mean = m, sd = s ) / 
    sum( dnorm( 1:12, mean = m, sd = s) )
}
lines( 1:12, prod_norm(2, 4), 
       lwd = 4, col = '#009E73' )

legend('topright', 
       legend = c('CSM',
                  'WMM',
                  'SAM'),
       cex = 1.5,
       lwd = 4, 
       bty = 'n',
       col = c('#E69F00','#009E73','#59B4E9'))

dev.off()




# set rstan options
rstan_options( auto_write = TRUE )
options( mc.cores = parallel::detectCores() )

# "pipeable" Reduce rbind
rbind_l <- function(x) Reduce(function(...) rbind(...), x)

# climate predictor, response, months back, max. number of knots
response  <- "grow"
clim_var  <- "precip"
m_back    <- 36    
st_dev    <- FALSE

# read data -----------------------------------------------------------------------------------------
lam       <- read.csv("all_demog_6tr.csv", stringsAsFactors = F)
m_info    <- read.csv("MatrixEndMonth_information.csv", stringsAsFactors = F)
clim      <- data.table::fread(paste0(clim_var,"_chelsa_hays.csv"),  stringsAsFactors = F)
# clim      <- data.table::fread(paste0(clim_var,"_fc_hays.csv"),  stringsAsFactors = F)

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
ii            <- 17
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
  clim_means = rowMeans(mod_data$climate),
  clim_yr = list( rowMeans(mod_data$climate[, 1:12]),
                  rowMeans(mod_data$climate[,13:24]),
                  rowMeans(mod_data$climate[,25:36]) ) %>% rbind_l,
  M       = 12,    # number of months in a year
  K       = ncol(mod_data$climate) / 12,
  S       = mod_data$resp$population %>% unique %>% length,
  site_i  = mod_data$resp$population %>% as.factor %>% as.numeric,
  expp_beta = expp_beta
)

# simulation parameters
sim_pars <- list(
  warmup = 1000, 
  iter = 4000, 
  thin = 2, 
  chains = 3
)

# year weights
fit_yr_weight <- stan(
  file = paste0("stan/",family,"_yr_dirichlet.stan"),
  data = dat_stan,
  pars = c('theta', 'alpha', 'beta', 'y_sd', 'log_lik'),
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains#,
  #control = list(adapt_delta = 0.999, stepsize = 0.001, max_treedepth = 20)
)

# Generalized extreme value 
fit_gev <- stan(
  file = paste0("stan/",family,"_gev.stan"),
  data = dat_stan,
  pars = c('loc', 'scale', "shape", 'alpha', 'beta', 'y_sd', 'log_lik'),
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains#,
  #control = list(adapt_delta = 0.999, stepsize = 0.001, max_treedepth = 20)
)

# Simplex - 24 months
dat_stan$clim <- t(mod_data$climate)
fit_24 <- stan(
  file = paste0("stan/",family,"_dirichlet.stan"),
  data = dat_stan,
  pars = c('theta', 'alpha', 'beta', 'y_sd', 'log_lik'),
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains#,
  #control = list(adapt_delta = 0.999, stepsize = 0.001, max_treedepth = 20)
)
dat_stan$clim <- mod_data$climate
 

# Nested models 
# update data list
dat_stan$clim         <- t(mod_data$climate)
dat_stan$clim1        <- t(mod_data$climate)[1:12 ,]
dat_stan$clim2        <- t(mod_data$climate)[13:24,]
dat_stan$clim3        <- t(mod_data$climate)[25:36,]

# Simplex nested
fit_24_nest <- stan(
  file = paste0("stan/",family,"_dirichlet_nest.stan"),
  data = dat_stan,
  pars = c('theta_y', 'theta_m', 'alpha', 'beta', 'y_sd', 'log_lik'),
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains#,
  #control = list(adapt_delta = 0.999, stepsize = 0.001, max_treedepth = 20)
)

# Generalized Extreme Value nested
fit_gev_nest <- stan(
  file = paste0("stan/",family,"_gev_nest.stan"),
  data = dat_stan,
  pars = c('loc', 'scale', "shape", 'theta_y', 'alpha', 'beta', 'y_sd', 'log_lik'),
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains#,
  #control = list(adapt_delta = 0.999, stepsize = 0.001, max_treedepth = 20)
)



# PARAMETRIC Models  ------------------------------------------------------------------
tiff(paste0("ms/methods/Models_Representations.tiff"),
     unit="in", width=6.3, height=6.3, res=600 ,compression="lzw")
tiff(paste0("C:/cloud/Dropbox/sApropos/conferences/ESA_symposium\presentation/Models_Representations.tiff"),
     unit="in", width=6.3, height=6.3, res=600 ,compression="lzw")


par(mfrow = c(2,2), mar = c(1.5,3.5,2,0.2), mgp = c(1.8,0.6,0), 
    oma = c(3,0,0,0), lwd=1 )

# A
plot(0:37,rep(1,38), type = 'n',
     ylim = c(0,0.16), main = "A. Non-nested models",
     ylab = 'Weight', 
     xlab = '' )
abline( h = 0, lty=2)

# three year mean
# segments(1,0,1,1/36, lwd = 3)
# segments(1,1/36,36,1/36, lwd = 3)
# segments(36,1/36,36,0, lwd = 3)

# year 2 mean
segments(13,0,13,1/12, lwd = 3, col = "blue")
segments(13,1/12,25,1/12, lwd = 3, col = "blue")
segments(25,1/12,25,0, lwd = 3, col = "blue")

# simplex
lines(1:36, summary(fit_24)$summary[,'mean'][paste0('theta[',1:36,']')], 
      lwd = 2, pch = 16, col = 'green')

# gev
gev_par   <- summary(fit_gev)$summary[,'mean'][c('loc','scale','shape')]
gev_w     <- dgev(1:36, loc = 5, scale = 5, shape = 0.2 )
gev_w_s   <- gev_w / sum(gev_w)
lines(1:36, gev_w_s, lwd = 2, col = 'red')

# legend
legend('topleft', c('Annual climate',
                    'Parametric moving window',
                    'Stochastic antecedent modeling'), cex = 0.9,
       lwd=2, lty=1, col = c('blue','red','green'), bty = 'n')


# B
plot(0:37,rep(1,38), type = 'n',
     ylim = c(0,0.16), 
     ylab = 'Weight', main = 'B. Nested models',
     xlab = '' )
abline( h = 0, lty=2)

# weighted year
yr_wgts   <- summary(fit_yr_weight)$summary[,'mean'][paste0('theta[',1:3,']')]*(1/12)


# year 2 mean
segments(1,0,1,yr_wgts[1], lwd = 3, col = "blue")
segments(1,yr_wgts[1],12.5,yr_wgts[1], lwd = 3, col = "blue")

segments(12.5,yr_wgts[1],12.5,yr_wgts[2], lwd = 3, col = "blue")
segments(12.5,yr_wgts[2],25.5,yr_wgts[2], lwd = 3, col = "blue")

segments(25.5,yr_wgts[2],25.5,yr_wgts[3], lwd = 3, col = "blue")
segments(25.5,yr_wgts[3],36,yr_wgts[3], lwd = 3, col = "blue")
segments(36,yr_wgts[3],36,0, lwd = 3, col = "blue")

# GEV nested
yr_wgts   <- yr_wgts*12
gev_w     <- dgev(1:12, loc = 5, scale = 5, shape = 0.2 )
gev_w_s   <- gev_w / sum(gev_w)
gev_yrs   <- summary(fit_gev_nest)$summary[,'mean'][paste0('theta_y[',1:3,']')]

lines(1:12, gev_w_s*yr_wgts[1] , lwd = 2, col = 'red')
lines(13:24, gev_w_s*yr_wgts[2] , lwd = 2, col = 'red')
lines(25:36, gev_w_s*yr_wgts[3] , lwd = 2, col = 'red')

# Simplex nested
simpl_m   <- summary(fit_24_nest)$summary[,'mean'][paste0('theta_m[',1:12,']')]
simpl_y   <- summary(fit_24_nest)$summary[,'mean'][paste0('theta_m[',1:3,']')]

lines(1:12,  simpl_m * yr_wgts[1], lwd = 2, col = 'green')
lines(13:24, simpl_m * yr_wgts[2], lwd = 2, col = 'green')
lines(25:36, simpl_m * yr_wgts[3], lwd = 2, col = 'green')



# climate predictor, months back, max. number of knots
clim_var  <- "precip"
gdd       <- T
m_back    <- 36    
knots     <- 9
response  <- 'log_lambda'

# read data -----------------------------------------------------------------------------------------
lam       <- read.csv("all_demog_6tr.csv", stringsAsFactors = F) 
m_info    <- read.csv("MatrixEndMonth_information.csv", stringsAsFactors = F)
clim      <- data.table::fread(paste0(clim_var,"_fc_hays.csv"),  stringsAsFactors = F)
spp       <- clim$species %>% unique


# add monthg info to lambda information
month_add <- m_info %>%
                mutate(SpeciesAuthor = trimws(SpeciesAuthor) ) %>%
                dplyr::select(SpeciesAuthor, MatrixEndMonth)
lam_add   <- subset(lam, SpeciesAuthor %in% month_add$SpeciesAuthor) %>%  
                dplyr::select(-MatrixEndMonth) %>%
                inner_join(month_add)
lam_min   <- subset(lam, SpeciesAuthor %in% setdiff(lam$SpeciesAuthor, m_info$SpeciesAuthor) )
lambdas   <- bind_rows(lam_min, lam_add) %>%
                subset( !is.na(MatrixEndMonth) ) %>%
                arrange( SpeciesAuthor )     

mod_sum_l <- list()
spp_list  <- lambdas$SpeciesAuthor %>% unique   
spp_list  <- spp_list[ -c(2,7,12,15,20,21,22,24) ] #,3,8,14,18,19,20)]


# NON PARAMETRIC models ------------------------------------------------------

ii            <- 3

# set species
spp_name      <- spp_list[ii]

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

# tests
expect_equal(length(spp_resp), length(clim_separate) )

# unique years, the basis of crossvalidation samples
unique_yr         <- mod_data$resp$year %>% unique

# format data for splines --------------------------------------------------------

# precipitation matrix
pmat      <- mod_data$climate %>% setNames(NULL) %>% as.matrix
pmean     <- apply(pmat, 1, mean, na.rm = T)

# set up lags
lags      <- matrix(0, nrow(mod_data$climate), ncol(mod_data$climate))
for(i in 1:ncol(lags)) lags[,i]=i
lagsm     <- as.matrix(lags)

# assemble final data frame
# pmat <- pmat/diff(range(pmat));
dat       <- cbind(mod_data$resp[,c('year','population',response)], pmat)
dat       <- as.data.frame(dat, colnames=TRUE)
dat$lags  <- lags
dat$pmat  <- pmat
dat$pmean <- pmean

# model fits -----------------------------------------------------------------------------------

# fit full model
if( length(unique(dat$population)) > 1 ){
  mod_full  <-gam(log_lambda ~ s(lags, k=knots, by=pmat, bs="cs"), # + population,
                  data=dat,method="GCV.Cp",gamma=1.4, na.action = na.omit)
}else{
  mod_full  <-gam(log_lambda ~ s(lags, k=knots, by=pmat, bs="cs"),
                  data=dat,method="GCV.Cp",gamma=1.4, na.action = na.omit)
}

# crossvalidation function
crxval_spline <- function(i, dat, pmat){

  # select excluded year
  yr          <- unique_yr[i]

  # train and test sets, null prediction
  train_set   <- subset( dat, year != yr )
  test_set    <- subset( dat, year == yr )

  # model: leave-one-year-out
  if( length(unique(dat$population)) > 1 ){
    mod_spline  <- gam(log_lambda ~ s(lags, k=knots, by=pmat, bs="cs",sp=mod_full$sp),# + population,
                       method="GCV.Cp",gamma=1.4, na.action = na.omit, data=train_set)
    mod_lm      <- lm(log_lambda ~ pmean, data=train_set)# + population
  }else{
    mod_spline  <- gam(log_lambda ~ s(lags, k=knots, by=pmat, bs="cs",sp=mod_full$sp),
                       method="GCV.Cp",gamma=1.4, na.action = na.omit, data=train_set)
    mod_lm      <- lm(log_lambda ~ pmean, data=train_set)
  }

  # predictions leave-one-out: spline, linear model
  pred_null   <- mean( train_set$log_lambda )
  pred_lm     <- predict(mod_lm,     newdata = test_set, type="response")
  pred_spline <- predict(mod_spline, newdata = test_set, type="response")

  # data frame of predictions and "design" (year-by-population)
  pred_df     <- data.frame(pred_null   = pred_null,
                            pred_lm     = pred_lm,
                            pred_spline = pred_spline,
                            stringsAsFactors = F
                            )
  design_pred <- dplyr::select(test_set, year, population)

  bind_cols(design_pred, pred_df)

}

# apply function across samples, and create data frame
crxval_l    <- lapply(1:length(unique_yr), crxval_spline, dat, pmat)
crxval_df   <- Reduce(function(...) rbind(...), crxval_l) %>%
                  left_join( dplyr::select(dat,year, population,log_lambda) ) # %>%
                  # subset( !is.na(pred_spline) )

# Mean squared error
dev0        <- calc.deviance(crxval_df$log_lambda, crxval_df$pred_null,
                             weights = rep(1, nrow(crxval_df) ),
                             family="gaussian", calc.mean = TRUE)
dev1        <- calc.deviance(crxval_df$log_lambda, crxval_df$pred_lm,
                             weights = rep(1, nrow(crxval_df) ),
                             family="gaussian", calc.mean = TRUE)
dev2        <- calc.deviance(crxval_df$log_lambda, crxval_df$pred_spline,
                             weights = rep(1, nrow(crxval_df) ),
                             family="gaussian", calc.mean = TRUE)

# plot results ---------------------------------------------------------------
plot(mod_full, 
     xlab = '', 
     ylab = 'Climate effect', 
     lwd = 2, main = 'C. Spline models', xaxt = 'n')
axis(1, at =c(0,10,20,30))

abline(h=0)

mtext('Months before demographic observation', 1, 2.5, at = 40, cex = 1.3)

dev.off()
