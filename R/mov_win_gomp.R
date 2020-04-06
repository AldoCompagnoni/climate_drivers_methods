setwd("C:/Users/Aldo/Dropbox/sAPROPOS project/")
library(climwin)
library(dplyr)
library(tidyr)
library(magrittr)
library(glmnet)
#source("Analysis/functions_gomperz.R")

in_dir <- 'C:/Users/ac22qawo/Dropbox/popler/LTER/Data/Cedar Creek'

# read data ------------------------------------------------------------
pop   <- read.csv( paste0(in_dir,"/cdr_popler.csv"),
                  stringsAsFactors = F)
clm   <- read.csv( paste0(in_dir,"/DailyClimateSummary.csv"),
                  stringsAsFactors = F)
fire  <- list(year = c(1982:2004),
              fire =  c(2005:2014))

# format climate ---------------------------------------
clm   <- clm %>% mutate(Date = as.Date(Date,"%m/%d/%Y")) %>%
           mutate(Date = format.Date(Date, "%d/%m/%Y") ) %>%
           mutate(mean_t = (MaxTemp.degF. + MinTemp.degF.)/2 ) %>%
           rename(prec = Precip.inches.) %>%
           select(Date, mean_t, prec) 

# monthly data
month <-  clm %>%
            separate(Date, into = c("day", "month", "year"), sep="/") %>%
            group_by(year,month) %>%
            summarise(mean_t = mean(mean_t, na.rm=T),
                      prec = sum(prec, na.rm=T)) %>%
            mutate(month = as.numeric(month))

# month array
years <- pop$year %>% unique

# array for 36 months before each sampling date.
prec_form <- function(x){
  
  id <- which(month$year == x & month$month == 8)
  r  <- c(id:(id-35))
  return(month[r,"prec"])
  
}

# calculate monthly precipitation values
prec_mat  <- lapply(years, prec_form) %>% 
               do.call(cbind, .) %>% 
               setNames(years) %>% 
               t %>% 
               as.data.frame %>% 
               tibble::add_column(.before=1, year = row.names(.) ) %>% 
               mutate( year = as.numeric(year) )



# format pop. data --------------------------------------
pop   <- pop %>%
           select(year, genus, species, spatial_replication_level_1,
                  spatial_replication_level_2,treatment_type_1,
                  count_observation) %>% 
           subset(treatment_type_1 == 1 &
                  !(spatial_replication_level_1 %in% "site_cdr_exp001_field_D") &
                  !(species %in% "litter") )


# gomperz form
pop   <- pop %>% rename(Nt1 = count_observation)

# data at time step Nt0
pop_1 <- pop %>%
          as.data.frame() %>%
          mutate(year = year + 1) %>% 
          rename( Nt0 = Nt1 )

# data at time step Nt1
pop_2 <- pop

# Gompertz model
gomp  <- merge(pop_1, pop_2) %>% 
          mutate(Xt0 = log(Nt0),
                 Xt1 = log(Nt1)) %>% 
          # fire
          mutate(fire = year) %>% 
          mutate(fire = replace(year, year < 2005, 0) ) %>% 
          mutate(fire = replace(fire, year>2004 ,1) ) %>% 
          left_join(prec_mat) %>% 
          # select Poa
          # subset( genus == 'Poa' & species == 'pratensis' ) %>% 
          # collapse two replication levels
          mutate( site = paste0(spatial_replication_level_1,
                                spatial_replication_level_2) ) %>% 
          select( -spatial_replication_level_1,
                  -spatial_replication_level_2 ) %>% 
          # let's forget about fire years for now
          subset( year < 2005 )
  
# isolate climate in its own data frame
clim_df <- gomp %>% select(V1:V36)

# glmnet -----------------------------------------------------------------

# prepare data
preds <- gomp %>% 
          select(V1:V24) %>% 
          apply(2, scale ) %>% 
          as.matrix
y_v   <- gomp$Xt1-gomp$Xt0

# crossvalidation
cvfit <- cv.glmnet(x = preds, y = y_v, 
                   family = "gaussian",n=15,alpha=c(1),keep=T)

beta_est <- coef(cvfit, s = cvfit$lambda.min)[,1][-1]

plot(beta_est, type='b')
plot(gomp$Nt1 ~ gomp$Nt0)



# stan models ------------------------------------------------------------

dat_stan <- list(
  n_time  = nrow(clim_df),
  n_lag   = ncol(clim_df),
  y       = gomp$Xt1,
  x0      = gomp$Xt0,
  # clim    = t(clim_df),
  clim1   = t(clim_df)[1:12 ,],
  clim2   = t(clim_df)[13:24,],
  clim3   = t(clim_df)[25:36,],
  M       = 12,    # number of months in a year
  K       = ncol(mod_data$climate) / 12,
  expp_beta = expp_beta
)


# Simplex nested
fit_36_nest <- stan(
  file = paste0("code/stan/gomp/normal_dirichlet_nest.stan"),
  data = dat_stan,
  pars = c('theta_y', 'theta_m', 'alpha', 'beta', 'dd', 'y_sd'),
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = 1 #sim_pars$chains
)


fit_36_nest %>% 
  rstan::extract() %>% 
  as.data.frame %>% 
  select(theta_y.1:theta_y.3) %>% 
  boxplot

fit_36_nest %>% 
  rstan::extract() %>% 
  as.data.frame %>% 
  select(theta_m.1:theta_m.12) %>% 
  boxplot

fit_36_nest %>% 
  rstan::extract() %>% 
  as.data.frame %>% 
  .$beta %>% hist(freq=F)

dens <- fit_36_nest %>% 
  rstan::extract() %>% 
  as.data.frame %>% 
  .$beta %>% density

lines(dens$x, dens$y, lwd=2)
