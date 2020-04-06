#bjtwd<-"C:/Users/admin_bjt162/Dropbox/A.Current/Ongoing_Collab_Research/sApropos project/"
rm(list=ls())
setwd("C:/cloud/Dropbox/sApropos/")
source("C:/CODE/moving_windows/format_data.R")
library(tidyverse)
library(dismo)
library(mgcv)
library(testthat)
library(rstan)
library(loo)
library(egg)

# set rstan options
rstan_options( auto_write = TRUE )
options( mc.cores = parallel::detectCores() )

# "pipeable" Reduce rbind
rbind_l <- function(x) Reduce(function(...) rbind(...), x)

# climate predictor, response, months back, max. number of knots
response  <- "fec"
clim_var  <- "airt"
m_back    <- 36    st_dev    <- FALSE

# read data -----------------------------------------------------------------------------------------
lam       <- read.csv("all_demog_6tr.csv", stringsAsFactors = F)
m_info    <- read.csv("MatrixEndMonth_information.csv", stringsAsFactors = F)
clim      <- data.table::fread(paste0(clim_var,"_chelsa_hays.csv"),  stringsAsFactors = F)
# clim      <- data.table::fread(paste0(clim_var,"_fc_hays.csv"),  stringsAsFactors = F)

dalgleish_spp <- c("Cirsium_undulatum", "Echinacea_angustifolia", 
                   "Hedyotis_nigricans","Lesquerella_ovalifolia", 
                   "Paronychia_jamesii", "Psoralea_tenuiflora",      
                   "Ratibida_columnifera", "Solidago_mollis", 
                   "Sphaeralcea_coccinea", "Thelesperma_megapotamicum") %>% 
                    paste0(collapse='|')


spp       <- lam$SpeciesAuthor %>% unique

# format data --------------------------------------------------------------------------------------

# set up model "family" based on response
if( response == "surv" | response == "grow" )             family = "beta" 
if( response == "fec" )                                   family = "gamma"
if( grepl("PreRep", response) | grepl("Rep", response) )  family = "beta"
if( response == "rho" | response == "react_fsa" )         family = "gamma"
if( response == "log_lambda")                             family = "normal"

expp_beta     <- 20


# series of response variable
resp_series   <- function(ii, response){
  
  # set species (I pick Sphaeraclea_coccinea)
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

  print(spp_name)
  
  mod_data$resp %>% 
    mutate( species = spp_name ) %>% 
    return
  
}


spp <- setdiff(spp, grep('Daphne',spp,value=T) )

# create series of 
all_surv  <- lapply(1:34, resp_series, 'surv') %>% rbind_l
all_grow  <- lapply(1:34, resp_series, 'grow') %>% rbind_l
all_fec   <- lapply(1:34, resp_series, 'fec') %>% rbind_l
all_lam   <- lapply(1:34, resp_series, 'log_lambda') %>% rbind_l
all_df    <- Reduce(function(...) full_join(...), 
                    list(all_surv, all_grow, all_fec, all_lam) ) %>% 
                dplyr::select(-year, -month) %>% 
                group_by(species) %>% 
                summarize_all( sd ) %>%
                ungroup %>% 
                mutate( species = as.factor(species) ) %>% 
                mutate( col = 'Other spp.',
                        siz = 0.5 ) %>%
                mutate( col = replace( col, 
                                       grepl(dalgleish_spp, species),
                                       'Adler spp.') ) 

# plots of standard deviations
p1 <- ggplot(all_df, aes(species, surv, col)) +
        geom_point( aes(color = col), size = 3 ) +
        ylab( "Survival (Standard deviation)") +
        xlab('') +
        ggtitle('Survival') +
        theme( axis.text.x  = element_blank(),
               plot.title   = element_text(angle = 0, hjust = 0.5, size = 15) )

p2 <- ggplot(all_df, aes(species, grow, col)) +
        geom_point( aes(color = col), size = 3) +
        ylab( "Growth (Standard deviation)") +
        xlab('') +
        ggtitle('Growth') +
        theme( axis.text.x  = element_blank(),
               plot.title   = element_text(angle = 0, hjust = 0.5, size = 15) )

p3 <- ggplot(all_df, aes(species, fec, col)) +
        geom_point( aes(color = col), size = 3) +
        ylab( "Fecundity (Standard deviation)") +
        xlab('') +
        ggtitle('Fecundity') +
        theme( axis.text.x  = element_blank(),
               plot.title   = element_text(angle = 0, hjust = 0.5, size = 15) )

p4 <- ggplot(all_df, aes(species, log_lambda, col)) +
        geom_point( aes(color = col), size = 3) +
        ylab( "Log Lambda (Standard deviation)") +
        xlab('') +
        ggtitle('Log Lambda') +
        theme( axis.text.x  = element_blank(),
               plot.title   = element_text(angle = 0, hjust = 0.5, size = 15) )

tiff( 'vr_st_dev.tiff',
      unit="in", width=6.3, height=6.3, res=600,compression="lzw" )

grid.arrange(nrow = 2, ncol = 2,
             grobs = list(p1, p2, p3, p4)
             )

dev.off()

# Eryngium_cuneifolium ----------------------------------------------
cuneif <- Reduce(function(...) full_join(...), 
            list(all_surv, all_grow, all_fec, all_lam) ) %>% 
            subset( grepl('Eryngium_cuneifolium',species) )

p1 <- ggplot(cuneif, aes(x=year,y=surv,group=population,colour=population)) +
        ggtitle('Survival') +
        geom_point()+
        geom_line(aes(lty=population)) +
        theme( plot.title   = element_text(angle = 0, hjust = 0.5, size = 15) )
        
p2 <- ggplot(cuneif, aes(x=year,y=grow,group=population,colour=population)) +
        ggtitle('Growth') +
        geom_point()+
        geom_line(aes(lty=population)) +
        theme( plot.title   = element_text(angle = 0, hjust = 0.5, size = 15) )

p3 <- ggplot(cuneif, aes(x=year,y=fec,group=population,colour=population)) +
        ggtitle('Fecundity') +
        geom_point()+
        geom_line(aes(lty=population)) +
        theme( plot.title   = element_text(angle = 0, hjust = 0.5, size = 15) )

p4 <- ggplot(cuneif, aes(x=year,y=log_lambda,group=population,colour=population)) +
        ggtitle('Log Lambda') +
        geom_point() +
        geom_line(aes(lty=population)) +
        theme( plot.title   = element_text(angle = 0, hjust = 0.5, size = 15) )

tiff( 'vr_Eryngium_cuneifolium.tiff',
      unit="in", width=6.3, height=6.3, res=600,compression="lzw" )
grid.arrange(nrow = 2, ncol = 2,
             grobs = list(p1, p2, p3, p4)
             )
dev.off()
