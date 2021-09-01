rm(list=ls())
source("R/format_data.R")
library(tidyverse)
library(mgcv)
library(testthat)
library(rstan)
library(loo)

# set up design data frame
design_df <- expand.grid( response = c('fec','grow','surv',
                                       'log_lambda'),
                          spp      = 1:39 )


# read data ----------------------------------------------------------------------------------------
lam       <- read.csv("data/all_demog_updt.csv", stringsAsFactors = F)
spp       <- lam$SpeciesAuthor %>% unique

# format data --------------------------------------------------------------------------------------

# calculate mean and variance for each dataset
all_vr_data <- function( ii ){
  
  print(ii)
  
  spp_i         <- design_df$spp[ ii ]
  resp          <- design_df$response[ ii ]
  spp_name      <- spp[spp_i]
  
  # lambda data
  spp_resp      <- format_species(spp_name, lam, resp) %>% bind_rows

  # put data out
  data.frame( vr    = resp,
              value = spp_resp[,3] )
  
}

# all vital rate data
all_data_df <- lapply( 1:nrow(design_df), all_vr_data) %>% bind_rows

write.csv( all_data_df,
           'results/prior_pc/all_vr_data.csv',
           row.names=F)

# Vr 
ggplot(all_data_df) +
  geom_histogram( aes(value) ) +
  facet_wrap( ~vr, scale = 'free' )

all_data_df %>% 
  subset( value < 10 ) %>% 
  ggplot() +
  geom_boxplot( aes(value) ) +
  facet_wrap( ~vr, scale = 'free' )


all_data_df %>% 
  subset( vr == 'fec' ) %>% 
  .$value %>% sd

log(0.55)
