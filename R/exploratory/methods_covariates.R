# bjtwd<-"C:/Users/admin_bjt162/Dropbox/A.Current/Ongoing_Collab_Research/sApropos project/"
rm(list=ls())
source("code/format_data.R")
library(dplyr)
library(tidyr)
library(mgcv)
library(testthat)
library(rstan)
library(loo)
options( stringsAsFactors = F )

# methods species
dir_synth  <- 'C:/CODE/synthesis_sapropos/'
method_df  <- read.csv("data/all_demog_updt.csv", stringsAsFactors = F)
synth_df   <- read.csv( paste0(dir_synth,'data/vitalRatesCollapsed.csv') )
method_spp <- method_df$SpeciesAuthor %>% unique
synth_spp  <- synth_df$SpeciesAuthor %>% unique

# Only differences is Astragalue scaphoides_6, which I separated
intersect(method_spp, synth_spp)
setdiff(method_spp, synth_spp)

# what studies need a co-variate
cov_df    <- read.csv( paste0(dir_synth,'data/species_author_compadre.csv') ) %>% 
               subset( !is.na( type ) ) %>% 
               subset( !(type == '') )

# Ciao
setdiff(   method_spp, cov_df$SpeciesAuthor )
intersect( method_spp, cov_df$SpeciesAuthor )
