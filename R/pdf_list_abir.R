# Create list of DOIs for abir to decide which is a fixed and what is a random effect
# 16.9.2019
rm(list=ls())
library(dplyr)

# set rstan options
rstan_options( auto_write = TRUE )
options( mc.cores = parallel::detectCores() )

# read data -----------------------------------------------------------------------------------------
lam       <- read.csv("data/all_demog_updt.csv", stringsAsFactors = F)

# put out the list
lam %>% 
  select(SpeciesAuthor, DOI.ISBN) %>% 
  unique %>% 
  write.csv('results/issues/dois_species_methods.csv', 
            row.names=F )

  