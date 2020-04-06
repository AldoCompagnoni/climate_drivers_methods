# Summary figures for the maniscript
rm(list=ls())
source("code/format_data.R")
library(tidyverse)
library(testthat)

# read data 
lam       <- read.csv("data/all_demog_updt.csv", stringsAsFactors = F)

# number of studies
n_stud    <- lam %>% 
                subset( !is.na(DOI.ISBN) ) %>% 
                dplyr::select(DOI.ISBN) %>% 
                unique %>% 
                nrow %>% 
                # +4 (studies witn DOI==NA) -1 (Astragalus_scaphoides_2 is split in 2 datsets)
                `+` (3)
  
# number of species 
n_spp     <- lam$SpeciesAccepted %>% unique %>% length

# number of species with spatial replication
n_spp_spat<- lam %>% 
                dplyr::select( SpeciesAccepted, MatrixPopulation) %>% 
                unique %>% 
                count( SpeciesAccepted ) %>% 
                subset( n > 1) %>% nrow

# number of populations
n_pop     <- lam %>% 
                select( SpeciesAccepted, MatrixPopulation ) %>% 
                unique %>% 
                nrow

# number of populations per species
pop_p_spp <- lam %>% 
               select( SpeciesAccepted, MatrixPopulation ) %>% 
               unique %>% 
               count( SpeciesAccepted )

# studies preceeding 1979 (not using CHELSA data)
n_pre_79  <- subset(lam, grepl('Dalgleish', Authors) ) %>% 
               select(SpeciesAuthor,Authors) %>% 
               unique %>% 
               nrow

# n of CHELSA species and CHELSA populations 
# (spp==pops, because Dalgleish et al. only had 1 population per species)
n_ch_pop   <- n_pop - n_pre_79
n_ch_spp   <- n_spp - n_pre_79


# average number of datapoints ("y") per dataset
n_data_p   <- count(lam, SpeciesAuthor)  

n_data_p$n %>% mean
n_data_p$n %>% median


# summary table on data used
spp_tab <- lam %>% 
              dplyr::select(SpeciesAuthor, Authors, Journal, 
                            YearPublication, DOI.ISBN) %>% 
              unique 

pop_n   <- lam %>% 
              dplyr::select(SpeciesAuthor, MatrixPopulation) %>% 
              unique %>% 
              count(SpeciesAuthor) %>% 
              rename( n_of_populations = n )
    
year_n  <- lam %>% 
              dplyr::select(SpeciesAuthor, MatrixStartYear) %>% 
              unique %>% 
              count(SpeciesAuthor) %>% 
              rename( n_of_years = n )

# format and store table
Reduce(function(...) full_join(...), 
       list(spp_tab, pop_n, year_n)) %>% 
  dplyr::select(SpeciesAuthor, n_of_populations, n_of_years, 
                Authors, Journal, YearPublication, DOI.ISBN) %>% 
  mutate( Authors = replace( Authors,
                             SpeciesAuthor == 'Cryptantha_flava_2',
                             'Salguero-Gómez' ) ) %>% 
  write.csv('results/summaries_ms/species_tab.csv',
            row.names=F)
