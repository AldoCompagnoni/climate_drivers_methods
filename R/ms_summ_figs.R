# Summary figures for the maniscript
rm(list=ls())
source("R/format_data.R")
library(tidyverse)
library(testthat)

# read data 
lam       <- read.csv("data/all_demog_updt.csv", stringsAsFactors = F) %>% 
              subset( !(SpeciesAuthor %in% c( 'Astragalus_scaphoides_6_site_rep',
                                              'Astragalus_scaphoides_2') ) 
                    ) %>% 
              mutate( SpeciesAccepted = replace( SpeciesAccepted,
                                                 is.na(SpeciesAccepted),
                                                 'Cryptantha flava') )

# total number of observations (for appendix on Prior PC)
expect_equal( lam$log_lambda %>% is.na %>% sum, 0 )
expect_equal( lam$surv %>% is.na %>% sum, 0 )
expect_equal( lam$grow %>% is.na %>% sum, 0 )
expect_equal( lam$fec %>% is.na %>% sum, 0 )
n_y_norm   <- lam$log_lambda %>% length
n_y_beta   <- (lam$surv %>% length) + (lam$grow %>% length)
n_y_gamma  <- lam$fec %>% length

# number of studies
n_stud    <- lam %>% 
                subset( !is.na(DOI.ISBN) ) %>% 
                dplyr::select(DOI.ISBN) %>% 
                unique %>% 
                nrow %>% 
                # +4 (studies witn DOI==NA) -1 (Astragalus_scaphoides_2 is split in 2 datsets)
                `+` (3)
  
# number of species 
n_spp     <- lam$SpeciesAccepted %>% 
                unique %>% 
                sort %>% 
                length

# number of studies 
n_studies <- lam$SpeciesAuthor %>% 
              unique %>% 
              sort %>% 
              length

# number of populations with spatial replication
n_spp_spat<- lam %>% 
                dplyr::select( SpeciesAuthor, MatrixPopulation) %>% 
                unique %>% 
                count( SpeciesAuthor ) %>% 
                subset( n > 1) %>% 
                nrow

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

# studies preceeding 1979 (not using CHELSA data)
n_pre_79_2<- lam %>% 
              subset( MatrixStartYear <= 1979) %>% 
              select( SpeciesAuthor, MatrixPopulation ) %>% 
              unique %>% 
              nrow

# Table of before 1979 vs. after 1979 
lam %>% 
  mutate( climate_data = 'CHELSA data' ) %>% 
  mutate( climate_data = replace( climate_data, 
                                  MatrixStartYear <= 1979,
                                  'Original study data') ) %>% 
  mutate( species = SpeciesAuthor ) %>% 
  mutate( species = gsub('_',' ',species) ) %>% 
  mutate( species = replace(species,
                            grepl('Astragalus scaphoides 6 long',species),
                            'Astragalus scaphoides') ) %>% 
  mutate( species = replace(species,
                            grepl('Cirsium pitcheri 8',species),
                            'Cirsium pitcheri (1)') ) %>% 
  mutate( species = replace(species,
                            grepl('Cirsium pitcheri 4',species),
                            'Cirsium pitcheri (2)') ) %>% 
  mutate( species = replace(species,
                            grepl('Cirsium pitcheri 6',species),
                            'Cirsium pitcheri (3)') ) %>% 
  mutate( species = gsub(' 2','',species) ) %>% 
  dplyr::select( species, climate_data ) %>% 
  unique %>% 
  write.csv( 'results/climate_data_origin.csv',
             row.names = F )
  
  
# min max years


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


# calculate the correlation of climate for multiple-site species ---------------

clim_var  <- 'precip'
m_back    <- 36


# names of species with multiple populations
spat_spp <- lam %>% 
              dplyr::select( SpeciesAuthor, MatrixPopulation) %>% 
              unique %>% 
              count( SpeciesAuthor ) %>% 
              subset( n > 1) %>% 
              .$SpeciesAuthor %>% 
              # remove Daphne rodriguezii for now
              Filter( function(x) x != "Daphne_rodriguezii", .)

# calculate correlations 
calc_sites_correlation <- function( ii, clim_var ){

  print(ii)
  
  # which species do we care about?
  spp_name      <- spat_spp[ii]
  
  # download climate
  clim <- data.table::fread(paste0('data/',
                                   clim_var,
                                   "_chelsa_prism_hays_2014.csv"),
                            stringsAsFactors = F)
  
  
  # lambda data
  spp_resp      <- format_species(spp_name, lam, 'log_lambda')
  
  # climate data
  clim_separate <- clim_list(spp_name, clim, spp_resp)
  clim_detrnded <- lapply(clim_separate, clim_detrend, clim_var)
  clim_mats     <- Map(clim_long, clim_detrnded, spp_resp, m_back)
  
  # model data
  mod_data      <- lambda_plus_clim(spp_resp, clim_mats, 'log_lambda')
  
  # site level climate
  site_clim     <- cbind( mod_data$resp, mod_data$climate ) %>% 
                    select( -month, -log_lambda ) %>% 
                    gather( month, value, V1:V36 ) %>% 
                    pivot_wider( names_from  = population,
                                 values_from = value ) %>% 
                    setNames( c('year','month',
                                paste0('site_',1:(ncol(.)-2)) ) 
                              ) %>% 
                    select( -year ) %>% 
                    mutate( month = as.factor(month) )
  
  # calculate correlation
  correlate <- function( x ){
    tmp <- cor(x)
    tmp[upper.tri(tmp)]
  } 
  
  # correlation values
  cor_v <- split(site_clim[-1], site_clim$month ) %>% 
            lapply( as.matrix ) %>% 
            lapply( correlate ) %>% 
            unlist

  data.frame( spp  = spp_name,
              corr = cor_v )
  
}

corr_prec <- lapply( 1:16, calc_sites_correlation, 'precip')
corr_airt <- lapply( 1:16, calc_sites_correlation, 'airt')

corr_prec_df <- corr_prec %>% bind_rows 
corr_airt_df <- corr_airt %>% bind_rows 

corr_prec_df$corr %>% mean(na.rm=T)
corr_prec_df$corr %>% median(na.rm=T)
corr_prec_df$corr %>% quantile(na.rm=T, prob = .05)

corr_airt_df$corr %>% mean(na.rm=T)
corr_airt_df$corr %>% median(na.rm=T)
corr_airt_df$corr %>% quantile(na.rm=T, prob = .05)
