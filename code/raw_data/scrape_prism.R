# semi-automatic download of PRISM climate data
library(dplyr)
library(tidyr)
# devtools::install_github('jimhester/archive')
library(archive)
library(raster)
library(prism)
library(testthat)
library(RCurl)

# 1. Download PRISM MOHTHLY data for the coordinates 
# 2. Download PRISM YEARLY data for the coordinates 
# 3. Harmonize PRISM data with rest of CHELSA data

# read up locations
site_coord <- read.csv('data/coord_extra_climate.csv') %>% 
                dplyr::select(Lon, Lat) %>% 
                as.matrix

site_meta  <- read.csv('data/coord_extra_climate.csv',
                       stringsAsFactors = F) 


# PRISM set up -------------------------------------------------------------------------------

# what do I need from CHELSA?
prism_df <- expand.grid( variable  = c('ppt','tmean'),
                          year     = c(1981:2013),
                          month    = c(paste0('0',1:9),paste0(10:12)),
                          stringsAsFactors = F) %>% 
                arrange(variable,year,month)


# set up reading path
read_dir  <- 'ftp://prism.oregonstate.edu/monthly/'

# produce file name based on index 'ii'
produce_file_name <- function(ii){
  
  if( prism_df$variable[ii] == 'ppt'){
    file_root  <- paste0(prism_df$variable[ii],'/',prism_df$year[ii],
                         '/PRISM_',prism_df$variable[ii],'_stable_4kmM3_',
                          prism_df$year[ii],prism_df$month[ii],'_bil.zip')
  }
  
  if( prism_df$variable[ii] == 'tmean'){
    file_root  <- paste0(prism_df$variable[ii],'/',prism_df$year[ii],
                         '/PRISM_',prism_df$variable[ii],'_stable_4kmM2_',
                          prism_df$year[ii],prism_df$month[ii],'_bil.zip')
  }
  
  return(file_root)
  
}

# get all file names 
file_names <- lapply(1:nrow(prism_df), produce_file_name) %>% unlist

# get all the file links (from file name)
file_links <- paste0(read_dir,file_names)

# produce file destinations (put it all under C:/)
file_dest  <- gsub("tmean/[0-9]{4}/|ppt/[0-9]{4}/","",file_names) %>% 
                paste0('C:/',.)

file_links_good <- file_links


# Extract PRISM MONTHLY DATA ------------------------------------------------------------------------

# placeholder for year and monthly data
climate_all <- list()

# extract year and monthly data
for(ii in 1:length(file_links_good)){
  
  # extac with archive 
  # devtools::install_github('jimhester/archive')
  file_path <- file_links_good[ii]
  
  download.file( file_path, destfile = file_dest[ii], mode = "wb")
  archive_extract( archive(file_dest[ii]), "temp_dir")
  # # extract with 7z directly. This does extract directly in getwd()
  # system('"C:\\Program Files\\7-Zip\\7z" x "C:\\cloud\\MEGA\\Projects\\sApropos\\analyses\\CHELSA_temp_1979_01.7z"')
  
  # get climate information ----------------------------------------------------------------
  
  # read raster
  raster_file <- grep('.bil$',list.files('temp_dir'), value=T)
  rast_stack  <- raster(paste0('temp_dir/',raster_file) )
  
  # extract climatic info 
  values_clim <- raster::extract(rast_stack, site_coord,layer=1) #, method = 'bilinear')
  clim_df     <- data.frame( variable = prism_df$variable[ii],
                             year     = prism_df$year[ii],
                             month    = prism_df$month[ii],
                             value    = values_clim,
                             stringsAsFactors = F ) %>% 
                  bind_cols( site_meta )
  
  file.remove( paste0('temp_dir/',list.files('temp_dir/')) )
  file.remove( file_dest[ii] )
  
  climate_all[[ii]] <- clim_df
  
  print(ii)
  
}

# # extract year and monthly data
# # Service function to check that links work/exist
# check_links <- function(ii){
# 
#   url.exists(file_links[ii])
# 
# }
# # check that all file links actually exist
# check_links <- sapply(1:20,check_links)

climate_all <- lapply(1:length(file_links_good), 
                      extract_year_month)

# put it all in one data frame
climate_df <- climate_all %>% bind_rows()

# test that we have all month/variable combinations
test_df   <- climate_df %>% 
              dplyr::select(variable,year,month) %>% 
              mutate( year = as.integer(year) ) %>% 
              unique %>% 
              arrange( variable,year,month ) %>% 
              as.data.frame 

# tests
expect_true( all(prism_df$variable == test_df$variable) )
expect_true( all(prism_df$year == test_df$year) )
expect_true( all(prism_df$month == test_df$month) )

# store the data frame!
write.csv(climate_df,
         'data/climate/1_792_prism.csv',
         row.names=F)


# Download PRISM YEARLY data for the coordinates ---------------------------------------

# what do I need from CHELSA?
prism_df <- expand.grid( variable = c('ppt','tmean'),
                         year     = c(1979:1980),
                         stringsAsFactors = F) %>% 
                arrange(variable,year)

# set up reading path
read_dir  <- 'ftp://prism.oregonstate.edu/monthly/'

# create links 
prism_df  <- prism_df %>% 
               mutate( link = paste0(read_dir,
                                     variable,'/',year,'/',
                                     'PRISM_',variable,
                                     '_stable_4kmM2_',
                                     year,'_all_bil.zip') )


# extract year data (by month)
extract_year <- function(ii){
  
  # extac with archive 
  # devtools::install_github('jimhester/archive')
  file_path <- prism_df$link[ii]
  
  download.file( file_path, destfile = file_dest[ii], mode = "wb")
  archive_extract( archive(file_dest[ii]), "temp_dir")
  # # extract with 7z directly. This does extract directly in getwd()
  # system('"C:\\Program Files\\7-Zip\\7z" x "C:\\cloud\\MEGA\\Projects\\sApropos\\analyses\\CHELSA_temp_1979_01.7z"')
  
  # get climate information ----------------------------------------------------------------
  
  # extract yr/month from the string containing link
  extract_yr_month <- function(link_string, response){
    
    # digits
    digs <- regmatches(link_string, 
                       gregexpr("[[:digit:]]{6}", 
                       link_string) )
    
    if(response == 'year')  return( substr(digs,1,4) )
    if(response == 'month') return( substr(digs,5,6) )
    
  }
  
  # read raster
  raster_file <- grep('[0-9]{6}_bil.bil$',list.files('temp_dir'), value=T)
  
  # extract information from the 12 raster files
  extact_12_rast <- function( rast_file ){
    
    # get 
    rast_stack  <- raster(paste0('temp_dir/',rast_file) )
    
    # extract climatic info 
    values_clim <- raster::extract(rast_stack, site_coord,layer=1) #, method = 'bilinear')
    clim_df     <- data.frame( variable = prism_df$variable[ii],
                               year     = extract_yr_month(rast_file, 'year'),
                               month    = extract_yr_month(rast_file, 'month'),
                               value    = values_clim,
                               stringsAsFactors = F) %>% 
                      bind_cols( site_meta )
    clim_df
    
  }
  
  clim_df <- lapply(raster_file, extact_12_rast) %>% bind_rows
  
  file.remove( paste0('temp_dir/',list.files('temp_dir/')) )
  file.remove( file_dest[ii] )
  
  print(ii)
  
  return(clim_df)

}


# get all climate data for 1979/1980
all_1979_80 <- lapply(1:nrow(prism_df), extract_year) %>% 
                  bind_rows %>% 
                  mutate( year  = as.numeric(year),
                          month = as.numeric(month) )

# all 1981 on 
all_1981    <- read.csv('data/climate/1_792_prism.csv',
                        stringsAsFactors = F)

# combine datasets
all_ppt     <- bind_rows( subset(all_1979_80, variable == 'ppt'),
                          subset(all_1981,    variable == 'ppt') )

all_tmean   <- bind_rows( subset(all_1979_80, variable == 'tmean'),
                          subset(all_1981,    variable == 'tmean') )

all_prism   <- bind_rows( all_ppt, all_tmean )


# "split" astragalus scaphoides in two "species": one with all temporal reps.
# one with all spatial reps (but less temporal reps.)
# Astragalus_scaphoides: all sites (but not all years)
prism_astr  <- all_prism %>% 
                  mutate( SpeciesAuthor = replace(SpeciesAuthor,
                                                  SpeciesAuthor == 'Astragalus_scaphoides_6',
                                                  'Astragalus_scaphoides_6_site_rep') 
                         )

# Astragalus_scaphoides: all years (but not all sites)
astr_long   <- prism_astr %>% 
                  subset( MatrixPopulation %in% c('McDevitt Creek', 
                                                  'Sheep Corral Gulch') ) %>% 
                  mutate( SpeciesAuthor = 'Astragalus_scaphoides_6_long' )


prism_out   <- bind_rows(prism_astr, astr_long) %>% 
                  arrange(variable,year,month,SpeciesAuthor)

# read PRISM climate data
write.csv( prism_out,
           'data/climate/prism_astragalus_cryptantha_pediocactus.csv',
           row.names=F)



# Harmonize PRISM data with rest of CHELSA data ------------------------------

# read CHELSA climate data
chelsa_df  <- read.csv('data/airt_chelsa_hays.csv',
                       stringsAsFactors = F)  

# read PRISM climate data
prism_raw  <- read.csv('data/climate/prism_astragalus_cryptantha_pediocactus.csv',
                       stringsAsFactors = F)  

# prism data + metadata
prism_df   <- prism_raw %>% 
                rename( species    = SpeciesAuthor,
                        population = MatrixPopulation) %>% 
                spread( key = variable, value = value )


# format climate in sApropos format (namely, expand "by day") -----------


# current precipitation/temperature data
airt_df   <- data.table::fread('data/airt_chelsa_hays.csv')
prec_df   <- data.table::fread('data/precip_chelsa_hays.csv')

# format day one
day_one   <- as.Date( paste0("1/1/", 1979), format="%d/%m/%Y" ) 

# ~ total number of days in 1979-2013 range
tot_days  <- (2014-1979)*366

# Hays species (Dalgleish et al. 2010)
hays_spp  <- c("Cirsium_undulatum", "Echinacea_angustifolia", 
               "Hedyotis_nigricans","Lesquerella_ovalifolia", 
               "Paronychia_jamesii", "Psoralea_tenuiflora",      
               "Ratibida_columnifera", "Solidago_mollis", 
               "Sphaeralcea_coccinea", "Thelesperma_megapotamicum")

# climate frame (frame of climate variables)
clim_fr   <- as.Date(1:tot_days, day_one-1) %>%
                as.character %>%
                as.data.frame(stringsAsFactors=F) %>%
                separate_(col=".",into=c("year","month","day"),sep="-") %>%
                subset( year != '2014' ) %>% 
                mutate( year  = as.numeric(year),
                        month = as.numeric(month),
                        day   = as.numeric(day) ) %>% 
                # join with CHELSA data
                left_join( prism_df ) %>% 
                # remove Hays species (add weather data later)
                subset( !(species %in% hays_spp) ) %>% 
                dplyr::select( -Lat, -Lon ) %>% 
                mutate( julian = as.POSIXlt(paste(year,month,day, sep = "-"), 
                                            format = "%Y-%m-%d")$yday + 1 ) %>% 
                # format to harmonize with pre-existing data frame
                dplyr::select(-day) %>% 
                rename( airt = tmean,
                        day  = julian ) %>% 
                dplyr::select(species, population, year, day, ppt, airt) %>% 
                arrange(species, population, year, day) 


# harmonize prism with chelsa data
airt_out <- bind_rows( airt_df,  
                       dplyr::select(clim_fr, -ppt) ) %>% 
                arrange(species, population, year, day)
prec_out <- bind_rows( prec_df,  
                       dplyr::select(clim_fr, -airt) ) %>% 
                arrange(species, population, year, day)


# finally put it out
write.csv(airt_out,'data/airt_chelsa_prism_hays.csv',  row.names=F)
write.csv(prec_out,'data/precip_chelsa_prism_hays.csv',row.names=F)
