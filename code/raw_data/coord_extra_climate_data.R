library(measurements)

# read in compadre and Astragalus data
compadre <- readRDS('data/compadre_unrel_v4.rds')
astraga  <- readRDS('data/Astragalus_scaphoides_6.rds')
crypt    <- read.csv('data/raw_mats/Cryptantha_flava.csv')

# convert lat/lon in decimal form
conv_plot_coord <- function(lat_in, lon_in, from_unit){
  
  coord_df <- data.frame( lat = conv_unit(lat_in,  
                                          from = from_unit, 
                                          to = 'dec_deg'),
                          lon = conv_unit(lon_in, 
                                          from = from_unit, 
                                          to = 'dec_deg'),
                          stringsAsFactors = F) %>% 
                mutate(   lat = as.numeric(lat),
                          lon = as.numeric(lon) )
  
  return(coord_df)
  
}

# cryptantha
crypta_lat_lon <- compadre$metadata %>% 
                    subset( grepl('Cryptantha',SpeciesAuthor) ) %>% 
                    dplyr::select( SpeciesAuthor, MatrixPopulation, Lat, Lon ) %>% 
                    unique %>% 
                    mutate( SpeciesAuthor = 'Cryptantha_flava_2')

# Astragalus
astr_lat_lon   <- astraga$metadata %>% 
                    dplyr::select(SpeciesAuthor, MatrixPopulation, Lat, Lon ) %>% 
                    unique()


# periocactus 
pedio  <- compadre$metadata %>% 
            subset( SpeciesAuthor == 'Pediocactus_bradyi' ) %>% 
            dplyr::select(SpeciesAuthor, MatrixPopulation, MatrixStartYear,
                   MatrixComposite, Lat, Lon ) %>% 
            mutate( Lat = replace(Lat,
                                  SpeciesAuthor == 'Pediocactus_bradyi',
                                  conv_plot_coord('36 48 56', 
                                                  '-111 38 16', 
                                                  'deg_min_sec')$lat),
                    Lon = replace(Lon,
                                  SpeciesAuthor == 'Pediocactus_bradyi',
                                  conv_plot_coord('36 48 56', 
                                                  '-111 38 16', 
                                                  'deg_min_sec')$lon)
                    )

pedio_lat_lon <- pedio %>% 
                   subset( MatrixComposite == 'Individual' ) %>% 
                   dplyr::select(SpeciesAuthor, MatrixPopulation, Lat, Lon) %>% 
                   unique


# lat lon for which we need extra climatic data
list( crypta_lat_lon,
      astr_lat_lon,
      pedio_lat_lon) %>% 
  bind_rows %>% 
  write.csv('data/coord_extra_climate.csv', 
            row.names=F)
