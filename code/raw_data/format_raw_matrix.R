library(readxl)
library(dplyr)
library(tidyr)
library(measurements)

# 1. get rows of separate matrices
# 2. get matrices matA, matU, matF, matC
# 3. get metadata
  #3a. convert Lat/Lon to DECIMAL.
  # IMPORTANT: need to understand whether Lon is negative or positive
# 4. get matrixClass
# 5. output a compadre-like object

# prove species name (SpeciesAuthor)
spp_nam <- 'Astragalus_scaphoides_6'

# raw matrix (these excel files should be in "Data" somewhere)
raw_mat <- read_xlsx('data/raw_mats/Astragalus_scaphoides_6.xlsx',
                     sheet='MatrixStacked')
raw_spp <- read_xlsx('data/raw_mats/Astragalus_scaphoides_6.xlsx',
                     sheet='SpeciesDescriptors')


# get rows of separate matrices ---------------------------------

# get indices of first row
first_row <- setdiff( c(1:nrow(raw_mat)),
                        which(raw_mat$EnteredBy == 'NA') )

# matrix dimension
if( is.na(first_row[2]) ){
  # if we only have 1 matrix
  mat_dim <- nrow(raw_mat)
}else{
  mat_dim   <- (first_row[2]-1)
}

# indices raws associated w/ separate matrices
mat_r_ids  <- lapply(first_row, function(x, add_rows) x:(x+add_rows),
                     mat_dim-1)


# get matrices matA, matU, matF, matC -------------------------------

# where do matrices start in the spreadsheet?
first_col <- which(names(raw_mat) == 'MatrixClassNumber') + 1

# indices to add to get to "next first column"
add_cols  <- mat_dim+1

# indices for columns associatedc w/ matA, matU, matF, and matC -----
first_cols<- first_col + c(0, add_cols, add_cols*2, add_cols*3)

# indices raws associated w/ separate matrices
mat_c_ids  <- lapply(first_cols, function(x, add_cols) x:(x+add_cols),
                     mat_dim-1) %>%
                # precise
                setNames( c('matA', 'matU', 'matF', 'matC') )

# get matrices
get_mats <- function(rows_ids, mat_c_ids, raw_mat){

  list( matA = raw_mat[rows_ids,mat_c_ids$matA],
        matU = raw_mat[rows_ids,mat_c_ids$matU],
        matF = raw_mat[rows_ids,mat_c_ids$matF],
        matC = raw_mat[rows_ids,mat_c_ids$matC] )

}

# matrices in a list
mat_l <- lapply(mat_r_ids, get_mats, mat_c_ids, raw_mat)


# get metadata -----------------------------------------------------

# load compadre database for reference
compadre <- readRDS('data/Compadre_unrel_v4.rds')

# data we do not have in the excel file
cmp_nams <- names(compadre$metadata)
raw_nams <- c(raw_spp %>% names,
              raw_mat %>% names)
miss_nam <- setdiff(cmp_nams, raw_nams)

# first line of metadata for matrices and species
meta_m_id<- which((raw_mat %>% names) %in% cmp_nams)
meta_s_id<- which((raw_spp %>% names) %in% cmp_nams)

# metadata for matrices
meta_mat  <- lapply(first_row,
                    function(x,raw_mat,meta_cols) raw_mat[x,meta_cols],
                    raw_mat, meta_m_id) %>% bind_rows

# metadata for species
meta_spp  <- raw_spp[,meta_s_id] %>% uncount(nrow(meta_mat))

# metadata in a data frame
meta_df   <- bind_cols(meta_mat,meta_spp) %>%
                mutate( Lat = NA,
                        Lon = NA,
                        MatrixDimension = NA,
                        SurvivalIssue   = NA ) %>%
                select( cmp_nams )



# calculate lat/lon information


# function: convert lat/lon in decimal form
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

# get metadata on coordinates
meta_coord<- which((raw_mat %>% names) %in%
                    c("LatDeg", "LatMin", "LatSec",
                      "LonDeg", "LonMin", "LonSec") )

# calculate coordinates in decimal degrees
meta_coord<- lapply(first_row,
                    function(x,raw_mat,meta_cols) raw_mat[x,meta_cols],
                    raw_mat, meta_coord) %>%
                bind_rows %>%
                # combine Degrees/Minutes/Seconds
                mutate( lat_deg_min_sec = paste(LatDeg, LatMin, LatSec,
                                                sep=' '),
                        lon_deg_min_sec = paste(LonDeg, LonMin, LonSec,
                                                sep=' ') ) %>%
                # add a '-' to Longitude if the species is in US or extreme W europe
                # unfortunately, I don't know how to automate this
                mutate( lon_deg_min_sec = paste0('-',lon_deg_min_sec) ) %>%
                # convert coordinates to decimal
                mutate( lat = conv_plot_coord(lat_deg_min_sec,
                                              lon_deg_min_sec,
                                              'deg_min_sec')$lat,
                        lon = conv_plot_coord(lat_deg_min_sec,
                                              lon_deg_min_sec,
                                              'deg_min_sec')$lon )

# update with Lat/Lon information
meta_df <- meta_df %>%
              mutate( Lat = meta_coord$lat,
                      Lon = meta_coord$lon )


# get matrixClass -----------------------------------------------------

# get ids of matrix class info
m_class_id <- which(names(raw_mat) %in% names(compadre$matrixClass[[1]]) )

# get matrix classes
get_mat_class <- function(rows_ids, m_class_id, raw_mat){

  raw_mat[rows_ids,m_class_id]

}

# matrix class information in a list
mat_class_l <- lapply(mat_r_ids,
                      get_mat_class,
                      m_class_id, raw_mat)


# output a compadre-like object ---------------------------------------

# compadre-like object
out_mat <- list( metadata    = meta_df,
                 matrixClass = mat_class_l,
                 mat         = mat_l )

# store compadre-like object
saveRDS(out_mat,
        paste0('Data/',spp_nam,'.rds') )
