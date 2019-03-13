# Script to format Cryptantha's .csv file
# cannot use "format_raw_matrix" because these are not matrices in "compadre format"
library(dplyr)
library(tidyr)
library(measurements)

# 1. get rows of separate matrices
# 2. get matrices matA, matU, matF, matC
# 3. get metadata
# 4. Don't do it, missing data: 'get matrixClass'
# 5. output a compadre-like object

# read in data
compadre  <- readRDS('data/Compadre_unrel_v4.rds')
raw_mat   <- read.csv('data/raw_mats/Cryptantha_flava.csv',
                      stringsAsFactors = F)[,-1]
spp_nam   <- 'Cryptantha_flava_2'

# get rows of separate matrices ---------------------------------

# get indices of first row
first_row <- setdiff( c(1:nrow(raw_mat)),
                        which(is.na(raw_mat$MatrixStartYear)) )

# matrix dimension
if( is.na(first_row[2]) ){
  # if we only have 1 matrix
  mat_dim   <- nrow(raw_mat)
}else{
  mat_dim   <- (first_row[2]-1)
}

# indices raws associated w/ separate matrices
mat_r_ids  <- lapply(first_row, function(x, add_rows) x:(x+add_rows),
                     mat_dim-1)


# get matrices matA, matU, matF, matC ------------------------------------

# where do matrices start in the spreadsheet?
first_col <- which(names(raw_mat) == 'MatrixTreatment') + 1

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


# get metadata --------------------------------------------------

# data we do not have in the excel file
cmp_nams <- names(compadre$metadata)
raw_nams <- raw_mat %>% names
miss_nam <- setdiff(cmp_nams, raw_nams)

# first line of metadata for matrices and species
meta_m_id<- which((raw_mat %>% names) %in% cmp_nams)

# metadata for matrices
meta_mat  <- lapply(first_row,
                    function(x,raw_mat,meta_cols) raw_mat[x,meta_cols],
                    raw_mat, meta_m_id) %>% 
                bind_rows 

# compadre already associated with Lucas et al. in COMPADRE
lucas_meta <- compadre$metadata %>% 
                subset( grepl('Cryptantha',SpeciesAuthor) ) %>% 
                dplyr::select(MatrixPopulation,Lat,Lon,Altitude,
                              Country, Continent, Ecoregion) %>% 
                unique %>% 
                uncount(nrow(meta_mat))


# metadata in a data frame
meta_df   <- bind_cols(meta_mat, lucas_meta) %>% 
                # these are all individual matrices
                mutate( MatrixComposite = 'Individual' )
                

# output a compadre-like object -----------------------------

# compadre-like object
out_mat <- list( metadata    = meta_df,
                 mat         = mat_l )

# store compadre-like object
saveRDS(out_mat,
        paste0('Data/',spp_nam,'.rds') )
