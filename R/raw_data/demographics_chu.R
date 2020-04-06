# extract demographic metrics from chu et al.'s matrices
library(dplyr)
library(purrr)
library(maps)
library(scales)
library(MASS)
library(popbio)
library(popdemo)
library(Matrix)
library(measurements)

# convert lat/lon ---------------------------------------------

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


# read in matrices ----------------------------------------------------------------

data_dir <- 'C:/cloud/Dropbox/sApropos/results/new_spp/adler/'

# load most current compadre database
# load('COMPADRE_v.X.X.X.4.RData')

# list folders
folders <- grep('^[[:upper:]]{4}$',
                list.files( data_dir ),
                value=T) %>% 
            c( 'BOER_NM', 'HECO_MT', 'POSE_MT')


# read matrices from chu et al. 
read_chu_mpms <- function(ii){
  
  do_spp  <- folders[ii]
  mats    <- read.csv( paste0(data_dir,do_spp,'/',do_spp,'_COMPADRE.csv'))
  mats
  
}

# download species-specific matrices.
mat_spp_l <- lapply(1:length(folders), read_chu_mpms)


# extract metavariables 
extract_metavar <- function(spp_mat_df){

  tot_rows  <- spp_mat_df %>% nrow /10 
  r_get     <- (0:(tot_rows-1)*10) + 1
  meta_var  <- spp_mat_df %>% 
                  dplyr::select(EnteredBy:MatrixClassOrganized) %>% 
                  names
  
  spp_mat_df[r_get, meta_var]
  
}

metadf <- lapply(mat_spp_l, extract_metavar) %>% 
            bind_rows %>% 
            # CORRECT degree Longitude sign
            mutate( LonDeg = -LonDeg )


# extract metavariables
extract_mats <- function(spp_mat_df){
  
  tot_rows  <- spp_mat_df %>% nrow /10 
  r_get     <- (0:(tot_rows-1)*10) + 1
  
  col_matA  <- paste0('A',1:10)
  col_matT  <- paste0('A',12:21)
  col_matF  <- paste0('A',23:32)

  # make a list of matA, matT, matF
  mat_ATF_list <- function(r_get, mat_df){
    
    out       <- list()
    out$matA  <- mat_df[r_get:(r_get+9), col_matA]
    out$matu  <- mat_df[r_get:(r_get+9), col_matT]
    out$matF  <- mat_df[r_get:(r_get+9), col_matF]
    out
    
  }
  
  lapply(r_get, mat_ATF_list, spp_mat_df)
  
}

# all Chu et al.'s matrices 
mat_l <- lapply(mat_spp_l, extract_mats) %>% flatten


# create demographics ----------------------------------------------------------------


# Function to calculate summary vital rates (not stage-dependent)
vitalRates_lambda <- function( list_of_mats ){
    
    matU <- list_of_mats$matu
    matF <- list_of_mats$matF
    matC <- matrix(0,10,10)
  
    matA=matU+matF+matC
    matDim=dim(matA)[1]
    
    out = data.frame("SurvivalSSD"=NA,
                     "ProgressionSSD"=NA,
                     "RetrogressionSSD"=NA,
                     "ReproductionSSD"=NA,
                     "ClonalitySSD"=NA)
    
    #Extracting SSD corrected vital rate values
    SSD=eigen.analysis(matA)$stable.stage
    f=colSums(matF)
    out$ReproductionSSD=mean(f*SSD)
    c=colSums(matC)
    out$ClonalitySSD=mean(c*SSD)
    
    #Preparing survival-independent matrix to calculate individual growth rates
    uDistrib=matrix(NA,ncol=matDim,nrow=matDim)
    u=colSums(matU)
    out$SurvivalSSD=mean(u*SSD)
  
    #Making matrix for transitions conditional on survival
    for (j in which(u>0)) uDistrib[,j]=matU[,j]/u[j]
    UPrime=uDistrib
    UPrime[is.na(UPrime)]=0
    CPrime=colSums(matC)
    CPrime[is.na(CPrime)]=0
    UCPrime=UPrime+CPrime
    
    #Extracting proxy to individual progressive growth rate
    UCPrimeGrowth=UCPrime
    UCPrimeGrowth[upper.tri(UCPrime, diag = T)]=NA
    UCPrimeGrowth[matDim,matDim]=UCPrime[matDim,matDim]  #Putting back the last element of stasis bc there is likely growth on the top of class
    out$ProgressionSSD=mean(colSums(UCPrimeGrowth,na.rm=T)*SSD)
    #Extracting proxy to individual retrogressive growth rate
    UCPrimeShrinkage=UCPrime
    UCPrimeShrinkage[lower.tri(UCPrime, diag = T)]=NA
    out$RetrogressionSSD=mean(colSums(UCPrimeShrinkage,na.rm=T)*SSD)
  
    out$lambda <- Re(eigen(list_of_mats$matA)$value[1])
    
    return(out)  
}

# get all vital rates 
vital_rates_df <- lapply(mat_l, vitalRates_lambda) %>% bind_rows

# all demographics
all_demog_chu  <- bind_cols(metadf, vital_rates_df) %>% 
                    # convert lat/lon to decimal
                    mutate( LatNS = conv_plot_coord(paste(LatDeg,LatMin,LatSec, sep=' '),
                                                    paste(LonDeg,LonMin,LonSec, sep=' '),
                                                    'deg_min_sec')$lat,
                            LonWE = conv_plot_coord(paste(LatDeg,LatMin,LatSec, sep=' '),
                                                    paste(LonDeg,LonMin,LonSec, sep=' '),
                                                    'deg_min_sec')$lon 
                            ) %>% 
                    # rename to lat/lon
                    rename( Lat = LatNS, 
                            Lon = LonWE )

# write down chu et al. ---------------------------------------------------
write.csv(all_demog_chu, 'results/demographics/chu_demog.csv',row.names=F)
