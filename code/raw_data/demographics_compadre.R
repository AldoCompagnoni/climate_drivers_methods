# Produce demographic metrics for 'new' species (Mar 12.3.2019)
# Pediocactus_bradyi, Astragalus_scaphoides_6, and Cryptantha_flava_2 
library(dplyr)
library(testthat)
library(measurements)
library(popbio)

# 1. get demographics for:
# Pediocactus_bradyi directly from compadre
# Astragalus_scaphoides_6, and 
# Cryptantha_flava_2 from compadre/like objects
# 2. harmonize results with all_demog_6tr.csv


# read data -----------------------------------------------------------------

compadre  <- readRDS('data/compadre_unrel_v4.rds')
astr      <- readRDS('data/Astragalus_scaphoides_6.rds')
crypt     <- readRDS('data/Cryptantha_flava_2.rds')
all_demo  <- read.csv('data/all_demog_6tr.csv', stringsAsFactors = F)


# get demographics for all species  -----------------------------------------


# Function to calculate summary vital rates (not stage-dependent)
vitalRates_lambda <- function( mats_l ){
    
    matU    <- mats_l[[1]]$matU %>% as.matrix
    matF    <- mats_l[[1]]$matF %>% as.matrix
    matC    <- mats_l[[1]]$matC %>% as.matrix
  
    matA    <- matU+matF+matC
    matDim  <- dim(matA)[1]
    
    out     <- data.frame( "SurvivalSSD"     =NA,
                           "ProgressionSSD"  =NA,
                           "RetrogressionSSD"=NA,
                           "ReproductionSSD" =NA,
                           "ClonalitySSD"    =NA,
                           'lambda'          =NA)
    
    if( any(is.na(matU)) | 
        any(is.na(matF)) | 
        any(is.na(matC)) ){
      return(out)  
    } else{
      #Extracting SSD corrected vital rate values
      SSD     <- eigen.analysis(matA)$stable.stage
      f       <- colSums(matF)
      out     <- mutate(out, ReproductionSSD = mean(f*SSD) )
      c       <- colSums(matC)
      out     <- mutate(out, ClonalitySSD    = mean(c*SSD) )
      
      #Preparing survival-independent matrix to calculate individual growth rates
      uDistrib<- matrix(NA, ncol=matDim, nrow=matDim)
      u       <- colSums(matU)
      out     <- mutate(out, SurvivalSSD = mean(u*SSD) )
    
      #Making matrix for transitions conditional on survival
      for (j in which(u>0) ) uDistrib[,j] <- matU[,j] / u[j]
      UPrime  <- uDistrib
      UPrime[is.na(UPrime)] <- 0
      CPrime  <- colSums(matC)
      CPrime[is.na(CPrime)] <- 0
      UCPrime <- UPrime+CPrime
      
      #Extracting proxy to individual progressive growth rate
      UCPrimeGrowth=UCPrime
      UCPrimeGrowth[upper.tri(UCPrime, diag = T)]=NA
      UCPrimeGrowth[matDim,matDim]=UCPrime[matDim,matDim]  #Putting back the last element of stasis bc there is likely growth on the top of class
      out$ProgressionSSD=mean(colSums(UCPrimeGrowth,na.rm=T)*SSD)
      #Extracting proxy to individual retrogressive growth rate
      UCPrimeShrinkage=UCPrime
      UCPrimeShrinkage[lower.tri(UCPrime, diag = T)]=NA
      out$RetrogressionSSD=mean(colSums(UCPrimeShrinkage,na.rm=T)*SSD)
    
      out     <- mutate(out, lambda = Re(eigen(matA)$value[1]) )
      
      return(out)  
    }
    
}


# 1. get demographics -------------------------------------------

# loop through matrices
vr_by_mat_spp <- function(spp_compadre, comp_obj){
  
  # list for outputs
  vr_l    <- list()
  
  # identify ids in compadre-like object
  mat_ids <- which(comp_obj$metadata$SpeciesAuthor == spp_compadre &
                   comp_obj$metadata$MatrixComposite == 'Individual' )

  # loop through matrices
  for(ii in 1:length(mat_ids) ){
    vr_l[[ii]] <- comp_obj$mat[mat_ids[ii]] %>% vitalRates_lambda
  }

  # combine metadata with vital rates
  cbind( comp_obj$metadata[mat_ids,],
          bind_rows( vr_l ) )
  
}

# get and format vital rates (add MatrixEndMonth by hand)
pedio_vr <- vr_by_mat_spp('Pediocactus_bradyi', compadre) %>% 
              mutate( MatrixEndMonth = 3 )
astr_vr  <- vr_by_mat_spp('Astragalus_scaphoides_6', astr) %>% 
              mutate( MatrixEndMonth = 6 )
crypt_vr <- vr_by_mat_spp('Cryptantha_flava_2', crypt) %>% 
              mutate( MatrixEndMonth = 5 )


# harmonize results with all_demog_6tr.csv ------------------------------------

# convert all columns to character (to allow bind_rows)
all_t0_char <- function(input_df){
  data.frame(lapply(input_df, as.character), 
             stringsAsFactors=FALSE)
}

# put it all together!
new_vr   <- list( pedio_vr,
                  astr_vr,
                  crypt_vr) %>% 
              # convert everything to character
              lapply(all_t0_char) %>% 
              bind_rows %>% 
              # change names to pair with all_demo data frame
              mutate( surv = SurvivalSSD,
                      grow = ProgressionSSD,
                      shri = RetrogressionSSD,
                      fec  = ReproductionSSD,
                      clo  = ClonalitySSD ) %>% 
              # remove NAs in Pediocactus
              subset( !is.na(lambda) )

# Which columns to convert to numeric (new_vr is all character)?
common_v   <- intersect(names(new_vr), names(all_demo))
id_num     <- which(sapply(dplyr::select(all_demo,common_v),class) == 'numeric' | 
                    sapply(dplyr::select(all_demo,common_v),class) == 'integer' )
conv_v     <- names(id_num)

# columns in new_vr converted to numeric
new_vr_num <- data.frame( lapply(dplyr::select(new_vr,conv_v), as.numeric) )

# substitute converted columns in original data frame
sub_id     <- which(names(new_vr) %in% names(new_vr_num))
new_vr[,
   sub_id] <- new_vr_num


# new data harmonized in all_demog file
bind_rows(all_demo, new_vr) %>% 
  write.csv('data/all_demog_updt.csv',
            row.names=F)

