rm(list=ls())
source("R/format_data.R")
library(tidyverse)
library(mgcv)
library(testthat)
library(rstan)
library(loo)

# set up design data frame
design_df <- expand.grid( response = c('fec','grow','surv'),
                          spp      = 1:39 )


# read data ----------------------------------------------------------------------------------------
lam       <- read.csv("data/all_demog_updt.csv", stringsAsFactors = F)
spp       <- lam$SpeciesAuthor %>% unique

# format data --------------------------------------------------------------------------------------

# calculate mean and variance for each dataset
mean_var <- function( ii ){
  
  print(ii)
  
  spp_i         <- design_df$spp[ ii ]
  resp          <- design_df$response[ ii ]
  spp_name      <- spp[spp_i]
  
  # lambda data
  spp_resp      <- format_species(spp_name, lam, resp) %>% bind_rows
    
  # mean-variance output
  out <- data.frame( species  = spp_name,
                     response = resp,
                     mean     = spp_resp[,3] %>% mean(na.rm=T),
                     variance = spp_resp[,3] %>% var(na.rm=T) )
  
  if( resp == 'fec' ){
    out <- out %>% 
            mutate( mean     = log(mean),
                    variance = log(variance) )
  }
  
  out
  
}

# mean variance relationships
meanvar_df <- lapply( 1:nrow(design_df), mean_var) %>% bind_rows

# plot mean variance relationship
ggplot(meanvar_df) +
  geom_point( aes(mean, variance) ) + 
  facet_wrap( ~response, scale = 'free' ) + 
  theme_minimal() + 
  ggsave( 'results/prior_pc/mean_var_relationship.tiff',
          width = 6.3, height = 3, compression = 'lzw')
  

# Estimate mean-variance relationships --------------------------------------------------

# Fecundity
fec_mmeanvar_mod  <- meanvar_df %>% 
                       subset( response == 'fec' ) %>% 
                       subset( !(mean == -Inf) ) %>% 
                       lm( variance ~ mean, data = .)

write.csv( summary(fec_mmeanvar_mod)$coefficients, 
           'results/prior_pc/fec_meanvar_coef.csv',
           row.names = F)

meanvar_df %>% 
  subset( response == 'fec' ) %>% 
  subset( !(mean == -Inf) ) %>% 
  ggplot() + 
  geom_point( aes(mean, variance) ) + 
  theme_minimal() + 
  title( 'Fecundity: Mean-Variance Relationship') + 
  labs( x = 'log(mean)',
        y = 'log(variance)' ) +
  theme( axis.title = element_text( size = 20),
         axis.text  = element_text( size = 20) ) +
  ggsave( 'results/prior_pc/fec_mean_var_relationship.tiff',
          width = 6.3, height = 6.3, compression = 'lzw')


# Growth
grow_mmeanvar_df  <- meanvar_df %>% 
                       subset( response == 'grow' ) %>% 
                       subset( !(mean == -Inf) )
  
# Fit nonlinear function (From Jongejan's 2010 Appendix)
fit_nl <- function(params, in_df ){
  
  alpha    <- params[1]
  theta    <- params[2]
  sigma    <- params[3]
  yhat     <- alpha * (in_df$mean - (in_df$mean^2))^theta
  NLL      <- -sum( dnorm(in_df$variance,
                          mean = yhat,
                          sd   = sigma,log=T) )
   
  
  null_par  <<- c("lam.", "b.", "size") 
  return(NLL)
  
}

mvgrow_mod <- optim( par=c(0.1,0.1,0.1), fn=fit_nl, gr=NULL, 
                     grow_mmeanvar_df, 
                     control=list(maxit=10000000) )

x_seq      <- seq(0.01,0.99,length.out=100)
y_pred     <- mvgrow_mod$par[1]*( x_seq - (x_seq^2))^mvgrow_mod$par[2]
resid      <- grow_mmeanvar_df$variance - (mvgrow_mod$par[1]*( grow_mmeanvar_df$mean - (grow_mmeanvar_df$mean^2))^mvgrow_mod$par[2])

tiff( 'results/prior_pc/grow_mean_var_relationship.tiff',
      width = 6.3, height = 6.3, unit = 'in', 
      res = 500, compression = 'lzw')
plot( variance ~ mean, data = grow_mmeanvar_df)
lines( x_seq, y_pred, lwd = 3 )
text( 0.25, 0.11, expression('y='*0.12*'(x - x'^2*')'^0.91) )
dev.off()


# residuals
plot( grow_mmeanvar_df$mean, resid )


data.frame( alpha = mvgrow_mod$par[1],
            theta = mvgrow_mod$par[2],
            sigma = mvgrow_mod$par[3] ) %>% 
  write.csv( 'results/prior_pc/grow_meanvar_coef.csv',
             row.names = F)
