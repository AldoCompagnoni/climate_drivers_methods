library(rstan)
library(dplyr)
library(tidyr)

# full models
mod_full  <- c('_null',          '_yr',
               '_gaus',          '_expp',
               '_gev_nest',      '_expp_nest', '_dirichlet_nest',
               '_movbeta_h',     '_movbeta',
               '_movbeta_h_nest','_movbeta_n') %>% paste0('.stan')

# crossval models
mod_cross <- c('_null',          '_yr',
               '_gaus',          '_expp',
               '_gev_nest',      '_expp_nest', '_dirichlet_nest',
               '_movbeta_h',     '_movbeta',
               '_movbeta_h_nest','_movbeta_n') %>% paste0('_crossval.stan')

# all models
mod_all   <- expand.grid( fam  = c('normal', 'beta', 'gamma'),
                          root = c(mod_full, mod_cross) ) %>% 
                mutate( model = paste0(fam, root) )
                          
# function: comile models and save RDS file
mod_compile_store <- function(mod_x){
  
  mod_c     <- stanc( file = paste0(mod_x) )
  mod       <- stan_model( stanc_ret = mod_c, save_dso = TRUE)
  out_file  <- gsub('.stan','',mod_x)
  saveRDS( mod, file = paste0(out_file,".RDS") )
  
}

# comile and save models
lapply(mod_all$model, mod_compile_store)
