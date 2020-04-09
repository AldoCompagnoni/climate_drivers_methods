library(rstan)

mod_dir   <- 'R/simulation/regularization/gamma_horse_brms.stan'

mod_dir   <- 'gamma_horse_brms.stan'
mod_c     <- stanc( file = mod_dir )
mod       <- stan_model( stanc_ret = mod_c, save_dso = TRUE)
out_file  <- gsub('.stan','',mod_dir)
saveRDS( mod, file = paste0(out_file,".RDS") )
