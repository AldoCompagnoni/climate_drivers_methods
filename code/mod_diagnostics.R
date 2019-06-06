# produce summary graphs
# 1. LOOIC results from 
# 2. LPPD results from supercomputer
source("code/format_data.R")
library(testthat)
library(dplyr)
library(tidyr)
library(viridis)
require(gridExtra)
library(egg)
options(stringsAsFactors = F )

# 1. LOOIC results from  -------------------------------------------------------
out_dir  <- 'C:/CODE/climate_drivers_methods/'
clim_var <- 'precip'
resp_l   <- c('log_lambda','surv','grow','fec')

# demograhpic information
lam      <- read.csv("data/all_demog_updt.csv", stringsAsFactors = F)
                     
# species list
spp      <- lam %>% .$SpeciesAuthor %>% unique 

# response
response <- resp_l[1]

# create information on total replication
create_rep_n <- function(x, lam, response){
  
  data.frame( species = x,
              rep_n   = format_species(x, lam, response) %>% bind_rows %>% nrow,
              rep_yr  = format_species(x, lam, response) %>% sapply(nrow) %>% max,
              rep_p   = format_species(x, lam, response) %>% length,
              stringsAsFactors = F) 
  
}

# data frame on total replication
rep_n_df <- lapply(spp, create_rep_n, lam, response) %>% bind_rows


# model order no control
mod_ordr <- c('ctrl1','yr1','yr2','yr3',
              'gaus',    'expp',
              'gaus_n',  'expp_n', 'simpl_n',
              'movb_h',  'movb',
              'movb_h_n','movb_n')

# read in files
sum_f_l   <- grep('mod_summ',
                  list.files(paste0('results/mod_sel/',
                                    clim_var)),
                  value=T) 

# extract function
extract_element <- function(x,patt){
  regmatches(x, 
             gregexpr(patt, x) )
}
  
# read all summaries
f_paths   <- paste0('results/mod_sel/',clim_var,'/',sum_f_l)
spp_v     <- gsub( paste0('mod_summaries_',resp_l,'_') %>% 
                   paste(collapse='|') %>% 
                   paste0('|.csv'),
                   '',sum_f_l)
resp_v    <- gsub( 'mod_summaries_', '', sum_f_l) %>%
                # delete ANYTHING less than 101 character long that follows:
                # either fec or surv or grow or log_lambda
                gsub(paste0('(?<=',paste(resp_l,collapse='|'),').{1,100}'), 
                     '', ., perl=T)

read_format_summ <- function(x,spp_x,resp_x){
  read.csv(x) %>% 
    tibble::add_column( species  = spp_x, .before = 2) %>% 
    tibble::add_column( response = resp_x, .before = 2)
}

# all model summaries in one data frame
diag_df <- Map(read_format_summ, 
               f_paths, spp_v, resp_v) %>% 
               bind_rows %>% 
               left_join( rep_n_df ) %>% 
               # Add "flags" for convergence issues
               mutate( diverg   = n_diverg > 0,
                       rhat     = rhat_high > 0 ) %>% 
               # good to "order" predictors
               mutate( model    = factor(model, levels = mod_ordr),
                       response = factor(response, 
                                         levels = c('surv',
                                                    'grow',
                                                    'fec',
                                                    'log_lambda') ),
                       sppcode1 = substr(species,1,2),
                       sppcode2 = substr(gsub('.{1,40}_', '', species, perl=T), 1, 2)
                       ) %>% 
               # make species code
               mutate( sppcode  = paste0(sppcode1,sppcode2) %>% toupper )


# focus on divergent transitions only
mod_div_df <- diag_df %>% 
                dplyr::select(model,sppcode,response,n_diverg) 

# fitted models 
fitted_df   <- diag_df %>% 
                  subset( n_diverg  == 0 & 
                          rhat_high == 0 & 
                          n_eff_low == 0 & 
                          mcse_high == 0 )

# divergent transitions
diverg_df  <- diag_df %>% 
                  subset( !(n_diverg  %in% 0) )

# Rhat < 1.1 
rhat_df    <- diag_df %>% 
                  subset( !(rhat_high  %in% 0) )


# divergent transition ----------------------------------------
p1 <- ggplot(mod_div_df, aes(x=model, y=n_diverg) ) +
        geom_jitter( aes(color = response) ) +
        scale_color_viridis( discrete = T ) + 
        theme( axis.text     = element_text( angle = 90),
               legend.text   = element_text(size=7),
               legend.title  = element_text(size=7),
               legend.margin = margin(0,0,0,0,unit='pt'))

p2 <- ggplot(mod_div_df, aes(x=sppcode, y=n_diverg) ) +
        geom_jitter( aes(color = response) ) +
        theme( axis.text = element_text( angle = 90),
               legend.position = 'none',
               legend.margin   = 1) +
        scale_color_viridis( discrete = T)

out_p <- grid.arrange(p1, p2, nrow=2, ncol=1) 
ggsave(file = 'results/converg/div_mod_spp.tiff',
       plot = out_p, width = 6.3, height = 6.3,
       compression="lzw")


# Nhat ---------------------------------------------------
p1 <- ggplot(diag_df, aes(x=model, y=rhat_high) ) +
        geom_jitter( aes(color = response) ) +
        scale_color_viridis( discrete = T ) + 
        theme( axis.text     = element_text( angle = 90),
               legend.text   = element_text(size=7),
               legend.title  = element_text(size=7),
               legend.margin = margin(0,0,0,0,unit='pt'))

p2 <- ggplot(diag_df, aes(x=sppcode, y=rhat_high) ) +
        geom_jitter( aes(color = response) ) +
        theme( axis.text = element_text( angle = 90),
               legend.position = 'none',
               legend.margin   = 1) +
        scale_color_viridis( discrete = T)

out_p <- grid.arrange(p1, p2, nrow=2, ncol=1) 
ggsave(file = 'results/converg/rhat_mod_spp.tiff',
       plot = out_p, width = 6.3, height = 6.3,
       compression="lzw")


# diagnostics versus rep_yr -----------------------------------
p1 <- diag_df %>% 
        subset( !(model %in% c('ctrl1','yr1','yr2','yr3')) ) %>% 
        ggplot( aes(x=rep_yr, y=n_diverg) ) +
        geom_jitter( alpha = 0.8,
          aes(color = model) ) +
        scale_color_viridis( discrete = T ) +
        ylab( expression('Number of divergent transitions') ) + 
        xlab( 'Temporal replication of dataset (years)' )

p2 <- diag_df %>% 
        subset( !(model %in% c('ctrl1','yr1','yr2','yr3')) ) %>% 
        ggplot( aes(x=rep_yr, y=rhat_high) ) +
        geom_jitter( alpha = 0.8,
          aes(color = model) ) +
        scale_color_viridis( discrete = T ) +
        ylab( expression('Number of '*hat(R)*' > 1.01') ) + 
        xlab( 'Temporal replication of dataset (years)' )

out_p <- grid.arrange(p1, p2, nrow=2, ncol=1) 
ggsave(file = 'results/converg/diag_vs_yr_rep.tiff',
       plot = out_p, width = 6.3, height = 6.3,
       compression="lzw")


# plots ------------------------------------------------------------
prop_div <- diag_df %>% 
              group_by( rep_yr) %>% 
              summarise( tot  = sum(diverg),
                         rep  = n() ) %>% 
              ungroup %>% 
              mutate( p_div_yr = tot / rep )

prop_hat <- diag_df %>% 
              group_by( rep_yr) %>% 
              summarise( tot  = sum(rhat),
                         rep  = n() ) %>% 
              ungroup %>% 
              mutate( p_rhat = tot / rep )

prop_mod_div <- diag_df %>% 
                    group_by( model ) %>% 
                    summarise( tot  = sum(diverg),
                               rep  = n() ) %>% 
                    ungroup %>% 
                    mutate( p_div = tot / rep )

prop_mod_rhat <- diag_df %>% 
                    group_by( model ) %>% 
                    summarise( tot  = sum(rhat),
                               rep  = n() ) %>% 
                    ungroup %>% 
                    mutate( p_rhat = tot / rep )



# proportions with Rhat issue
p1 <- ggplot(prop_mod_rhat, aes(model, p_rhat)) +
      geom_jitter(aes(size = '0.7') ) +
      # geom_jitter(aes(color = response) ) +
      scale_color_viridis(discrete = T) + 
      ylim(0 ,1 ) + 
      ylab( 'Prop. Rhat issue' )+
      xlab( 'Model' ) +
      theme( axis.text.x  = element_text(angle = 90, hjust = 1), 
             axis.title      = element_text( size = 15),
             legend.position="none" )

p2 <- ggplot(prop_mod_div, aes(model, p_div) ) +
      geom_jitter(aes(size = '0.7')) +
      # geom_jitter(aes(color = response)) +
      scale_color_viridis(discrete = T) + 
      ylim(0 ,1 ) +
      ylab( 'Prop. divergence issue' )+
      xlab( 'Model' ) +
      theme( axis.text.x     = element_text(angle = 90, hjust = 1),
             axis.title      = element_text( size = 15),
             legend.position = "none" )

p3 <- ggplot(prop_hat, aes(rep_yr, p_rhat) ) +
      geom_jitter(aes(size = '0.7')) +
      # geom_jitter(aes(color = response)) +
      scale_color_viridis(discrete = T) + 
      ylim(0 ,1 ) +
      ylab( 'Prop. Rhat issue' ) +
      xlab( 'Years study' ) +
      theme( axis.text.x     = element_text(angle = 90, hjust = 1),
             axis.title      = element_text( size = 15),
             legend.position = "none" )

p4 <- ggplot(prop_div, aes(rep_yr, p_div_yr) ) +
      geom_jitter(aes(size = '0.7')) +
      # geom_jitter(aes(color = response)) +
      scale_color_viridis(discrete = T) + 
      ylim(0 ,1 ) +
      ylab( 'Prop. Divergence issue' ) +
      xlab( 'Years study' ) +
      theme( axis.text.x     = element_text(angle = 90, hjust = 1),
             axis.title      = element_text( size = 15) )


# proportion of models with convergence issues
out_p <- grid.arrange(p1, p2, p3, p4, nrow=2, ncol=2) 
ggsave(file = 'results/converg/prop_conv_issues.tiff',
       plot = out_p, width = 6.3, height = 6.3,
       compression="lzw")



# absolute number of issues
p1 <- ggplot(diag_df, aes(model, n_diverg)) +
      geom_point(aes(size = '0.7') ) +
      ylab( 'N. of divergent transitions' ) + 
      xlab( 'Model' ) +
      theme( axis.text.x     = element_text(angle = 90, hjust = 1 ),
             axis.title      = element_text( size = 15 ),
             legend.position = "none" )

p2 <- ggplot(diag_df, aes(model, rhat_high) ) +
      geom_point(aes(size = '0.7') ) +
      ylab( 'N. of Rhat above 1.1' ) + 
      xlab( 'Model' ) +
      theme( axis.text.x     = element_text(angle = 90, hjust = 1),
             axis.title      = element_text( size = 15),
             legend.position = "none" )

# absolute number of convergence issues
out_p <- grid.arrange(p1,p2,nrow=2)
ggsave(file = 'results/converg/conv_issues_nums.tiff',
       plot = out_p, width = 6.3, height = 6.3,
       compression="lzw")



# absolute number of issues BY SPECIES
p1 <- ggplot(diag_df, aes(species, n_diverg)) +
      geom_point(aes(size = '0.7') ) +
      ylab( 'N. of divergent transitions' ) + 
      xlab( 'Species' ) +
      theme( axis.text.x     = element_blank(),
             axis.title.x    = element_blank(),
             axis.title      = element_text( size = 15 ),
             legend.position = "none" )

p2 <- mutate(diag_df, 
             species = substr(species,1,5) ) %>%  
      ggplot(aes(species, rhat_high) ) +
      geom_point(aes(size = '0.7') ) +
      ylab( 'N. of Rhat above 1.1' ) + 
      xlab( 'Species' ) +
      theme( axis.text.x     = element_text(angle = 90, hjust = 1),
             axis.title      = element_text( size = 15),
             legend.position = "none" )

# 
out_p <- grid.arrange(p1,p2,nrow=2)
ggsave(file = 'results/converg/conv_issues_spp.tiff',
       plot = out_p, width = 6.3, height = 6.3,
       compression="lzw")



# other gibberish 

# summarize convergence problems
mods <- c('ctrl1','ctrl2','yr1','yr2','yr3','yr_bet','yr_wgt',
          'expp','gaus','gev','gev_n','expp_n','simpl_n')

# 
diverg_df$model    %>% table %>% .[mods]
diverg_df$response %>% table
diverg?loo_df$clim_var %>% table
diverg_df$rep_n    %>% table
diverg_df$rep_yr   %>% table

rhat_df$model    %>% table %>% .[mods]
rhat_df$response %>% table 
rhat_df$clim_var %>% table 
rhat_df$rep_n    %>% table
rhat_df$rep_yr   %>% table









# 2. LPPD results from supercomputer  ------------------------------------------

today_date    <- gsub("-","_",Sys.time() %>% as.Date %>% as.character)
m_back        <- 36
interval      <- NULL
pre_chelsa    <- NULL # '_pre_chelsa'


# Summarize moving windows results by climate variable -------------------------------
mod_perform <- function(ii){
  
  clim_var <- input_df$clim_var[ii]
  response <- input_df$response[ii]
  interval <- input_df$interval[ii]
  resp_clim<- paste0("_",response,"_",clim_var)
  
  # read lambda/clim data 
  lam     <- read.csv("C:/cloud/Dropbox/sAPROPOS/all_demog_6tr.csv", 
                      stringsAsFactors = F) #%>%
                #subset( SpeciesAuthor != "Purshia_subintegra" )

  # summary info 
  
  # result folder name
  res_folder<- paste0("C:/cloud/Dropbox/sAPROPOS/results/moving_windows/",
                      response,
                      "/",
                      clim_var,pre_chelsa,interval) 
  # summary files names
  sum_files <- list.files(res_folder)[grep("mod_summaries_", list.files(res_folder) )] %>% 
                  stringr::str_subset(resp_clim)
  # read files
  mod_summ  <- lapply(sum_files, function(x) read.csv(paste0(res_folder,"/",x)) ) %>%
                  setNames( gsub("mod_summaries_", "", sum_files ) ) %>%
                  setNames( gsub(paste0(resp_clim,".csv"), "", names(.) ) )
  # all model selection summaries
  all_sums  <- Map(function(x,y) tibble::add_column(x, species = y, .before = 1), 
                   mod_summ, names(mod_summ) ) %>% 
                  lapply(function(x) 
                    dplyr::select(x,species, model, 
                                  n_diverg, rhat_high, 
                                  n_eff_low, mcse_high) ) %>% 
                  bind_rows
          
  # species
  spp           <- names(mod_summ)
 
  # create information on total replication
  create_rep_n <- function(x, lam, response){
    
    data.frame( species = x,
                rep_n   = format_species(x, lam, response) %>% bind_rows %>% nrow,
                rep_yr  = format_species(x, lam, response) %>% sapply(nrow) %>% max,
                rep_p   = format_species(x, lam, response) %>% length,
                stringsAsFactors = F) 
    
  }
  rep_n_df <- lapply(spp, create_rep_n, lam, response) %>% bind_rows
  
  # spit it out
  all_sums %>% 
    left_join( rep_n_df ) %>% 
    mutate( clim_var  = clim_var,
            response  = response )

}

# all models results
input_df    <- expand.grid( clim_var = c("precip","airt"),
                            response = c("surv","grow","fec","log_lambda"),
                            interval = "",
                            stringsAsFactors = F)

# diagnostics df
diag_df     <- lapply(1:nrow(input_df), mod_perform) %>% 
                  bind_rows %>% 
                  # Add "flags" for convergence issues
                  mutate( diverg  = n_diverg > 0,
                          rhat    = rhat_high > 0 ) %>% 
                  mutate( species = gsub('_var._gnaphalifolium_2','',species) )

# model convergence overview --------------------------------------------

# fitted models 
fitted_df   <- diag_df %>% 
                  subset( n_diverg  == 0 & 
                          rhat_high == 0 & 
                          n_eff_low == 0 & 
                          mcse_high == 0 )
           
# divergent transitions
diverg_df  <- diag_df %>% 
                  subset( !(n_diverg  %in% 0) )


# Rhat < 1.1 
rhat_df    <- diag_df %>% 
                  subset( !(rhat_high  %in% 0) )


# test predictors of convergence issues as a binomial process

# divergent transitions
div_m  <- glm (diverg ~ model,    data = diag_df, family=binomial)
div_r  <- glm (diverg ~ response, data = diag_df, family=binomial)
div_c  <- glm (diverg ~ clim_var, data = diag_df, family=binomial)
div_y  <- glm (diverg ~ rep_yr,   data = diag_df, family=binomial)
div_n  <- glm (diverg ~ rep_p,    data = diag_df, family=binomial)
div_ym <- glm (diverg ~ model + rep_p, data = diag_df, family=binomial)

# model type is the biggest predictor, yr_replication runner up
AIC(div_m, div_r, div_c, div_y, div_n, div_ym)


# divergent transitions
hat_m  <- glm (rhat ~ model,    data = diag_df, family=binomial)
hat_r  <- glm (rhat ~ response, data = diag_df, family=binomial)
hat_c  <- glm (rhat ~ clim_var, data = diag_df, family=binomial)
hat_y  <- glm (rhat ~ rep_yr,   data = diag_df, family=binomial)
hat_n  <- glm (rhat ~ rep_p,    data = diag_df, family=binomial)
hat_ym <- glm (rhat ~ model + rep_p, data = diag_df, family=binomial)

# Same story: model type is the biggest predictor, yr_replication runner up
AIC(hat_m, hat_r, hat_c, hat_y, hat_n, hat_ym)



# test predictors of convergence issues as a POISSON process

# divergent transitions
div_p_m  <- glm (n_diverg ~ model,    data = diag_df, family=poisson)
div_p_r  <- glm (n_diverg ~ response, data = diag_df, family=poisson)
div_p_c  <- glm (n_diverg ~ clim_var, data = diag_df, family=poisson)
div_p_y  <- glm (n_diverg ~ rep_yr,   data = diag_df, family=poisson)
div_p_n  <- glm (n_diverg ~ rep_p,    data = diag_df, family=poisson)
div_p_ym <- glm (n_diverg ~ model + rep_p, data = diag_df, family=poisson)

# model type is the biggest predictor, yr_replication runner up
AIC(div_p_m, div_p_r, div_p_c, div_p_y, div_p_n, div_p_ym)


# divergent transitions
hat_p_m  <- glm (rhat_high ~ model,    data = diag_df, family=poisson)
hat_p_r  <- glm (rhat_high ~ response, data = diag_df, family=poisson)
hat_p_c  <- glm (rhat_high ~ clim_var, data = diag_df, family=poisson)
hat_p_y  <- glm (rhat_high ~ rep_yr,   data = diag_df, family=poisson)
hat_p_n  <- glm (rhat_high ~ rep_p,    data = diag_df, family=poisson)
hat_p_ym <- glm (rhat_high ~ model + rep_p, data = diag_df, family=poisson)

# Same story: model type is the biggest predictor, yr_replication runner up
AIC(hat_p_m, hat_p_r, hat_p_c, hat_p_y, hat_p_n, hat_p_ym)


# plots ------------------------------------------------------------

# would be nice to give model name along x-axis
coef(hat_m)[-1] %>% sort %>% plot(xlab='Model',ylab='Coefficient Rhat')
coef(div_m)[-1] %>% sort %>% plot(xlab='Model',ylab='Coefficient Divergent')


prop_div <- diag_df %>% 
              group_by( rep_yr) %>% 
              summarise( tot  = sum(diverg),
                         rep  = n() ) %>% 
              ungroup %>% 
              mutate( p_div_yr = tot / rep )

prop_hat <- diag_df %>% 
              group_by( rep_yr) %>% 
              summarise( tot  = sum(rhat),
                         rep  = n() ) %>% 
              ungroup %>% 
              mutate( p_rhat = tot / rep )

prop_mod_div <- diag_df %>% 
                    group_by( model ) %>% 
                    summarise( tot  = sum(diverg),
                               rep  = n() ) %>% 
                    ungroup %>% 
                    mutate( p_div = tot / rep )

prop_mod_rhat <- diag_df %>% 
                    group_by( model ) %>% 
                    summarise( tot  = sum(rhat),
                               rep  = n() ) %>% 
                    ungroup %>% 
                    mutate( p_rhat = tot / rep )


par( mfrow=c(2,2), mar=c(3,3,0.1,0.1), mgp=c(1.5,0.8,0))

plot(p_rhat   ~ rep_yr,   data=prop_hat, 
     ylim=c(0,1), pch=16)
plot(p_div_yr ~ rep_yr, data=prop_div, 
     ylim=c(0,1), pch=16 )


# proportions with Rhat issue
p1 <- ggplot(prop_mod_rhat, aes(model, p_rhat)) +
      geom_point(aes(size = '0.7') ) +
      ylim(0 ,1 ) +
      ylab( 'Prop. Rhat issue' )+
      xlab( 'Model' ) +
      theme( axis.text.x  = element_text(angle = 90, hjust = 1), 
             axis.title      = element_text( size = 15),
             legend.position="none" )

p2 <- ggplot(prop_mod_div, aes(model, p_div) ) +
      geom_point(aes(size = '0.7') ) +
      ylim(0 ,1 ) +
      ylab( 'Prop. divergence issue' )+
      xlab( 'Model' ) +
      theme( axis.text.x     = element_text(angle = 90, hjust = 1),
             axis.title      = element_text( size = 15),
             legend.position = "none" )

p3 <- ggplot(prop_hat, aes(rep_yr, p_rhat) ) +
      geom_point(aes(size = '0.7') ) +
      ylim(0 ,1 ) +
      ylab( 'Prop. Rhat issue' ) +
      xlab( 'Years study' ) +
      theme( axis.text.x     = element_text(angle = 90, hjust = 1),
             axis.title      = element_text( size = 15),
             legend.position = "none" )

p4 <- ggplot(prop_div, aes(rep_yr, p_div_yr) ) +
      geom_point(aes(size = '0.7') ) +
      ylim(0 ,1 ) +
      ylab( 'Prop. Divergence issue' ) +
      xlab( 'Years study' ) +
      theme( axis.text.x     = element_text(angle = 90, hjust = 1),
             axis.title      = element_text( size = 15),
             legend.position = "none" )
# plot(p_rhat   ~ rep_yr,   data=prop_hat, 
#      ylim=c(0,1), pch=16)
# plot(p_div_yr ~ rep_yr, data=prop_div, 
#      ylim=c(0,1), pch=16 )

tiff( 'results/conv_issues.tiff', 
      unit="in", width=6.3, height=6.3, res=600,compression="lzw" )
grid.arrange(p1, p2, 
             p3, p4, nrow=2, ncol=2)
dev.off()


# absolute number of issues
p1 <- ggplot(diag_df, aes(model, n_diverg)) +
      geom_point(aes(size = '0.7') ) +
      ylab( 'N. of divergent transitions' ) + 
      xlab( 'Model' ) +
      theme( axis.text.x     = element_text(angle = 90, hjust = 1 ),
             axis.title      = element_text( size = 15 ),
             legend.position = "none" )

p2 <- ggplot(diag_df, aes(model, rhat_high) ) +
      geom_point(aes(size = '0.7') ) +
      ylab( 'N. of Rhat above 1.1' ) + 
      xlab( 'Model' ) +
      theme( axis.text.x     = element_text(angle = 90, hjust = 1),
             axis.title      = element_text( size = 15),
             legend.position = "none" )

tiff( 'results/conv_issues_pois.tiff', 
      unit="in", width=3.15, height=6.3, res=600,compression="lzw" )
grid.arrange(p1,p2,nrow=2)
dev.off()



# absolute number of issues BY SPECIES
p1 <- ggplot(diag_df, aes(species, n_diverg)) +
      geom_point(aes(size = '0.7') ) +
      ylab( 'N. of divergent transitions' ) + 
      xlab( 'Species' ) +
      theme( axis.text.x     = element_text(angle = 90, hjust = 1 ),
             axis.title      = element_text( size = 15 ),
             legend.position = "none" )

p2 <- ggplot(diag_df, aes(species, rhat_high) ) +
      geom_point(aes(size = '0.7') ) +
      ylab( 'N. of Rhat above 1.1' ) + 
      xlab( 'Species' ) +
      theme( axis.text.x     = element_text(angle = 90, hjust = 1),
             axis.title      = element_text( size = 15),
             legend.position = "none" )

tiff( 'results/conv_issues_pois.tiff', 
      unit="in", width=3.15, height=6.3, res=600,compression="lzw" )
grid.arrange(p1,p2,nrow=2)
dev.off()



# other gibberish -----------------------------------------------

# summarize convergence problems
mods <- c('ctrl1','ctrl2','yr1','yr2','yr3','yr_bet','yr_wgt',
          'expp','gaus','gev','gev_n','expp_n','simpl_n')

# 
diverg_df$model    %>% table %>% .[mods]
diverg_df$response %>% table
diverg_df$clim_var %>% table
diverg_df$rep_n    %>% table
diverg_df$rep_yr   %>% table

rhat_df$model    %>% table %>% .[mods]
rhat_df$response %>% table 
rhat_df$clim_var %>% table 
rhat_df$rep_n    %>% table
rhat_df$rep_yr   %>% table
