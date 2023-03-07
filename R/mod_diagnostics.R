# produce diagnostics for all models
# 1. Set up data
# 2. Proportion models with issues
# 3. Proportion of parameters with issues
#       (not doable, too complicated)
# 4. Rhat ~ n_eff
source("R/format_data.R")
library(testthat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggthemes)
library(viridis)
require(gridExtra)
library(egg)
library(stringr)
library(cowplot)
options(stringsAsFactors = F )

# quote a series of bare names
quote_bare <- function( ... ){
  substitute( alist(...) ) %>% 
    eval( ) %>% 
    sapply( deparse )
}

# set theme for plotting
theme_set( theme_minimal() )

# 1. Set up data  -------------------------------------------------------
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
mod_ordr <- quote_bare( ctrl1, 
                        yr1,    yr2,     yr3,
                        gaus1,  gaus2,   gaus3,
                        simpl1, simpl2,  simpl3,
                        ridge1, ridge2,  ridge3 )
                        # gaus,   simpl_n, ridge )

# data frame containing inputs
input_df  <- expand.grid( resp     = c('log_lambda','surv','grow','fec'),
                          clim_var = c('precip','airt'),
                          stringsAsFactors = F ) %>% 
                mutate( resp_clim = paste(resp, clim_var, sep='_') )

# where are the files?
# mod_dir     <- 'C:/Users/ac22qawo/sapropos_main/out_22.4.2020/'
# mod_dir     <- 'C:/Users/ac22qawo/sapropos_main/out_2021/'
# mod_dir     <- 'C:/Users/ac22qawo/sapropos_main/out_2021.2/'
mod_dir     <- 'C:/Users/ac22qawo/sapropos_main/out_2021.8.12/'

# file list
sum_f_l   <- grep('mod_summ',
                  list.files(mod_dir),
                  value=T) #%>% 
               # grep(resp_clim, ., value=T)

# path list
path_l    <- paste0(mod_dir, sum_f_l)

# extract species names
spp_names <- gsub("mod_summaries_", "", sum_f_l ) %>%
                gsub('array_vr-[0-9]{7}-[0-9]{1,3}-[0-9]{1,2}_', '', . ) %>% 
                gsub('.csv', "", . ) %>% 
                gsub( paste(input_df$resp_clim, collapse='|'),
                      '', .) %>% 
                gsub('_$','', .)

# species, response, and climate
spp_resp_clim <- gsub("mod_summaries_", "", sum_f_l ) %>%
                   gsub('array_vr-[0-9]{7}-[0-9]{1,3}-[0-9]{1,2}_', "", . ) %>% 
                   gsub('.csv', "", . )

# get response and climate variable (without spp info)
get_res_clim <- function(ii){
  gsub(spp_names[ii], '', spp_resp_clim[ii]) %>% 
    gsub('^_','',.)
}

# response and climate
resp_clim <- sapply(1:length(spp_resp_clim), get_res_clim)

# download files serially
create_sum_df <- function(path_l, 
                          spp_names,
                          resp_clim){
  
  read.csv(path_l, stringsAsFactors = F) %>% 
    tibble::add_column(species = spp_names,.before=1) %>% 
    tibble::add_column(resp_clim = resp_clim, .before=4)
  
}

# get all model summaries
mod_sum_df  <- Map(create_sum_df, path_l,
                   spp_names, resp_clim) %>% 
                  bind_rows %>% 
                  mutate( resp_clim = gsub('log_lambda', 'loglambda', resp_clim) ) %>% 
                  separate( resp_clim, c('response','clim'),sep='_') %>% 
                  mutate( response = gsub('loglambda', 'log_lambda', response) ) %>% 
                  subset( model %in% mod_ordr ) %>% 
                  subset( !(species %in% c('Astragalus_scaphoides_6_site_rep',
                                           'Astragalus_scaphoides_2')) )

# # temporarily store n. of divergent transitions
# select(mod_sum_df, 
#        species, model, response, 
#        clim, n_diverg) %>% 
#   subset( model == 'mb_h') %>% 
#   subset( n_diverg > 0) %>% 
#   write.csv('C:/Users/ac22qawo/Desktop/n_div.csv',
#             row.names=F)

# # read in files
# sum_f_l   <- grep('mod_summ',
#                   list.files(paste0('results/mod_sel/',
#                                     clim_var)),
#                   value=T) 
# 
# # extract function
# extract_element <- function(x,patt){
#   regmatches(x, 
#              gregexpr(patt, x) )
# }
#   
# # read all summaries
# f_paths   <- paste0('results/mod_sel/',clim_var,'/',sum_f_l)
# spp_v     <- gsub( paste0('mod_summaries_',resp_l,'_') %>% 
#                    paste(collapse='|') %>% 
#                    paste0('|.csv'),
#                    '',sum_f_l)
# resp_v    <- gsub( 'mod_summaries_', '', sum_f_l) %>%
#                 # delete ANYTHING less than 101 character long that follows:
#                 # either fec or surv or grow or log_lambda
#                 gsub(paste0('(?<=',paste(resp_l,collapse='|'),').{1,100}'), 
#                      '', ., perl=T)
# 
# # 
# read_format_summ <- function(x,spp_x,resp_x){
#   read.csv(x) %>% 
#     tibble::add_column( species  = spp_x,  .before = 2) %>% 
#     tibble::add_column( response = resp_x, .before = 2)
# }

# all model summaries in one data frame
# diag_df <- Map(read_format_summ, 
#                f_paths, spp_v, resp_v) %>% 
#                bind_rows %>% 

# Introduce model acronyms -----------------------------------------------------

# labels for models to plot on x-axis
mod_labs    <- quote_bare( ctrl1, 
                           yr1,    yr2,     yr3,
                           gaus1,  gaus2,   gaus3,
                           simpl1, simpl2,  simpl3,
                           ridge1, ridge2,  ridge3 )

new_labs    <- c( 'NM', 
                  'CSM 1', 'CSM 2', 'CSM 3',
                  'WMM 1', 'WMM 2', 'WMM 3',
                  'SAM 1', 'SAM 2', 'SAM 3',
                  'FHM 1', 'FHM 2', 'FHM 3' )

labs_df     <- data.frame( model     = mod_labs,
                           new_model = new_labs )



# Diagnostics data frame -----------------------------------------------------------

# diagnostics data frame
diag_df <- mod_sum_df %>% 
               left_join( rep_n_df ) %>% 
               # Add "flags" for convergence issues
               mutate( diverg   = n_diverg > 9,
                       rhat     = rhat_high > 0,
                       n_eff    = n_eff_low > 0,
                       mcse     = mcse_high > 0 ) %>% 
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
               mutate( sppcode  = paste0(sppcode1,sppcode2) %>% toupper ) %>% 
               # Correct model acronyms
               left_join( labs_df ) %>% 
               subset( model %in% mod_labs ) %>% 
               dplyr::select( -model ) %>% 
               rename( model = new_model ) %>% 
               mutate( model = factor(model, levels = new_labs) ) 

# focus on divergent transitions only
mod_div_df <- diag_df %>% 
                dplyr::select(model,sppcode,response,n_diverg) 

# models with no convergence issues
fitted_df   <- diag_df %>% 
                  subset( n_diverg  == 0 & 
                          rhat_high == 0 & 
                          n_eff_low == 0 & 
                          mcse_high == 0 )


# 2. Proportion models with issues ------------------------------------------


# proportion of issues by model

# proportions of ALL issues
diag_df %>% 
  select( species, model, clim, response, 
          diverg,  rhat,  n_eff ) %>% 
  mutate( issue = diverg + rhat + n_eff ) %>% 
  mutate( issue = issue > 0 ) %>% 
  group_by( model ) %>% 
  summarise( tot  = sum(issue),
             rep  = n() ) %>% 
  ungroup %>% 
  mutate( p_all = tot / rep )

# prop divergent transitions
div_df     <- diag_df %>% 
                group_by( model ) %>% 
                summarise( tot  = sum(diverg),
                           rep  = n() ) %>% 
                ungroup %>% 
                mutate( p_div = tot / rep )

# prop of issues with rhat
rhat_df     <- diag_df %>% 
                group_by( model ) %>% 
                summarise( tot  = sum(rhat),
                           rep  = n() ) %>% 
                ungroup %>% 
                mutate( p_rhat = tot / rep )

# prop mcse issue
mcse_df     <- diag_df %>% 
                group_by( model ) %>% 
                summarise( tot  = sum(mcse),
                           rep  = n() ) %>% 
                ungroup %>% 
                mutate( p_mcse = tot / rep )

# prop msce issue
n_eff_df    <- diag_df %>% 
                group_by( model ) %>% 
                summarise( tot  = sum(n_eff),
                           rep  = n() ) %>% 
                ungroup %>% 
                mutate( p_n_eff = tot / rep )


# proportions of problems by model
p1 <- ggplot(div_df, aes(model, p_div)) +
        geom_point(aes(size = '0.7') ) +
        ylim(0 ,1 ) +
        geom_text( aes(x=1.5,y=1,label='A)') ) +
        ylab( 'Prop. div issue' )+
        xlab( 'Model' ) +
        theme( axis.text.x  = element_text(angle = 90, 
                                           hjust = 1,vjust = 0.5), 
               axis.title      = element_text( size = 15),
               legend.position="none" )

p2 <- ggplot(rhat_df, aes(model, p_rhat)) +
        geom_point(aes(size = '0.7') ) +
        geom_text( aes(x=1.5,y=1,label='B)') ) +
        ylim(0 ,1 ) +
        ylab( 'Prop. rhat issue' )+
        xlab( 'Model' ) +
        theme( axis.text.x  = element_text(angle = 90, 
                                           hjust = 1,vjust = 0.5), 
               axis.title      = element_text( size = 15),
               legend.position="none" )

p3 <- ggplot(mcse_df, aes(model, p_mcse)) +
        geom_point(aes(size = '0.7') ) +
        geom_text( aes(x=1.5,y=1,label='C)') ) +
        ylim(0 ,1 ) +
        ylab( 'Prop. mcse issue' )+
        xlab( 'Model' ) +
        theme( axis.text.x  = element_text(angle = 90, 
                                           hjust = 1,vjust = 0.5), 
               axis.title      = element_text( size = 15),
               legend.position="none" )

p4 <- ggplot(n_eff_df, aes(model, p_n_eff)) +
        geom_point(aes(size = '0.7') ) +
        geom_text( aes(x=1.5,y=1,label='D)') ) +
        ylim(0 ,1 ) +
        ylab( 'Prop. n_eff issue' )+
        xlab( 'Model' ) +
        theme( axis.text.x  = element_text(angle = 90, 
                                           hjust = 1,vjust = 0.5), 
               axis.title      = element_text( size = 15),
               legend.position="none" )

out_p <- grid.arrange(p1, p2, p3, p4, nrow=2, ncol=2) 
ggsave(file = 'results/converg/prop_issue_mod.tiff',
       plot = out_p, width = 6.3, height = 6.3,
       compression="lzw")





# proportion of issues by model AND VITAL RATE 

# update the response variable
update <- data.frame( response = c('surv','grow','fec','log_lambda'),
                      Response = c('Survival',
                                   'Development',
                                   'Reproduction',
                                   'log(\u03BB)') ) %>% 
            mutate(   Response = factor(Response, levels = c('Survival',
                                                             'Development',
                                                             'Reproduction',
                                                             'log(\u03BB)') ) 
                  )
            

# prop divergent transitions
div_df     <- diag_df %>% 
                left_join( update ) %>% 
                group_by( model, Response ) %>% 
                summarise( tot  = sum(diverg),
                           rep  = n() ) %>% 
                ungroup %>% 
                mutate( p_div = tot / rep )

# prop divergent transitions
rhat_df     <- diag_df %>% 
                left_join( update ) %>% 
                group_by( model, Response ) %>% 
                summarise( tot  = sum(rhat),
                           rep  = n() ) %>% 
                ungroup %>% 
                mutate( p_rhat = tot / rep )

# prop mcse issue
mcse_df     <- diag_df %>% 
                left_join( update ) %>% 
                group_by( model, Response ) %>% 
                summarise( tot  = sum(mcse),
                           rep  = n() ) %>% 
                ungroup %>% 
                mutate( p_mcse = tot / rep )

# prop msce issue
n_eff_df    <- diag_df %>% 
                left_join( update ) %>% 
                group_by( model, Response ) %>% 
                summarise( tot  = sum(n_eff),
                           rep  = n() ) %>% 
                ungroup %>% 
                mutate( p_n_eff = tot / rep ) %>% 
                mutate( response = as.character(response) )

# proportions of problems by model
p1 <- ggplot(div_df, aes(model, p_div)) +
        geom_point(aes( color = Response),
                        size = 2,
                        position = position_jitter(w=0.1,h=0) )  +
        ylim(0 ,1 ) +
        ylab( 'Prop. div issue' )+
        xlab( 'Model' ) +
        scale_color_colorblind() +
        theme_bw() + 
        theme( axis.text.x  = element_text(angle = 90, 
                                           hjust = 1,
                                           vjust = 0.5), 
               axis.title      = element_text( size = 12) )

p2 <- ggplot(rhat_df, aes(model, p_rhat)) +
        geom_point(aes( color = Response),
                        size = 2,
                        position = position_jitter(w=0.1,h=0) ) +
        ylim(0 ,1 ) +
        ylab( 'Prop. rhat issue' )+
        xlab( 'Model' ) +
        scale_color_colorblind() +
        theme_bw() + 
        theme( axis.text.x  = element_text(angle = 90, 
                                           hjust = 1,vjust = 0.5), 
               axis.title   = element_text( size = 12),
               legend.position="none" )

p3 <- ggplot(mcse_df, aes(model, p_mcse)) +
        geom_point(aes( color = Response),
                        size = 2,
                        position = position_jitter(w=0.1,h=0) ) +
        ylim(0 ,1 ) +
        ylab( 'Prop. mcse issue' )+
        xlab( 'Model' ) +
        scale_color_colorblind() +
        theme_bw() + 
        theme( axis.text.x  = element_text(angle = 90, 
                                           hjust = 1,vjust = 0.5), 
               axis.title      = element_text( size = 12),
               legend.position="none" )

p4 <- ggplot(n_eff_df, aes(model, p_n_eff)) +
        geom_point(aes( color = Response ),
                        size = 2,
                        position = position_jitter(w=0.1,h=0) ) +
        ylim(0 ,1 ) +
        ylab( 'Prop. n_eff issue' )+
        xlab( 'Model' ) +
        scale_color_colorblind() +
        theme_bw() + 
        theme( axis.text.x  = element_text(angle = 90, 
                                           hjust = 1,vjust = 0.5), 
               axis.title      = element_text( size = 12),
               legend.position="none" )


plot_4x4 <- plot_grid(
  p1 + theme(legend.position="none" ), p2,
  p3, p4, 
  align = 'vh',
  labels = c("A", "B", "C", "D"),
  hjust = -0.5,
  nrow = 2
)

# produce legend
legend <- cowplot::get_legend(
  # create some space to the left of the legend
  p1 + theme(legend.box.margin = margin(0, 2, 0, 12))
)

out_p <- plot_grid(plot_4x4, legend, rel_widths = c(3, 0.7))

ggsave(file = 'results/converg/prop_issue_mod_vr.tiff',
       plot = out_p, width = 6.3, height = 5,
       compression="lzw")


# proportion of issues by replication of years 

# prop divergent transitions
div_df     <- diag_df %>% 
                group_by( rep_yr ) %>% 
                summarise( tot  = sum(diverg),
                           rep  = n() ) %>% 
                ungroup %>% 
                mutate( p_div = tot / rep )

# prop divergent transitions
rhat_df     <- diag_df %>% 
                group_by( rep_yr ) %>% 
                summarise( tot  = sum(rhat),
                           rep  = n() ) %>% 
                ungroup %>% 
                mutate( p_rhat = tot / rep )

# prop mcse issue
mcse_df     <- diag_df %>% 
                group_by( rep_yr ) %>% 
                summarise( tot  = sum(mcse),
                           rep  = n() ) %>% 
                ungroup %>% 
                mutate( p_mcse = tot / rep )

# prop msce issue
n_eff_df    <- diag_df %>% 
                group_by( rep_yr ) %>% 
                summarise( tot  = sum(n_eff),
                           rep  = n() ) %>% 
                ungroup %>% 
                mutate( p_n_eff = tot / rep )

# proportions of problems by model
p1 <- ggplot(div_df, aes(rep_yr, p_div)) +
        geom_point(aes(size = '0.7') ) +
        # ylim(0 ,1 ) +
        ylab( 'Prop. div issue' )+
        xlab( 'Number of years' ) +
        theme_bw() +
        theme( axis.text.x  = element_text(angle = 90, 
                                           hjust = 1,vjust = 0.5), 
               axis.title      = element_text( size = 15),
               legend.position="none" )

cor(div_df$rep_yr, div_df$p_div, 
    method = 'pearson') %>% 
  summary

p2 <- ggplot(rhat_df, aes(rep_yr, p_rhat)) +
        geom_point(aes(size = '0.7') ) +
        # ylim(0 ,1 ) +
        ylab( 'Prop. rhat issue' )+
        xlab( 'Number of years' ) +
        theme_bw() +
        theme( axis.text.x  = element_text(angle = 90, 
                                           hjust = 1,vjust = 0.5), 
               axis.title      = element_text( size = 15),
               legend.position="none" )

cor(rhat_df$rep_yr, rhat_df$p_rhat, 
    method = 'pearson') %>% 
  summary

p3 <- ggplot(mcse_df, aes(rep_yr, p_mcse)) +
        geom_point(aes(size = '0.7') ) +
        # ylim(0 ,1 ) +
        ylab( 'Prop. mcse issue' )+
        xlab( 'Number of years' ) +
        theme_bw() +
        theme( axis.text.x  = element_text(angle = 90, 
                                           hjust = 1,vjust = 0.5), 
               axis.title      = element_text( size = 15),
               legend.position="none" )

p4 <- ggplot(n_eff_df, aes(rep_yr, p_n_eff)) +
        geom_point(aes(size = '0.7') ) +
        # ylim(0 ,1 ) +
        ylab( 'Prop. n_eff issue' ) +
        xlab( 'Number of years' ) +
        theme_bw() +
        theme( axis.text.x  = element_text(angle = 90, 
                                           hjust = 1,vjust = 0.5), 
               axis.title      = element_text( size = 15),
               legend.position="none" )

cor(n_eff_df$rep_yr, n_eff_df$p_n_eff)

out_p <- grid.arrange(p1, p2, p3, p4, nrow=2, ncol=2) 
out_p <- plot_grid(
  p1 + theme(legend.position="none" ), p2,
  p3, p4, 
  align = 'vh',
  labels = c("A", "B", "C", "D"),
  hjust = -0.5,
  nrow = 2
)


ggsave(file = 'results/converg/prop_issue_rep_yr.tiff',
       plot = out_p, width = 6.3, height = 5,
       compression="lzw")


       
# 3. Proportion of parameters with issues ----------------------

# NOT DOABLE WITH 

# get all names of diagnostics
par_n <- diag_df %>% names 

# model parameters (NOT yhat, log_lik, etc.)
params<- Filter(function(x) !grepl('yhat|median|log_lik',x),
                par_n)[c(3,12:67)] 

# numbers of parameters per model
diag_df %>% 
  select(species, model, response, clim) %>% 
  mutate(n_pars = apply(diag_df[,params],1,function(x) sum(!is.na(x)) ) ) %>% 
  select(model,n_pars) %>% 
  unique

# this clearly includes non-model params. Need to work on raw data.
diag_df$rhat_high %>% unique



# 4. Rhat ~ n_eff ----------------------------------------

# data frame containing inputs
input_df  <- expand.grid( resp     = c('log_lambda','surv',
                                       'grow','fec'),
                          clim_var = c('precip','airt'),
                          stringsAsFactors = F ) %>% 
                mutate( resp_clim = paste(resp, clim_var, sep='_') )

# file list
sum_f_l   <- grep('diagnostics_',
                  list.files(mod_dir),
                  value=T) #%>% 
               # grep(resp_clim, ., value=T)

# path list
path_l    <- paste0(mod_dir, sum_f_l)

# extract species names
spp_names <- gsub("diagnostics_", "", sum_f_l ) %>%
                gsub('array_vr-[0-9]{7}-[0-9]{1,3}-[0-9]{1,2}_', "", . ) %>%
                gsub('.csv', "", . ) %>%
                gsub( paste(input_df$resp_clim, collapse='|'),
                      '', .) %>%
                gsub('_$','', .)

# species, response, and climate
spp_resp_clim <- gsub("diagnostics_", "", sum_f_l ) %>%
                   gsub('array_vr-[0-9]{7}-[0-9]{1,3}-[0-9]{1,2}_', "", . ) %>%
                   gsub('.csv', "", . )

# get response and climate variable (without spp info)
get_res_clim <- function(ii){
  gsub(spp_names[ii], '', spp_resp_clim[ii]) %>% 
    gsub('^_','',.)
}

# response and climate
resp_clim <- sapply(1:length(spp_resp_clim), get_res_clim)

# download files serially
create_sum_df <- function(path_l, 
                          spp_names,
                          resp_clim){
  
  read.csv(path_l, stringsAsFactors = F) %>% 
    tibble::add_column(species = spp_names,.before=1) %>% 
    tibble::add_column(resp_clim = resp_clim, .before=4)
  
}

# a common format for both diagnost_df and plot_l_df
format_species_names <- function( x ){
  
  x %>% 
    # remove "double" populations
    subset( !(species %in% c('Astragalus_scaphoides_6_site_rep',
                             'Astragalus_scaphoides_2')) ) %>% 
    mutate( species   = gsub('_',' ',species) ) %>% 
    mutate( species = replace(species,
                              grepl('Eriogonum',species),
                              'Eriogonum longifolium...') ) %>% 
    mutate( species = replace(species,
                              grepl('Astragalus scaphoides 6 long',species),
                              'Astragalus scaphoides') ) %>% 
    mutate( species = replace(species,
                              grepl('Cirsium pitcheri 8',species),
                              'Cirsium pitcheri (1)') ) %>% 
    mutate( species = replace(species,
                              grepl('Cirsium pitcheri 4',species),
                              'Cirsium pitcheri (2)') ) %>% 
    mutate( species = replace(species,
                              grepl('Cirsium pitcheri 6',species),
                              'Cirsium pitcheri (3)') ) %>% 
    # remove the "floating 2" left
    mutate( species = gsub(' 2','',species) ) 
  
}


# get all model summaries
diagnost_df  <- Map(create_sum_df, path_l,
                    spp_names, resp_clim) %>% 
                  bind_rows %>% 
                  mutate( resp_clim = gsub('log_lambda', 'loglambda', resp_clim) ) %>% 
                  separate( resp_clim, c('response','clim'),sep='_') %>% 
                  mutate( response = gsub('loglambda', 'log_lambda', response) ) %>%
                  subset( !(model %in% c('gaus','simpl_n','ridge')) ) %>% 
                  mutate( model = factor(model, levels=mod_ordr) ) %>% 
                  format_species_names %>% 
                  # Introduce "correct" model names
                  left_join( labs_df ) %>% 
                  subset( model %in% mod_labs ) %>% 
                  dplyr::select( -model ) %>% 
                  rename( model = new_model ) %>% 
                  mutate( model = factor(model, levels = c(new_labs[2:13],new_labs[1])) )

# Plot

# data frame containing inputs
plot_l_df  <- expand.grid( resp     = c('log_lambda','surv',
                                       'grow','fec'),
                           clim_var = c('precip','airt'),
                           species  = spp,
                           stringsAsFactors = F ) %>% 
                  mutate( resp_clim = paste(resp, clim_var, sep='_') ) %>% 
                  format_species_names
                  

# plot rhat vs. neff for each species/climate variable/response BY SPECIES
plot_rhat_neff_spp <- function(ii){
  
  print(ii)
  
  plot_df <- diagnost_df %>% 
    subset( response == plot_l_df$resp[ii] &
              clim     == plot_l_df$clim_var[ii] &
              species  == plot_l_df$species[ii] )
  
  title_spp <- plot_df$species %>% unique
  
  if( nrow(plot_df) > 1 ){
    
    plot_out <- 
      ggplot(plot_df) +
      geom_point( aes(n_eff,Rhat) ) +
      facet_wrap( ~ model, ncol = 3) + 
      labs( title = title_spp,
            x     = expression('N'['eff']),
            y     = expression(hat(R)) ) +
      theme_bw() +
      theme( strip.text.y  = element_text( size = 10,
                                           margin = margin(0.5,0.5,0.5,0.5,
                                                           'mm') ),
             strip.text.x  = element_text( size = 10,
                                           margin = margin(0.5,0.5,0.5,0.5,
                                                           'mm') ),
             strip.switch.pad.wrap = unit('0.5',unit='mm'),
             panel.spacing = unit('0.5',unit='mm'),
             plot.title    = element_text( hjust = 0.5,
                                           size  = 20),
             axis.title    = element_text( vjust = 1,
                                           hjust = 0.5,
                                           size  = 17)
            ) 
    
      # store plot
      ggsave( paste0(getwd(),
                    '/results/converg/rhat_neff/',
                    substr(unique(plot_df$clim),1,4),'/',
                    unique(plot_df$response),'/',
                    unique(plot_df$species),'_neff_rhat',
                    '.png'),
             plot = plot_out,
             width=6.3,height=9)
  }
  
}

# Plot those out
lapply(1:nrow(plot_l_df), plot_rhat_neff_spp)



# plot rhat vs. neff ACROSS species
plot_rhat_neff <- function(ii){
  
  plot_df <- diagnost_df %>% 
                subset( response == input_df$resp[ii] &
                        clim     == input_df$clim_var[ii]  )
  
  # big plot across all models
  ggplot(plot_df) +
    geom_point( aes(n_eff,
                    Rhat),
                alpha = 0.2) +
    geom_hline( aes(yintercept = 1.1), 
                lty = 2) +
    ylim(1, 1.5) +
    facet_wrap( ~ model) + 
    theme( strip.text.y  = element_text( size = 20,
                                     margin = margin(0.5,0.5,0.5,0.5,
                                                     'mm') ),
           strip.text.x  = element_text( size = 10,
                                         margin = margin(0.5,0.5,0.5,0.5,
                                                         'mm') ),
           strip.switch.pad.wrap = unit('0.5',unit='mm'),
           panel.spacing = unit('0.5',unit='mm') 
           ) +
    ggsave(paste0('results/converg/rhat_neff/',
                  substr(unique(plot_df$clim),1,4),'/',
                  unique(plot_df$response),'_',
                  'neff_rhat.tiff'),
           width=6.3,height=9,compression='lzw')

}

# store all Rhat vs. n_eff across all species
lapply(1:nrow(input_df), plot_rhat_neff)


# PRECIPITATION: big plot across all models and vital rates
diagnost_df %>% 
  subset( clim == 'precip' ) %>% 
  ggplot() +
  geom_point( aes(n_eff,
                  Rhat,
                  color = response),
              alpha = 0.2) +
  geom_hline( aes(yintercept = 1.1), 
              lty = 2) +
  ylim(1, 1.15) +
  facet_wrap( ~ model) + 
  scale_color_viridis_d() +
  ylab( expression(hat(R)) ) +
  xlab( 'Effective sample size' ) +
  theme( strip.text.y  = element_text( size = 20,
                                   margin = margin(0.5,0.5,0.5,0.5,
                                                   'mm') ),
         strip.text.x  = element_text( size = 10,
                                       margin = margin(0.5,0.5,0.5,0.5,
                                                       'mm') ),
         strip.switch.pad.wrap = unit('0.5',unit='mm'),
         panel.spacing = unit('0.5',unit='mm') 
         ) 
  ggsave(paste0('results/converg/rhat_neff/prec_neff_rhat.tiff'),
         width=6.3,height=6.3,compression='lzw')

# AIR TEMPERATURE: big plot across all models and vital rates
diagnost_df %>% 
  subset( clim == 'airt' ) %>% 
  ggplot() +
  geom_point( aes(n_eff,
                  Rhat,
                  color = response),
              alpha = 0.2) +
  geom_hline( aes(yintercept = 1.1), 
              lty = 2) +
  ylim(1, 1.15) +
  facet_wrap( ~ model) + 
  scale_color_viridis_d() +
  ylab( expression(hat(R)) ) +
  xlab( 'Effective sample size' ) +
  theme( strip.text.y  = element_text( size = 20,
                                   margin = margin(0.5,0.5,0.5,0.5,
                                                   'mm') ),
         strip.text.x  = element_text( size = 10,
                                       margin = margin(0.5,0.5,0.5,0.5,
                                                       'mm') ),
         strip.switch.pad.wrap = unit('0.5',unit='mm'),
         panel.spacing = unit('0.5',unit='mm') 
         ) +
  ggsave('results/converg/rhat_neff/airt_neff_rhat.tiff',
         width=6.3,height=6.3,compression='lzw')


# ACROSS temperature and precipitation 
diagnost_df %>% 
  left_join( update ) %>%
  # left_join( leg_lab_df ) %>% 
  ggplot() +
  geom_point( aes(n_eff,
                  Rhat,
                  color = Response),
              alpha = 0.4) +
  geom_hline( aes(yintercept = 1.1), 
              lty = 2) +
  ylim(1, 1.15) +
  facet_wrap( ~ model, nrow = 6) + 
  scale_color_viridis_d() +
  labs( y     = expression(hat(R)),
        x     = 'Effective sample size',
        color = 'Response') +
  theme_bw() +
  theme( strip.text.y  = element_text( size = 20,
                                   margin = margin(0.5,0.5,0.5,0.5,
                                                   'mm') ),
         strip.text.x  = element_text( size = 10,
                                       margin = margin(0.5,0.5,0.5,0.5,
                                                       'mm') ),
         strip.switch.pad.wrap = unit('0.5',unit='mm'),
         strip.background = element_blank(),
         panel.spacing = unit('0.5',unit='mm'),
         plot.margin = margin( r = 2, l = 0.5 ),
         axis.title.y = element_text( angle = 0, vjust = 0.5 ),
         axis.text.x  = element_text( angle = 60, hjust = 1 ) ) +
  ggsave( 'results/converg/rhat_neff/neff_rhat.tiff',
          width=6.3,height=8,compression='lzw' )





# plots ------------------------------------------------------------
prop_div <- diag_df %>% 
              group_by( model ) %>% 
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
