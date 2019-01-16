# produce summary graphs
source("code/format_data.R")
library(tidyverse)
library(magrittr)
library(testthat)
library(viridis)
library(egg)
options(stringsAsFactors = F )

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
                  # selec ONLY these model selection variables
                  lapply(function(x) x %>% 
                                      dplyr::select(species, 
                                                     model, 
                                                     waic, 
                                                     looic, 
                                                     mse, 
                                                     elpd)) %>% 
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
    mutate( mean_elpd = elpd / rep_n,
            clim_var  = clim_var,
            response  = response )

}

# all models results
input_df    <- expand.grid( clim_var = c("precip","airt"),
                            response = c("surv","grow","fec","log_lambda"),
                            interval = "",
                            stringsAsFactors = F)

# ALL model information
mod_perf_df  <- lapply(1:nrow(input_df), mod_perform) %>%
                  bind_rows %>% 
                  arrange(species, mse)


# keep only species that are always present -------------------------------------------

# species that are always present
get_spp <- function(ii){
  mod_perf_df %>% 
    subset( clim_var == input_df$clim_var[ii] & 
            response == input_df$response[ii]) %>% 
    .$species %>% 
    unique
}

# list of vectors
spp_l_v <- lapply(1:nrow(input_df), get_spp )
spp_all <- spp_l_v %>% unlist %>% unique

# test all combinations 
all_intersect <- function(spp_all, list_vec){
  sapply(spp_all,
         function(x) sapply(list_vec, 
                            function(y) x %in% y) 
  )
}

# vector of species common to all results files
spp_common_v <- all_intersect(spp_all,spp_l_v) %>% 
                  apply(2,sum) %>% 
                  Filter( function(x) x == 8, .) %>% 
                  names

# keep spp for which data is available for every response X clim var combination
mod_perf_df  <- mod_perf_df %>% 
                  subset( species %in% spp_common_v )


# Multiple Tileplots: differences in absolute lppd ---------------------------------------

# species order - based on rep_yr, rep_n, species, model... 
spp_df        <- mod_perf_df %>% 
                    dplyr::select(species, model, rep_yr, rep_n) %>% 
                    unique %>% 
                    arrange(rep_yr, rep_n, species, model)

# labels for models to plot on x-axis
mod_labs    <- c('ctrl1','ctrl2','yr1','yr2','yr3','yr_wgt',
                 'yr_bet', 'gaus','expp','gev',
                 'expp_n','gev_n','simpl_n')

# resp    <- 'surv'
# clim_v  <- 'airt'

# format the differences in lppd
form_diff_lppd_df <- function(resp, clim_v){
  
  perf_by_spp <- mod_perf_df %>%
                    subset( response == resp & clim_var == clim_v ) %>%
                    dplyr::select(species, model, rep_yr, rep_n, elpd) %>% 
                    mutate( elpd = replace(elpd, elpd == -Inf, NA)) %>%
                    spread(model, elpd)
  perf_mat    <- perf_by_spp %>% 
                    dplyr::select(-species, -rep_yr, -rep_n) %>% t
  
  # get minimum of maximum values
  get_max     <- function(ii) which(perf_mat[,ii] == max(perf_mat[,ii],na.rm=T))
  
  # species sequence
  spp_seq     <- 1:ncol(perf_mat)
  
  # matrix of benchmark ids
  max_ids     <- sapply(spp_seq, get_max) %>%
                        cbind( spp_seq )
  
  # benchmark values
  max_val     <- perf_mat[max_ids]
  
  # differences matrix
  diff_mat    <- sweep(perf_mat, 2,  max_val, FUN='-') 
  
  # differences data frame
  long_df <- diff_mat %>%
                as.data.frame %>%
                setNames( perf_by_spp$species ) %>%
                tibble::add_column(model = row.names(.), .before=1) %>%
                gather(species,elpd,Actaea_spicata:Thelesperma_megapotamicum) %>% 
                full_join( spp_df ) %>% 
                mutate( model = factor(model, levels = mod_labs) ) %>%
                left_join( dplyr::select(perf_by_spp, species, rep_yr, rep_n) ) %>% 
                arrange( rep_yr, rep_n, species, model )  %>% 
                mutate( species = replace(species,
                                          grepl('Eriogonum',species),
                                          'Eriogonum_longifolium...') ) %>% 
                mutate( species = factor(species, levels = unique(species)) )
  
  spp_v  <- long_df$species %>% unique %>% as.character
  
  # find_best_mod <- function(x){
  #   long_df %>%
  #     subset( species == x) %>%
  #     subset( elpd == max(elpd) ) %>%
  #     select( model, species )
  # }
  
  model_rank <- function(x){
    long_df %>%
      subset( species == x) %>%
      arrange( desc(elpd) ) %>% 
      mutate( mod_rank = c(1:nrow(.)) ) %>% 
      mutate( mod_rank = replace(mod_rank,
                                 mod_rank > 3,
                                 NA) ) %>% 
      mutate( mod_rank = as.character(mod_rank) ) %>% 
      select( model, species, mod_rank )
  }

  # pick the best models
  mod_sel_df <- lapply(spp_v, model_rank) %>% 
                    bind_rows %>% 
                    right_join( long_df )
  
  return(mod_sel_df)
    
}


form_diff_lppd_df('log_lambda', 'precip') %>% 
  subset( grepl('Eryngium',species) ) %>% head

# four tile plot 
four_tile_plot <- function(format_function, fill_var, clim_v, var_lim, 
                           expr_fill, file_title){
  
  # plot it out
  p1 <- ggplot(format_function('surv', clim_v), aes(model, species)) +
        geom_tile(aes_string(fill = fill_var), color = "white") +
        geom_point(aes(size = '0.7',
                       shape = mod_rank) ) +
        scale_fill_viridis( guide = F, limits = var_lim ) + #
        ggtitle('Survival') +
        theme(title        = element_text(angle = 0, hjust = 0.5, size = 10), 
              legend.title = element_text(size = 10),
              legend.text  = element_text(size = 12),
              axis.title   = element_blank(),
              legend.position ="none",
              axis.text.x  = element_text(angle = 90, hjust = 1),
              legend.key   = element_blank()) +
        labs(fill = element_blank() )
    
  p2 <- ggplot(format_function('grow', clim_v), aes(model, species)) +
        geom_tile(aes_string(fill = fill_var), color = "white") +
        geom_point(aes(size = '0.7',
                       shape = mod_rank) ) +
        scale_fill_viridis( limits = var_lim ) + #
        guides(fill = F ) +
        ggtitle('Growth') +
        theme(title        = element_text(angle = 0, hjust = 0.5, size = 10), 
              legend.title = element_blank(), 
              legend.position ="none",
              axis.text.y  = element_blank(), 
              legend.text  = element_blank(),
              axis.title   = element_blank(),
              axis.text.x  = element_text(angle = 90, hjust = 1),
              legend.key   = element_blank())
  
  p3 <- ggplot(format_function('fec',clim_v), aes(model, species)) +
        geom_tile(aes_string(fill = fill_var), color = "white") +
        geom_point(aes(size = '0.7',
                       shape = mod_rank) ) +
        scale_fill_viridis( limits = var_lim ) + #
        guides(fill = F ) +
        ggtitle('Fecundity') +
        theme(title        = element_text(angle = 0, hjust = 0.5, size = 10), 
              legend.title = element_blank(), 
              axis.text.y  = element_blank(), 
              legend.text  = element_blank(),
              legend.position ="none",
              axis.title   = element_blank(),
              axis.text.x  = element_text(angle = 90, hjust = 1),
              legend.key   = element_blank()) +
        labs(fill = element_blank() )
  
  p4 <- ggplot(format_function('log_lambda',clim_v), aes(model, species)) +
        geom_tile(aes_string(fill = fill_var), color = "white") +
        geom_point(aes(size = '0.7',
                       shape = mod_rank) ) +
        scale_fill_viridis( limits = var_lim ) + #
        ggtitle('Log Lambda') +
        theme(title        = element_text(angle = 0, hjust = 0.5, size = 10),
              legend.title = element_text(size = 10),
              legend.text = element_text(size = 12),
              axis.title = element_blank(),
              axis.text.y  = element_blank(),
              axis.text.x = element_text(angle = 90, hjust = 1)) +
        labs(fill = expr_fill )
  
  # plot 
  tiff( file_title, 
        unit="in", width=10, height=6.3, res=600,compression="lzw" )
  grid.arrange(nrow = 1,
               grobs = list(p1, p2, p3, p4),
               widths = c(1.8, 0.9, 0.9, 1.4),
               layout_matrix = rbind(c(1, 2, 3, 4)) )
  dev.off()

}

# differences in absolute lppd plots
four_tile_plot(form_diff_lppd_df, 'elpd', 'airt',    c(-30,0), 
               expression(Delta*" LPPD"), 'results/lppd_diff_vr_airt.tiff')
four_tile_plot(form_diff_lppd_df, 'elpd', 'precip',  c(-30,0), 
               expression(Delta*" LPPD"), 'results/lppd_diff_vr_precip.tiff')


# Code below is no more current -------------------------------------------------------

# format_function <- form_diff_lppd_df
# fill_var    <- 'elpd'
# clim_v      <- 'airt'
# var_lim     <- c(-30,0)
# expr_fill   <- expression(Delta*" LPPD")
# file_title  <- 'results/lppd_diff_vr_airt.tiff'
# 
# 
# # mean lppd
# range_m_lppd  <- mod_perf_df$mean_elpd %>% range(na.rm=T)
# four_tile_plot(form_mean_lppd_df, 'mean_elpd', 'airt',   c(-15,2), #range_m_lppd,
#                "Mean LPPD",       'mean_lppd_vr_airt.tiff')
# four_tile_plot(form_mean_lppd_df, 'mean_elpd', 'precip', c(-15,2), #range_m_lppd, 
#                "Mean LPPD",        'mean_lppd_vr_precip.tiff')
# 
# 
# # rank tile plots
# four_tile_plot(form_rank_lppd_df, 'rank', 'airt',   c(1,13),  
#                'Model rank',    'rank_lppd_vr_airt.tiff')
# four_tile_plot(form_rank_lppd_df, 'rank', 'precip', c(1,13),  
#                'Model rank',    'rank_lppd_vr_precip.tiff')
# 
# 
# # format the data as we want
# form_mean_lppd_df <- function(resp, clim_v){
#   
#   # data for specific response and climate variable
#   elpd_mod_df <- mod_perf_df %>%
#                     subset( response == resp & clim_var == clim_v ) %>%
#                     dplyr::select(species, model, rep_yr, rep_n, mean_elpd, mod_rank) %>% 
#                     mutate( mean_elpd = replace(mean_elpd, mean_elpd == -Inf, NA) )
#       
#   # long df with model ranks!
#   mean_elpd_df <- elpd_mod_df %>% 
#                     full_join( spp_df ) %>% 
#                     mutate( model = factor(model, levels = mod_labs) ) %>%
#                     arrange( rep_yr, rep_n, species, model )  %>% 
#                     mutate( species = replace(species,
#                                               grepl('Eriogonum',species),
#                                               'Eriogonum_longifolium...') ) %>% 
#                     mutate( species = factor(species, levels = unique(species)) ) 
#   
#   return(mean_elpd_df)
#     
# }
# 
# 
# 
# # format the data as we want
# form_rank_lppd_df <- function(resp, clim_v){
#   
#   # data for specific response and climate variable
#   elpd_mod_df <- mod_perf_df %>%
#                     subset( response == resp & clim_var == clim_v ) %>%
#                     dplyr::select(species, model, rep_yr, rep_n, elpd) %>% 
#                     mutate( elpd = replace(elpd, elpd == -Inf, NA) ) %>% 
#                     arrange(species, elpd )
#       
#   # split data by species (to ricompose later)
#   spp_mod_l   <- split(elpd_mod_df, elpd_mod_df$species)
#   
#   # long df with model ranks!
#   mod_rank_df <- lapply(spp_mod_l, function(x) mutate(x, rank = 1:13) ) %>% 
#                     # ricompose in one data frame!
#                     bind_rows %>% 
#                     full_join( spp_df ) %>% 
#                     mutate( model = factor(model, levels = mod_labs) ) %>%
#                     arrange( rep_yr, rep_n, species, model )  %>% 
#                     mutate( species = replace(species,
#                                               grepl('Eriogonum',species),
#                                               'Eriogonum_longifolium...') ) %>% 
#                     mutate( species = factor(species, levels = unique(species)) ) %>% 
#                     mutate( rank = replace(rank, is.na(elpd), NA) )
#   
#   return(mod_rank_df)
#     
# }
# 