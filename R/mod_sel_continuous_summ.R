# produce summary graphs
source("R/format_data.R")
library(tidyverse)
library(testthat)
library(viridis)
library(egg)
library(ggthemes)
library(gridExtra)
library(ggpubr)
library(quantreg)

options(stringsAsFactors = F )

today_date    <- gsub("-","_",Sys.time() %>% as.Date %>% as.character)
m_back        <- 36
interval      <- NULL
pre_chelsa    <- NULL # '_pre_chelsa'


# quote a series of bare names
quote_bare <- function( ... ){
  substitute( alist(...) ) %>% 
    eval( ) %>% 
    sapply( deparse )
}


# Summarize moving windows results by climate variable -------------------------------
mod_perform <- function(ii, res_folder){
  
  clim_var <- input_df$clim_var[ii]
  response <- input_df$response[ii]
  interval <- input_df$interval[ii]
  resp_clim<- paste0("_",response,"_",clim_var)
  
  # read lambda/clim data 
  lam      <- read.csv("data/all_demog_updt.csv", 
                       stringsAsFactors = F) #%>%
                #subset( SpeciesAuthor != "Purshia_subintegra" )

  # summary info 
  
  # summary files names
  sum_files <- list.files(res_folder)[grep("mod_summaries_", list.files(res_folder) )] %>% 
                  stringr::str_subset(resp_clim)
  
  # read files
  mod_summ  <- lapply(sum_files, function(x) read.csv(paste0(res_folder,"/",x)) ) %>%
                  setNames( gsub("mod_summaries_", "", sum_files ) ) %>%
                  setNames( gsub(paste0(resp_clim,".csv"), "", names(.) ) ) %>% 
                  setNames( gsub('array_vr-[0-9]{7}-[0-9]{1,3}-[0-9]{1,3}_', "", names(.) ) )
  
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
                                                     elpd)
                          ) %>% 
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

# model results in 
# dir          <- 'C:/Users/ac22qawo/sapropos_main/out_22.4.2020'
# dir          <- 'C:/Users/ac22qawo/sapropos_main/out_2021'
# dir          <- 'C:/Users/ac22qawo/sapropos_main/out_2021.2'
dir          <- 'C:/Users/ac22qawo/sapropos_main/out_2021.8.12'

# ALL model information
mod_perf_df  <- lapply(1:8, mod_perform, dir) %>%
                  # lapply(1:nrow(input_df), mod_perform, dir) %>%
                  bind_rows %>% 
                  arrange(species, mse) %>% 
                  # PROVISIONAL: remove mistake "simpl"
                  subset( !(model %in% 'simpl') )


# check missing simulations ------------------------------------------------------

# species 
spp          <- read.csv("data/all_demog_updt.csv", 
                         stringsAsFactors = F) %>% 
                  .$SpeciesAuthor %>% 
                  unique

# family
family_df    <- data.frame( response = c("surv","grow","fec","log_lambda"),
                            family   = c('beta','beta','gamma','normal'),
                            stringsAsFactors = F )

# species df
spp_df       <- data.frame( species = spp,
                            ii      = 1:39 )

# missing species and numbers
spp_all_df   <- expand.grid( clim_var = c("precip","airt"),
                             response = c("surv","grow","fec","log_lambda"),
                             species  = spp,
                             stringsAsFactors = F) %>% 
                  arrange( clim_var, response, species ) %>% 
                  mutate( job_n = 1:nrow(.) )

# Missing simulations
miss_df      <- anti_join( spp_all_df, mod_perf_df ) %>% 
                  full_join( family_df ) %>% 
                  subset( !grepl('Daphne',species) )

# write.csv(miss_df,
#           'C:/Users/ac22qawo/sapropos_main/miss_sims.csv',
#           row.names=F)


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
                  Filter( function(x) x == nrow(input_df), .) %>% 
                  names

setdiff(spp_all, spp_common_v)


# keep spp for which data is available for every response X clim var combination
mod_perf_df  <- mod_perf_df %>% 
                  subset(   species %in% spp_common_v ) %>% 
                  subset( !(species %in% c('Astragalus_scaphoides_6_site_rep',
                                         'Astragalus_scaphoides_2')) )


# Multiple Tileplots: differences in absolute lppd ---------------------------------------

# labels for models to plot on x-axis
mod_labs    <- quote_bare( ctrl1, 
                           yr1,    yr2,     yr3,
                           gaus1,  gaus2,   gaus3,
                           simpl1, simpl2,  simpl3,
                           ridge1, ridge2,  ridge3 )
                           # gaus,   simpl_n, ridge )

new_labs    <- c( 'NM', 
                  'CSM 1', 'CSM 2', 'CSM 3',
                  'WMM 1', 'WMM 2', 'WMM 3',
                  'SAM 1', 'SAM 2', 'SAM 3',
                  'FHM 1', 'FHM 2', 'FHM 3' )

labs_df     <- data.frame( model     = mod_labs,
                           new_model = new_labs )

# mod_labs    <- quote_bare( ctrl1, 
#                            yr1,   gaus1, simpl1, ridge1,
#                            yr2,   gaus2, simpl2, ridge2,
#                            yr3,   gaus3, simpl3, ridge3 )

mod_no_r    <- quote_bare( ctrl1, 
                           yr1,    yr2,     yr3,
                           gaus1,  gaus2,   gaus3,
                           simpl1, simpl2,  simpl3 )
                           # gaus,   simpl_n )

  
mod_perf_df  <- mod_perf_df %>% 
                  subset( model %in% mod_labs ) %>% 
                  left_join( labs_df ) %>% 
                  dplyr::select( -model ) %>% 
                  rename( model = new_model )

mod_no_r_df  <- mod_perf_df %>% 
                  subset( model %in% mod_no_r )

# species order - based on rep_yr, rep_n, species, model... 
spp_df        <- mod_perf_df %>% 
                  dplyr::select(species, model, rep_yr, rep_n) %>% 
                  unique %>% 
                  arrange(rep_yr, rep_n, species, model)

clim_v  = 'precip'

# format the differences in lppd
form_diff_lppd_df <- function(resp, clim_v, mod_df){
  
  perf_by_spp <- mod_df %>%
                    subset( response == resp & clim_var == clim_v ) %>%
                    dplyr::select(species, model, rep_yr, rep_n, elpd) %>% 
                    mutate( elpd = replace(elpd, elpd == -Inf, NA)) %>%
                    spread(model, elpd) %>% 
                    mutate( species = gsub('_',' ',species) )
  
  spp_df      <- spp_df %>% 
                    mutate( species = gsub('_',' ',species) )
  
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
                gather(species,elpd, setdiff(names(.),'model') ) %>% 
                full_join( spp_df ) %>% 
                mutate( model = factor(model, levels = new_labs) ) %>%
                left_join( dplyr::select(perf_by_spp, species, rep_yr, rep_n) ) %>% 
                arrange( rep_yr, rep_n, species, model )  %>% 
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
                mutate( species = gsub(' 2','',species) ) %>% 
                mutate( species = factor(species, levels = unique(species)) )
  
  spp_v  <- long_df$species %>% unique %>% as.character
  
  # find_best_mod <- function(x){
  #   long_df %>%
  #     subset( species == x) %>%
  #     subset( elpd == max(elpd) ) %>%
  #     select( model, species )
  # }
  
  # model rank
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


# # test --------------------------------------------------------
# resp    <- 'surv'
# clim_v  <- 'prec'
# mod_df  <- mod_perf_df
# 
# all_demo = form_diff_lppd_df(resp,'precip', mod_perf_df) %>% 
#               subset( grepl('Eryngium_cuneifolium',species ) )
# 
# nor_demo = form_diff_lppd_df(resp,'precip', mod_no_r_df) %>% 
#               subset( grepl('Eryngium_cuneifolium',species ) )
# 
# all_demo
# nor_demo
# 
# subset(mod_perf_df, species == 'Eryngium_cuneifolium') %>% 
#   subset( clim_var == 'precip' ) %>% 
#   subset( response == 'grow' ) %>% 
#   dplyr::select( model, elpd, clim_var, response )

# ----------------------------------------------------

# four tile plot 
four_tile_plot <- function( fill_var,  clim_v, 
                            mod_df, var_lim, expr_fill, file_title){
  
  # plot it out
  p1 <- ggplot(form_diff_lppd_df('surv', clim_v, mod_df), aes(model, species)) +
        geom_tile(aes_string(fill = fill_var), color = "white") +
        geom_point(aes(shape = mod_rank),
                   size = 2) +
        geom_vline( xintercept = c(1.5,4.5,7.5,10.5),
                    col = 'red', lwd=0.5, lty=2) + 
        scale_fill_viridis( guide = F, limits = var_lim ) + #
        ggtitle('A) Survival') +
        theme(title        = element_text(angle = 0, hjust = 0.5, size = 10), 
              legend.title = element_text(size = 10),
              legend.text  = element_text(size = 12),
              axis.title   = element_blank(),
              legend.position ="none",
              axis.text.y  = element_text(face = 'italic'),
              axis.text.x  = element_text(angle = 90, 
                                          hjust = 1,
                                          vjust = 0.5),
              legend.key   = element_blank()) +
        labs(fill = element_blank() )
    
  p2 <- ggplot(form_diff_lppd_df('grow', clim_v, mod_df), aes(model, species)) +
        geom_tile(aes_string(fill = fill_var), color = "white") +
        geom_point(aes(shape = mod_rank),
                   size = 2) +
        geom_vline( xintercept = c(1.5,4.5,7.5,10.5),
                    col = 'red', lwd=0.5, lty=2) + 
        scale_fill_viridis( limits = var_lim ) + #
        ggtitle('B) Development') +
        guides(fill = F ) +
        theme(title        = element_text(angle = 0, hjust = 0.5, size = 10), 
              legend.title = element_blank(), 
              legend.position ="none",
              axis.text.y  = element_blank(),
              legend.text  = element_blank(),
              axis.title   = element_blank(),
              axis.text.x  = element_text(angle = 90, 
                                          hjust = 1,
                                          vjust = 0.5),
              legend.key   = element_blank())
  
  p3 <- ggplot(form_diff_lppd_df('fec',clim_v, mod_df), aes(model, species)) +
        geom_tile(aes_string(fill = fill_var), color = "white") +
        geom_point(aes(shape = mod_rank),
                   size = 2) +
        geom_vline( xintercept = c(1.5,4.5,7.5,10.5),
                    col = 'red', lwd=0.5, lty=2) + 
        scale_fill_viridis( limits = var_lim ) + #
        ggtitle('C) Reproduction') +
        guides(fill = F ) +
        theme(title        = element_text(angle = 0, hjust = 0.5, size = 10), 
              legend.title = element_blank(), 
              axis.text.y  = element_blank(),
              legend.text  = element_blank(),
              legend.position ="none",
              axis.title   = element_blank(),
              axis.text.x  = element_text(angle = 90, 
                                          hjust = 1,
                                          vjust = 0.5),
              legend.key   = element_blank()) +
        labs(fill = element_blank() )
  
  p4 <- ggplot(form_diff_lppd_df('log_lambda',clim_v, mod_df), aes(model, species)) +
        geom_tile(aes_string(fill = fill_var), color = "white") +
        geom_point(aes(shape = mod_rank),
                   size = 2) +
        geom_vline( xintercept = c(1.5,4.5,7.5,10.5),
                     col = 'red', lwd=0.5, lty=2) +
        scale_fill_viridis( limits = var_lim ) +
        ggtitle( expression('D) Log('*lambda*')') ) +
        theme(title        = element_text(angle = 0, hjust = 0.5, size = 10),
              legend.title = element_text(size = 10),
              legend.text = element_text(size = 12),
              axis.title = element_blank(),
              axis.text.y  = element_blank(),
              axis.text.x = element_text(angle = 90,
                                         hjust = 1,
                                         vjust = 0.5) ) +
        labs(fill = expr_fill ) +
        scale_shape_discrete( labels = c('1','2','3','>3') )

  # plot 
  tiff( file_title, 
        unit="in", width=10, height=6.3, res=600,compression="lzw" )
  grid.arrange(nrow   = 1,
               # grobs  = list(p1, p2, p3),
               # widths = c(1.8, 0.9, 1.4),
               # layout_matrix = rbind(c(1, 2, 3))
               grobs  = list(p1, p2, p3, p4),
               widths = c(1.9, 0.9, 0.9, 1.4),
               layout_matrix = rbind(c(1, 2, 3, 4))
               )
  dev.off()

}

# differences in absolute lppd plots
four_tile_plot('elpd', 'airt', mod_perf_df,
               c(-20,0),  expression(Delta*" LPPD"), 
               'results/lppd_diff_airt_2022.2_order.tiff')
four_tile_plot('elpd', 'precip', mod_perf_df,
               c(-20,0), expression(Delta*" LPPD"), 
               'results/lppd_diff_precip_2022.2_order.tiff')

# # differences in absolute lppd plots
# four_tile_plot(form_diff_lppd_df, 'elpd', 'airt', mod_no_r_df,
#                c(-60,0),  expression(Delta*" LPPD"), 
#                'results/lppd_diff_airt_NO_RIDGE.tiff')
# four_tile_plot(form_diff_lppd_df, 'elpd', 'precip', mod_no_r_df,
#                c(-60,0), expression(Delta*" LPPD"), 
#                'results/lppd_diff_precip_NO_RIDGE.tiff')


# best model tables/tests ---------------------------------

# best models tables
best_mod_tab <- function( resp, clim_v ){
  
  form_diff_lppd_df( resp, clim_v, mod_perf_df) %>% 
    subset( mod_rank == 1) %>% 
    mutate( model = gsub('[0-9]','',model) ) %>% 
    count( model ) %>% 
    arrange( desc(n) ) %>% 
    mutate( resp = resp )
  
}

# put its ALL together
best_mods_out <- list( 
      best_mod_tab( 'surv', 'precip' ),
      best_mod_tab( 'grow', 'precip' ),
      best_mod_tab( 'fec', 'precip' ),
      best_mod_tab( 'log_lambda', 'precip' ),
      best_mod_tab( 'surv', 'airt' ),
      best_mod_tab( 'grow', 'airt' ),
      best_mod_tab( 'fec', 'airt' ),
      best_mod_tab( 'log_lambda', 'airt' )
      ) %>% 
  bind_rows %>% 
  group_by( model ) %>% 
  summarise( n = sum(n) ) %>% 
  ungroup %>% 
  arrange( desc(n) ) %>% 
  mutate( prop_best = n / sum(n) )

# write it
write.csv( best_mods_out, 'results/best_mods_prop.csv',
           row.names=F )


# Vital rates breakdown 
best_mods_out <- list( 
  best_mod_tab( 'surv', 'precip' ),
  best_mod_tab( 'grow', 'precip' ),
  best_mod_tab( 'fec', 'precip' ),
  best_mod_tab( 'log_lambda', 'precip' ),
  best_mod_tab( 'surv', 'airt' ),
  best_mod_tab( 'grow', 'airt' ),
  best_mod_tab( 'fec', 'airt' ),
  best_mod_tab( 'log_lambda', 'airt' ) ) %>% 
  bind_rows %>% 
  group_by( model, resp ) %>% 
  summarise( n = sum(n) ) %>% 
  ungroup %>% 
  mutate( prop_best = n / sum(n) ) %>% 
  arrange( resp, desc(prop_best) ) %>% 
  write.csv(  'results/best_vr_mods_prop.csv',
              row.names=F )


# best models tables
best_mod_prop <- function( resp, clim_v ){
  
  form_diff_lppd_df( resp, clim_v, mod_perf_df) %>% 
    subset( mod_rank == 1) %>% 
    mutate( model = gsub('[0-9]','',model) ) %>% 
    count( model ) %>% 
    arrange( desc(n) ) %>% 
    mutate( prop = n / sum(n),
            type = paste0(resp,'_',clim_v) )
  
}

# put its ALL together
mods_props <- list( 
  best_mod_prop( 'surv', 'precip' ),
  best_mod_prop( 'grow', 'precip' ),
  best_mod_prop( 'fec', 'precip' ),
  best_mod_prop( 'log_lambda', 'precip' ),
  best_mod_prop( 'surv', 'airt' ),
  best_mod_prop( 'grow', 'airt' ),
  best_mod_prop( 'fec', 'airt' ),
  best_mod_prop( 'log_lambda', 'airt' ) ) %>% 
  bind_rows %>% 
  mutate( prop_tras = asin(sqrt(prop)) ) 
  
# Tukey test on proportion of models selected
mods_props %>% 
  subset( !(model %in% c('ctrl','yr')) ) %>% 
  aov( prop_tras ~ model , data=.) %>% 
  TukeyHSD %>% 
  .$model %>% 
  write.csv( 'results/best_mods_tukey.csv' )

mods_props %>% 
  subset( !(model %in% c('ctrl','yr')) ) %>% 
  boxplot( prop ~ model, data=.)

mods_props %>% 
  subset( !(model %in% c('ctrl','yr')) ) %>% 
  group_by( model ) %>% 
  summarise( mean_p   = mean(prop),
             median_p = median(prop) )


# Continuous tests ---------------------------------------

# format the differences in lppd
delta_lppd_pos_df <- function(resp, clim_v, mod_df){
  
  perf_by_spp <- mod_df %>%
                  subset( response == resp & clim_var == clim_v ) %>%
                  dplyr::select(species, model, rep_yr, rep_n, elpd) %>% 
                  mutate( elpd = replace(elpd, elpd == -Inf, NA)) %>%
                  spread(model, elpd )
  
  perf_mat    <- perf_by_spp %>% 
                    dplyr::select(-species, -rep_yr, -rep_n) %>% t
  
  # get minimum of maximum values
  get_max     <- function(ii) which(perf_mat[,ii] == max(perf_mat[,ii],na.rm=T))
  
  # species sequence
  spp_seq     <- 1:ncol(perf_mat)
  
  # matrix of benchmark ids
  max_ids     <- sapply(spp_seq, get_max) %>% cbind( spp_seq )
  
  # benchmark values
  max_val     <- perf_mat[max_ids]
  
  # differences matrix
  diff_mat    <- sweep(perf_mat, 2,  max_val, FUN='-') 
  diff_mat    <- -diff_mat
  
  # differences data frame
  long_df <- diff_mat %>%
    as.data.frame %>%
    setNames( perf_by_spp$species ) %>%
    tibble::add_column(model = row.names(.), .before=1) %>%
    gather(species,elpd, setdiff(names(.),'model') ) %>% 
    full_join( spp_df ) %>% 
    mutate( model = factor(model, levels = mod_labs) ) %>%
    left_join( dplyr::select(perf_by_spp, species, rep_yr, rep_n) ) %>% 
    arrange( rep_yr, rep_n, species, model )  %>% 
    mutate( species = replace(species,
                              grepl('Eriogonum',species),
                              'Eriogonum_longifolium...') ) %>% 
    mutate( species = factor(species, levels = unique(species)) )
  
  long_df %>% 
    mutate( resp     = resp,
            clim_var = clim_v )
  
}

# lppd weight make
lppd_weight_make <- function(resp, clim_v, mod_df){
  
  perf_by_spp <- mod_df %>%
                  subset( response == resp & clim_var == clim_v ) %>%
                  dplyr::select(species, model, rep_yr, rep_n, elpd) %>% 
                  mutate( elpd = replace(elpd, elpd == -Inf, NA)) %>%
                  spread(model, elpd ) %>% 
                  mutate( species = gsub('_',' ',species) )
  
  spp_df      <- spp_df %>% 
                  mutate( species = gsub('_',' ',species) )
  
  perf_mat    <- perf_by_spp %>% 
    dplyr::select(-species, -rep_yr, -rep_n) %>% t
  
  # get minimum of maximum values
  get_max     <- function(ii) which(perf_mat[,ii] == max(perf_mat[,ii],na.rm=T))
  
  # species sequence
  spp_seq     <- 1:ncol(perf_mat)
  
  # matrix of benchmark ids
  max_ids     <- sapply(spp_seq, get_max) %>% cbind( spp_seq )
  
  # benchmark values
  max_val     <- perf_mat[max_ids]
  
  # differences matrix
  diff_mat    <- sweep( perf_mat, 2,  max_val, FUN='-' )
  weight_exp  <- function(x) exp(0.5*x)
  exp_mat     <- apply( diff_mat, 2, weight_exp )
  sum_mat     <- apply( exp_mat, 2, sum, na.rm=T )
  weight_mat  <- sweep( exp_mat, 2, sum_mat, FUN = '/' )
  
  # differences data frame
  long_df <- weight_mat %>%
    as.data.frame %>%
    setNames( perf_by_spp$species ) %>%
    tibble::add_column(model = row.names(.), .before=1) %>%
    gather(species, lppd_weight, setdiff(names(.),'model') ) %>% 
    full_join( spp_df ) %>% 
    mutate( model = factor(model, levels = new_labs) ) %>%
    # left_join( dplyr::select(perf_by_spp, species, rep_yr, rep_n) ) %>% 
    arrange( rep_yr, rep_n, species, model )  %>% 
    mutate( species = replace(species,
                              grepl('Eriogonum',species),
                              'Eriogonum_longifolium...') ) %>% 
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
    mutate( species = factor(species, levels = unique(species)) )
  
  spp_v  <- long_df$species %>% unique %>% as.character
  
  # model rank
  model_rank <- function(x){
    long_df %>%
      subset( species == x) %>%
      arrange( desc(lppd_weight) ) %>% 
      mutate( mod_rank = c(1:nrow(.)) ) %>% 
      mutate( mod_rank = replace(mod_rank,
                                 mod_rank > 3,
                                 NA) ) %>% 
      # mutate( mod_rank = mod_rank + 20 ) %>%
      mutate( mod_rank = as.character(mod_rank) ) %>% 
      select( model, species, mod_rank ) 
      
  }
  
  # pick the best models
  mod_sel_df <- lapply(spp_v, model_rank) %>% 
    bind_rows %>% 
    right_join( long_df )
  
  # add response variable and climate variable
  mod_sel_df %>% 
    mutate( resp     = resp,
            clim_var = clim_v )
  
}

# single data frame with positive delta lppd logged
pos_delta_lppd <- list( delta_lppd_pos_df( 'surv',       'precip', mod_perf_df ),
                        delta_lppd_pos_df( 'grow',       'precip', mod_perf_df ),
                        delta_lppd_pos_df( 'fec',        'precip', mod_perf_df ),
                        delta_lppd_pos_df( 'log_lambda', 'precip', mod_perf_df ),
                        delta_lppd_pos_df( 'surv',       'airt',   mod_perf_df ),
                        delta_lppd_pos_df( 'grow',       'airt',   mod_perf_df ),
                        delta_lppd_pos_df( 'fec',        'airt',   mod_perf_df ),
                        delta_lppd_pos_df( 'log_lambda', 'airt',   mod_perf_df ) ) %>% 
              bind_rows %>% 
  mutate( model    = as.character(model) ) %>% 
  mutate( model    = gsub('[0-9]','',model) ) %>% 
  mutate( model    = as.factor(model) ) %>% 
  mutate( elpd     = replace(elpd, elpd == 0, 0.00003554812) ) %>% 
  mutate( log_elpd = log(elpd) )
        

# single data frame with positive delta lppd logged
response_df <- data.frame( resp     = c('surv', 'grow', 'fec', 'log_lambda'),
                           Response = c('Survival',     'Development', 
                                        'Reproduction', 'Log(\u03BB)') )

# Format the LPPD weights data frame
lppd_weight_df <- list( lppd_weight_make( 'surv',       'precip', mod_perf_df ),
                        lppd_weight_make( 'grow',       'precip', mod_perf_df ),
                        lppd_weight_make( 'fec',        'precip', mod_perf_df ),
                        lppd_weight_make( 'log_lambda', 'precip', mod_perf_df ),
                        lppd_weight_make( 'surv',       'airt',   mod_perf_df ),
                        lppd_weight_make( 'grow',       'airt',   mod_perf_df ),
                        lppd_weight_make( 'fec',        'airt',   mod_perf_df ),
                        lppd_weight_make( 'log_lambda', 'airt',   mod_perf_df ) ) %>% 
  bind_rows %>% 
  mutate( lppd_weight_asin = asin(sqrt(lppd_weight)) ) %>% 
  mutate( model_type = 'mw' ) %>% 
  mutate( model_type = replace( model_type, grepl('yr',model), 'yr') ) %>% 
  mutate( model_type = replace( model_type, grepl('ctrl',model), 'ctrl') ) %>%
  mutate( model_type = as.factor(model_type) ) %>% 
  mutate( mv_type    = as.character(model) ) %>% 
  mutate( mv_type    = gsub(' [0-9]','',model) ) %>% 
  mutate( mv_type    = factor(mv_type, levels = c( 'NM',
                                                   'CSM',
                                                   'WMM',
                                                   'SAM',
                                                   'FHM') )) %>% 
  mutate( resp       = factor(resp,  levels = c( 'surv',
                                                 'grow',
                                                 'fec',
                                                 'log_lambda') )) %>% 
  mutate( mod_clim   = as.character(model) ) %>% 
  mutate( mod_clim   = replace( mod_clim, mod_clim != 'ctrl1', 'climate') ) %>% 
  mutate( mod_clim   = factor(mod_clim, levels = c('ctrl1','climate')) ) %>% 
  mutate( mod_fam    = gsub(' [0-9]{1}','',model) ) %>% 
  mutate( mod_fam    = factor(mod_fam, levels=c('NM','CSM','WMM','SAM','FHM')) ) %>% 
  mutate( year       = gsub('[A-Z]{3} |NM','',model) ) %>% 
  left_join( response_df ) %>% 
  # Order output
  mutate( Response   = factor(Response, levels = c( 'Survival', 
                                                    'Development', 
                                                    'Reproduction', 
                                                    'Log(\u03BB)') ))  
      

# calculate average lppd_weight for precip/airt
plot_lppd_weight_clim_resp_mod <- function( clim_var_chr ){

  p1 <- lppd_weight_df %>% 
    group_by( model, resp, clim_var ) %>%
    summarise( mean_lppd_weight = mean(lppd_weight, na.rm=T) ) %>% 
    ungroup %>% 
    mutate( resp = factor(resp, levels=c('surv', 
                                         'grow',
                                         'fec' , 
                                         'log_lambda')) ) %>% 
    subset( clim_var == clim_var_chr ) %>% 
    arrange( resp, model ) %>% 
    ggplot() +
    geom_line( aes( x     = model,
                    y     = mean_lppd_weight,
                    group = resp),
               lwd = 1 ) +
    geom_point( aes( x     = model,
                     y     = mean_lppd_weight,
                     group = resp) ) +
    geom_vline( xintercept = c(1.5,4.5,7.5,10.5),
                col = 'red', lwd=0.5, lty=2) +
    facet_wrap( ~ resp, ncol = 4 ) +
    theme_bw() +
    ylab( 'Average model weight' ) +
    # ylim(0, 0.11) +
    scale_color_colorblind() +
    scale_x_discrete(position = "top") + 
    theme( strip.background = element_blank(),
           strip.text.x     = element_blank(),
           panel.spacing    = unit(1.35, "lines"),
           axis.text.x      = element_blank(),
           # axis.ticks.y   = element_blank(),
           axis.title     = element_blank() ) 
  
    ggsave( paste0('results/mod_weights_',clim_var_chr,'.tiff'),
            plot = p1, width=10, height = 2.5, compression = 'lzw' )
    
}

# Store plots
plot_lppd_weight_clim_resp_mod( 'precip' )
plot_lppd_weight_clim_resp_mod( 'airt' )

# model type X response variable
lppd_weight_df %>% 
  group_by( model, resp, clim_var ) %>%
  summarise( mean_lppd_weight = mean(lppd_weight, na.rm=T) ) %>% 
  ungroup %>% 
  # subset( clim_var == 'precip' ) %>% 
  ggplot() +
  geom_line( aes( x     = model,
                   y     = mean_lppd_weight,
                  group = resp,
                  color = resp),
              lwd = 1 ) +
  geom_point( aes( x     = model,
                  y     = mean_lppd_weight,
                  group = resp,
                  color = resp) ) +
  geom_vline( xintercept = c(1.5,4.5,7.5,10.5),
              col = 'red', lwd=0.5, lty=2) +
  theme_bw() +
  scale_color_colorblind() +
  facet_wrap( ~clim_var ) +
  labs( y = 'Average model weight',
        color = 'Response' ) +
  theme( axis.text.x = element_text(angle = 90,
                                    hjust = 1,
                                    vjust = 0.5) ) +
  ggsave( 'results/weights_by_mod_resp_clim.tiff',
          width = 6.3, height = 3.15, 
          compression = 'lzw')

lppd_weight_df %>% 
  ggplot() +
  geom_boxplot( aes( x = mv_type,
                     y = lppd_weight ),
                outlier.shape = NA) +
  ylim(0,0.15) +
  labs( x = 'Model type',
        y = 'Model weight' ) +
  theme_minimal() +
  theme( axis.text.x = element_text(angle = 90,
                                    hjust = 1,
                                    vjust = 0.5,
                                    size  = 15),
         axis.title  = element_text( size = 20 ) ) +
  ggsave( 'results/weights_by_mod_boxplot.tiff',
          width = 6.3, height = 6.3, 
          compression = 'lzw')

# model type X response variable
lppd_weight_df %>% 
  group_by( mv_type ) %>%
  summarise( mean_lppd_weight = mean(lppd_weight, na.rm=T) ) %>% 
  ungroup %>% 
  ggplot() +
  geom_line( aes( x = mv_type,
                  y = mean_lppd_weight,
                  group = 1),
             lwd = 2 ) +
  geom_point( aes( x     = mv_type,
                   y     = mean_lppd_weight), size = 5 ) +
  theme_bw() +
  scale_color_colorblind() +
  theme( axis.text.x = element_text(angle = 90,
                                    hjust = 1,
                                    vjust = 0.5,
                                    size  = 15),
         axis.title  = element_text( size = 20 ) ) +
  labs( x = 'Model type',
        y = 'Average model weight' ) +
  ylim( 0, 0.10 ) +
  ggsave( 'results/weights_by_mod_line.tiff',
          width = 6.3, height = 6.3, 
          compression = 'lzw')

# Boxplot
p_box <- lppd_weight_df %>% 
  ggplot() +
  geom_jitter( aes( x = mv_type,
                    y = lppd_weight_asin,
                    group = Response,
                    color = Response),
               alpha = 0.6,
               width = 0.3 ) +
  geom_boxplot( aes( x = mv_type, 
                     y = lppd_weight_asin ),
                outlier.shape = NA,
                lwd = 1,
                alpha = 0.5) +
  theme_few() +
  scale_color_colorblind() +
  labs( x = 'Model',
        y = 'arcsin(LPPD weight)',
        color = 'Response') 

ggsave( 'results/weight_by_mod_box_dot.tiff',
        p_box,
        width = 6.3, height = 4, 
        compression = 'lzw')
 


# Boxplot
p_box_prec <- lppd_weight_df %>% 
  subset( clim_var == 'precip' ) %>% 
  ggplot() +
  geom_jitter( aes( x = mv_type,
                    y = lppd_weight_asin,
                    group = Response,
                    color = Response),
               size = 0.5,
               alpha = 0.6,
               width = 0.3 ) +
  geom_boxplot( aes( x = mv_type, 
                     y = lppd_weight_asin ),
                outlier.shape = NA,
                lwd = 0.5,
                alpha = 0.5) +
  theme_few() +
  scale_color_colorblind() +
  ylim( 0, 0.75 ) +
  labs( x = 'Model',
        y = 'arcsin(LPPD weight)',
        color = 'Response',
        title = 'Precipitation' )+
  theme( legend.position = 'none',
         plot.title  = element_text( hjust = 0.5,
                                     size  = 12 ),
         title       = element_text( size  = 8 ),
         axis.text   = element_text( size  = 8 ),
         legend.text = element_text( size  = 8 ) ) 

p_box_airt <- lppd_weight_df %>% 
  subset( clim_var == 'airt' ) %>% 
  ggplot() +
  geom_jitter( aes( x = mv_type,
                    y = lppd_weight_asin,
                    group = Response,
                    color = Response),
               size = 0.5,
               alpha = 0.6,
               width = 0.3 ) +
  geom_boxplot( aes( x = mv_type, 
                     y = lppd_weight_asin ),
                outlier.shape = NA,
                lwd = 0.5,
                alpha = 0.5) +
  theme_few() +
  scale_color_colorblind() +
  ylim( 0, 0.75 ) +
  labs( x = 'Model',
        y = 'arcsin(LPPD weight)',
        color = 'Response',
        title = 'Air temperature' ) +
  theme( plot.title  = element_text( hjust = 0.5,
                                     size  = 12 ),
         title       = element_text( size  = 8 ),
         axis.text   = element_text( size  = 8 ),
         legend.text = element_text( size  = 8 )
         ) 
  

p_box_airt_prec <- grid.arrange( p_box_prec, p_box_airt, 
                                 nrow = 1, widths = c(0.65,1)  )

ggsave( 'results/weight_by_mod_and_climvar_box_dot.tiff',
        p_box_airt_prec,
        width = 6.3, height = 3, 
        compression = 'lzw')


# Precipitation with the legend
p_box_prec <- lppd_weight_df %>% 
  subset( clim_var == 'precip' ) %>% 
  ggplot() +
  geom_jitter( aes( x = mv_type,
                    y = lppd_weight_asin,
                    group = Response,
                    color = Response),
               size = 2,
               alpha = 0.6,
               width = 0.3 ) +
  geom_boxplot( aes( x = mv_type, 
                     y = lppd_weight_asin ),
                outlier.shape = NA,
                lwd = 0.5,
                alpha = 0.5) +
  theme_few() +
  scale_color_colorblind() +
  ylim( 0, 0.6 ) +
  labs( x = 'Model',
        y = 'arcsin(LPPD weight)',
        color = 'Response:',
        title = 'Precipitation' )+
  theme( plot.title  = element_blank(),
         title       = element_text( size  = 15 ),
         axis.text   = element_text( size  = 15 ),
         legend.text = element_text( size  = 15 ),
         legend.position = 'bottom' ) 

ggsave( 'results/weight_by_mod_and_precip_box_dot_BES2022.tiff',
        p_box_prec,
        width = 8, height = 5, 
        compression = 'lzw')



# Split by response variable and model  
mod_resp_lab_df <-
  expand.grid( mod_fam  = c('NM', 'CSM', 'WMM', 'SAM', 'FHM'),
               resp     = c('surv','grow','fec','log_lambda'),
               stringsAsFactors = F ) %>% 
  left_join( response_df ) %>% 
  mutate( resp  = factor(resp, levels = c('surv','grow','fec','log_lambda')) ) %>% 
  arrange( resp ) %>% 
  mutate( resp  = as.character(resp) ) %>% 
  mutate( order = 1:nrow(.) ) %>% 
  mutate( mod_resp = paste0(mod_fam,' ',Response) )
  # .$mod_resp_label

# Preliminary plot
p_box <- lppd_weight_df %>% 
  left_join( response_df ) %>% 
  subset( clim_var == 'airt' ) %>% 
  mutate( mod_resp = paste0( mod_fam,' ',Response ) ) %>% 
  mutate( Response = factor( Response, levels =  c('Survival','Development','Reproduction',
                                                   'Log(\u03BB)')) ) %>% 
  # format symbol
  mutate( mod_resp = factor( mod_resp, levels = mod_resp_lab_df$mod_resp) ) %>% 
  ggplot() +
  geom_jitter( aes( x = mod_resp,
                    y = lppd_weight_asin,
                    group = Response,
                    color = Response ),
               alpha = 0.75,
               width = 0.3 ) +
  geom_boxplot( aes( x = mod_resp, 
                     y = lppd_weight_asin ),
                outlier.shape = NA,
                lwd = 1,
                alpha = 0.5 ) +
  theme_few() +
  scale_color_colorblind() +
  theme( axis.text = element_text( angle=90 ) ) +
  labs( x = 'Model',
        y = 'arcsin(LPPD weight)',
        color = 'Response' )

ggsave( 'results/weight_by_mod_resp_box_dot.tiff',
        p_box, 
        width = 6.3, height = 4, 
        compression = 'lzw')

  

# Precipitation plot
p_box_prec <- lppd_weight_df %>% 
  left_join( response_df ) %>% 
  subset( clim_var == 'precip' ) %>% 
  mutate( mod_resp = paste0( mod_fam,' ',Response ) ) %>% 
  mutate( Response = factor( Response, levels =  c('Survival','Development','Reproduction',
                                                   'Log(\u03BB)')) ) %>% 
  # format symbol
  mutate( mod_resp = factor( mod_resp, levels = mod_resp_lab_df$mod_resp) ) %>% 
  ggplot() +
  geom_jitter( aes( x = mod_resp,
                    y = lppd_weight_asin,
                    group = Response,
                    color = Response ),
               size  = 0.5,
               alpha = 0.85,
               width = 0.3 ) +
  geom_boxplot( aes( x = mod_resp, 
                     y = lppd_weight_asin ),
                outlier.shape = NA,
                lwd   = 0.5,
                alpha = 0.5 ) +
  theme_few() +
  scale_color_colorblind() +
  ylim( 0, 0.75 ) +
  labs( x = 'Model',
        y = 'arcsin(LPPD weight)',
        color = 'Response',
        title = 'Precipitation' ) +
  theme( axis.text.x = element_blank(),
         plot.title  = element_text( hjust = 0.5 ) )

# Air temperature
p_box_airt <- lppd_weight_df %>% 
  left_join( response_df ) %>% 
  subset( clim_var == 'airt' ) %>% 
  mutate( mod_resp = paste0( mod_fam,' ', Response ) ) %>% 
  mutate( Response = factor( Response, levels =  c('Survival','Development','Reproduction',
                                                   'Log(\u03BB)')) ) %>% 
  # format symbol
  mutate( mod_resp = factor( mod_resp, levels = mod_resp_lab_df$mod_resp) ) %>% 
  ggplot() +
  geom_jitter( aes( x = mod_resp,
                    y = lppd_weight_asin,
                    group = Response,
                    color = Response ),
               size  = 0.5,
               alpha = 0.85,
               width = 0.3 ) +
  geom_boxplot( aes( x = mod_resp, 
                     y = lppd_weight_asin ),
                outlier.shape = NA,
                lwd   = 0.5,
                alpha = 0.5 ) +
  theme_few() +
  scale_color_colorblind() +
  ylim( 0, 0.75 ) +
  labs( x = 'Model',
        y = 'arcsin(LPPD weight)',
        color = 'Response',
        title = 'Air temperature' ) +
  theme( axis.text.x = element_text( angle = 90,
                                     vjust = 0.5,
                                     hjust = 1,
                                     margin = margin( l = 0)),
         legend.text = element_text( size = 20),
         plot.title = element_text( hjust = 0.5 ) )


p_out <- ggarrange( p_box_prec, p_box_airt, heights = c(0.75,1),
                    ncol=1, nrow=2, common.legend = TRUE, legend="right")

ggsave( 'results/weight_by_mod_resp_box_dot_prec-temp.tiff',
        p_out, 
        width = 6.3, height = 8, 
        compression = 'lzw')




lppd_weight_df %>% 
  mutate( mod_resp = paste0( mod_fam,'_',resp ) ) %>% 
  # left_join( mod_resp_lab_df ) %>% 
  mutate( mod_resp = factor(mod_resp, levels = mod_resp_lab_df$mod_resp) ) %>%
  .$mod_resp  
    

cia <- c('NM', 'CSM', 'WMM', 'SAM', 'FHM') %>% sort
factor( cia, levels = c('NM', 'CSM', 'WMM', 'SAM', 'FHM') ) %>% as.numeric


# model weight by RESPONSE
lppd_weight_df %>% 
  subset( mv_type == 'NM' ) %>% 
  ggplot() +
  geom_jitter( aes( x = Response,
                    y = lppd_weight_asin),
               alpha = 0.75,
               width = 0.3 ) +
  geom_boxplot( aes( x = Response, 
                     y = lppd_weight_asin ),
                outlier.shape = NA,
                lwd = 1,
                alpha = 0.5) +
  theme_few() +
  theme( axis.text  = element_text( size = 15),
         axis.title = element_text( size = 20) ) +
  scale_color_colorblind() +
  labs( x = 'Model',
        y = 'arcsin(LPPD weight)',
        color = 'Response') +
  ggsave( 'results/weight_nullmodel_by_response_dot.tiff',
          width = 6.3, height = 6.3, 
          compression = 'lzw')


# model weight by RESPONSE
p_noNM_box <- lppd_weight_df %>% 
  subset( mv_type == 'NM' ) %>% 
  ggplot() +
  geom_jitter( aes( x = Response,
                    y = lppd_weight_asin),
               size  = 2,
               alpha = 0.75,
               width = 0.3 ) +
  geom_boxplot( aes( x = Response, 
                     y = lppd_weight_asin ),
                outlier.shape = NA,
                lwd = 1,
                alpha = 0.5) +
  theme_few() +
  theme( axis.text  = element_text( size = 15),
         axis.title = element_text( size = 20) ) +
  scale_color_colorblind() +
  ylim( 0, 0.6 ) +
  labs( x = 'Model',
        y = 'arcsin(LPPD weight)',
        color = 'Response',
        title = 'Predictive ability of null models (no climate)' ) +
  theme( plot.title = element_text( hjust = 0.5,
                                    size  = 20 ) )

ggsave( 'results/weight_nullmodel_by_response_dot_BES2022.tiff',
        p_noNM_box,
        width = 6.3, height = 6.3, 
        compression = 'lzw')


# model weight by rep and model type
lppd_weight_df %>% 
  group_by( species, rep_yr, mv_type ) %>%
  summarise( mean_lppd_weight = mean(lppd_weight, na.rm=T) ) %>% 
  ungroup %>% 
  ggplot() +
  geom_point( aes( x    = rep_yr,
                   y     = mean_lppd_weight) ) +
  theme_bw() +
  scale_color_colorblind() +
  labs( x = 'Replication (years)',
        y = 'Model weight' ) +
  ylim(0, 0.15) +
  facet_wrap( ~mv_type, ncol = 2 ) +
  theme( axis.title = element_text( size = 25),
         axis.text  = element_text( size = 20),
         strip.text = element_text( size = 20) ) +
  ggsave( 'results/weight_by_rep.tiff',
          width = 6.3, height = 9.9, 
          compression = 'lzw')


# model weight by rep and model type  
lppd_weight_df %>% 
  # left_join( nm_weight ) %>% 
  mutate( delta_w = lppd_weight - nm_weight ) %>% 
  subset( !(model == 'NM') ) %>% 
  group_by( species, resp, mv_type, clim_var, rep_yr ) %>% 
  summarise( lppd_weight = max(lppd_weight) ) %>% 
  ungroup %>% 
  right_join( expand.grid( rep_yr = c(6:35),
                           mv_type  = c('CSM','WMM','SAM','FHM') ) ) %>% 
  subset( !is.na(mv_type) ) %>% 
  ggplot() +
  geom_boxplot( aes( x     = as.factor(rep_yr),
                     y     = lppd_weight ) ) +
  theme_bw() +
  scale_color_colorblind() +
  labs( x = 'Replication (years)',
        y = 'Model weight' ) +
  ylim(0, 0.15) +
  facet_wrap( ~mv_type, ncol = 2 ) +
  theme( axis.title = element_text( size = 20),
         axis.text  = element_text( size = 10),
         strip.text = element_text( size = 15) ) +
  scale_x_discrete( breaks =  c(6,12,18,24,30,35) ) +
  ggsave( 'results/weight_by_rep_boxplot.tiff',
          width = 6.3, height = 6.3, 
          compression = 'lzw' )
  
# QUANTILE REGRESSIONS ---------------------------------------------------------

# Separate files for quantile regression
csm_df      <- lppd_weight_df %>% subset( mv_type == 'CSM' ) 
wmm_df      <- lppd_weight_df %>% subset( mv_type == 'WMM' ) 
sam_df      <- lppd_weight_df %>% subset( mv_type == 'SAM' ) 
fhm_df      <- lppd_weight_df %>% subset( mv_type == 'FHM' ) 

# fit quantile regressions
csm_mod     <- rq(lppd_weight_asin ~ rep_yr, tau = seq(0.1,0.9,by=0.1), data = csm_df ) 
wmm_mod     <- rq(lppd_weight_asin ~ rep_yr, tau = seq(0.1,0.9,by=0.1), data = wmm_df ) 
sam_mod     <- rq(lppd_weight_asin ~ rep_yr, tau = seq(0.1,0.9,by=0.1), data = sam_df ) 
fhm_mod     <- rq(lppd_weight_asin ~ rep_yr, tau = seq(0.1,0.9,by=0.1), data = fhm_df ) 

# hardcode viridis scale
colors      <- viridis_pal(begin = 0.1, end = 0.9)(9)

# self explanatory
exrtact_tau_estimates <- function( mod_name ){
  
  if( mod_name == 'CSM' ) mod_df = csm_mod
  if( mod_name == 'WMM' ) mod_df = wmm_mod
  if( mod_name == 'SAM' ) mod_df = sam_mod
  if( mod_name == 'FHM' ) mod_df = fhm_mod
  
  mod_df$coefficients %>% 
    as.data.frame %>% 
    mutate( coef = c('b0','b1') ) %>% 
    pivot_longer( 1:9, names_to = 'Quantile', values_to = 'b0b1' ) %>% 
    mutate( Quantile = gsub('tau= ','',Quantile) ) %>% 
    pivot_wider( names_from  = 'coef',
                 values_from = 'b0b1' ) %>% 
    mutate( mv_type = mod_name )
  
}

# Plot quantile regressions against data
plot_quantile_bivariate <- function( mod_name ){
  
  lppd_weight_df %>% 
    subset( mv_type == mod_name ) %>% 
    ggplot() +
    geom_point( aes( x  = rep_yr,
                     y  = lppd_weight_asin) ) +
    geom_abline( data = exrtact_tau_estimates( mod_name ),
                 aes( intercept=b0, 
                      slope = b1,
                      group = Quantile,
                      color = Quantile),
                 lwd = 2) +
    scale_color_manual( values = colors,
                        labels = exrtact_tau_estimates( mod_name )$Quantile ) +
    theme_bw() +
    labs( x = 'Replication (years)',
          y = 'Arcsin(LPPD weight)',
          title = mod_name ) +
    theme( axis.title = element_text( size  = 13),
           axis.text  = element_text( size  = 10),
           strip.text = element_text( size  = 10),
           plot.title = element_text( hjust = 0.5) )
  
}

# Put it all together
csm_p <- plot_quantile_bivariate( 'CSM' )
wmm_p <- plot_quantile_bivariate( 'WMM' )
sam_p <- plot_quantile_bivariate( 'SAM' )
fhm_p <- plot_quantile_bivariate( 'FHM' )
p_out <- ggarrange( wmm_p, sam_p, fhm_p,
                    ncol=1, nrow=3, common.legend = TRUE, legend="right")

ggsave( 'results/weight_by_rep_quantile.tiff',
        p_out,
        width = 4, height = 9, 
        compression = 'lzw')




  
# lppd weight by generation time
gen_times <- read.csv('C:/CODE/plant_generation_time_and_climate/data/generation_times.csv') %>% 
               # rename to allow the join down the line
               rename( species = SpeciesAuthor ) %>% 
               mutate( species = gsub('_',' ',species) ) %>% 
               mutate( species = gsub('Eriogonum longifolium var. gnaphalifolium 2',
                                      'Eriogonum_longifolium...',
                                      species) ) %>% 
               mutate( species = replace(species,
                                         grepl('Cirsium pitcheri 8',species),
                                         'Cirsium pitcheri (1)') ) %>% 
               mutate( species = replace(species,
                                         grepl('Cirsium pitcheri 4',species),
                                         'Cirsium pitcheri (2)') ) %>% 
               mutate( species = replace(species,
                                         grepl('Cirsium pitcheri 6',species),
                                         'Cirsium pitcheri (3)') )

# test overlap (35 over 38 species)
intersect( unique(lppd_weight_df$species),
           unique(gen_times$species) )


-# model weight by rep and model type
lppd_weight_df %>% 
  group_by( species, rep_yr, mv_type ) %>%
  summarise( mean_lppd_weight = mean(lppd_weight, na.rm=T) ) %>% 
  ungroup %>% 
  left_join( gen_times ) %>% 
  ggplot() +
  geom_point( aes( x     = log_gen_t,
                   y     = mean_lppd_weight) ) +
  theme_bw() +
  scale_color_colorblind() +
  labs( x = 'log(Generation time)',
        y = 'Model weight' ) +
  facet_wrap( ~mv_type, ncol = 2 ) +
  theme( axis.title = element_text( size = 25),
         axis.text  = element_text( size = 20),
         strip.text = element_text( size = 20) ) +
  ggsave( 'results/weight_by_genT.tiff',
          width = 6.3, height = 9.9, 
          compression = 'lzw')



# Hypothesis 1: do moving window models perform better than controls?
h1_aov <- aov( lppd_weight_asin ~ model_type, 
               data = lppd_weight_df)
  
# Tukey test
TukeyHSD(h1_aov) %>% 
  .$model_type %>% 
  write.csv( 'results/mod_type_comparison_tukey.csv' )

# Hypothesis 2: compare among moving window models
aov( lppd_weight_asin ~ mv_type, 
                 data = subset(lppd_weight_df, !(model_type == 'ctrl' | model_type == 'yr') ) 
                 ) %>% 
  TukeyHSD %>% 
  .$mv_type %>% 
  write.csv( 'results/movwin_type_comparison_tukey.csv' )

# Hypothesis 3: effect of years of replication 
h3_rep_yr <- lm( lppd_weight_asin ~ rep_yr, 
                 data = subset(lppd_weight_df, 
                               !(model_type == 'ctrl' | model_type == 'yr') ) 
                )

# store results
h3_rep_yr %>% 
  summary %>% 
  .$coefficients %>% 
  write.csv( 'results/rep_yr_and_predperf_mw.csv' )

# Hypothesis 4: are vital rates more predictable?
lm( lppd_weight_asin ~ mod_clim + mod_clim:resp, 
    data = lppd_weight_df) %>% 
  summary %>% 
  .$coefficients %>% 
  write.csv( 'results/vital_rates_predperf.csv' )



# black out models that did not converge --------------------------------
mod_diag <- function(ii){
  
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

# load model diagnostics
diag_df     <- lapply(1:nrow(input_df), mod_diag) %>% 
                  bind_rows %>% 
                  # Add "flags" for convergence issues
                  mutate( diverg  = n_diverg > 0,
                          rhat    = rhat_high > 0 ) %>% 
                  # select only variables useful for merging
                  select( species, model, response, clim_var,
                          diverg, rhat )

# NAs where one or more rhats > 1.1
mod_rhat_df <- diag_df %>% 
                  right_join( mod_perf_df) %>% 
                  # black out elpd if one rhat > 1.1
                  mutate( elpd = replace( elpd,
                                          rhat,
                                          NA) )

# NAs where one or more divergent transition
mod_dive_df <- diag_df %>% 
                  right_join( mod_perf_df) %>% 
                  # black out elpd if one rhat > 1.1
                  mutate( elpd = replace( elpd,
                                          diverg,
                                          NA) )
# differences in absolute lppd plots
four_tile_plot(form_diff_lppd_df, 'elpd', 'airt',   mod_rhat_df,
               c(-30,0), expression(Delta*" LPPD"), 
               'results/lppd_diff_rhat_airt.tiff')
four_tile_plot(form_diff_lppd_df, 'elpd', 'precip', mod_dive_df,
               c(-30,0), expression(Delta*" LPPD"), 
               'results/lppd_diff_div_precip.tiff')

# final 8-panel tile plots --------------------------------------

# 8-panel tile plot 
lppdW_tile_plot <- function( clim_var ){
  
  # tile plot function
  upper_tile_plot <- function( df, title_lab ){
    
    ggplot( df, aes(model, species) ) +
      geom_tile( aes(fill = lppd_weight), color = "white") +
      geom_point(aes(shape = mod_rank),
                 size   = 2,
                 colour = "grey" ) +
      geom_vline( xintercept = c(1.5,4.5,7.5,10.5),
                  col = 'red', lwd=1, lty=2) + 
      ggtitle( title_lab ) +
      theme( title        = element_text(angle = 0, hjust = 0.5, size = 10),
             axis.title   = element_blank(),
             axis.title.x = element_blank(),
             # axis.text.x  = element_text(angle = 90, 
             #                             hjust = 1,
             #                             vjust = 0.5),
             axis.text.x  = element_blank() )
    
  }
  
  # single tile-plots
  p1 <- upper_tile_plot( lppd_weight_make('surv', clim_var, mod_perf_df),
                         'A) Survival' ) +
    scale_fill_viridis(guide = F, limits = c(0,0.82) ) +
    theme(legend.title = element_text(size = 10),
          legend.text  = element_text(size = 12),
          legend.position ="none",
          legend.key   = element_blank()) +
    labs(fill = element_blank() )
  p2 <- upper_tile_plot( lppd_weight_make('grow', clim_var, mod_perf_df),
                         'B) Growth' ) +
    scale_fill_viridis(guide = F, limits = c(0,0.82) ) +
    theme(legend.title = element_blank(), 
          legend.position ="none",
          axis.text.y  = element_blank(),
          legend.text  = element_blank(),
          legend.key   = element_blank())
  p3 <- upper_tile_plot( lppd_weight_make('fec',  clim_var, mod_perf_df),
                         'C) Fecundity' ) +
    scale_fill_viridis(guide = F, limits = c(0,0.82) ) +
    theme(legend.title = element_blank(), 
          axis.text.y  = element_blank(),
          legend.text  = element_blank(),
          legend.position ="none",
          legend.key   = element_blank()) +
    labs(fill = element_blank() )
  p4 <- upper_tile_plot( lppd_weight_make('log_lambda', clim_var, mod_perf_df),
                         expression('D) Log('*lambda*')') ) +
    scale_fill_viridis( limits = c(0,0.82) ) +
    theme(legend.title = element_text(size = 10),
          legend.text  = element_text(size = 12),
          axis.text.y  = element_blank() ) +
    labs( fill = 'LPPDw',
          size = element_blank(),
          shape= 'Mod Rank')
  
  # bottom plot
  bottom_tile_plot <- function( resp, clim_var, mod_df ){           
    
    lppd_weight_make( resp, clim_var, mod_perf_df) %>% 
      group_by( model ) %>% 
      summarise( lppd_weight = mean(lppd_weight, na.rm=T) ) %>% 
      ungroup %>% 
      mutate( average = 'Average') %>% 
      ggplot( aes(model, average) ) +
      geom_tile( aes(fill = lppd_weight), color = "white") +
      geom_vline( xintercept = c(1.5,4.5,7.5,10.5),
                  col = 'red', lwd=1, lty=2) +
      theme_minimal() +
      labs( y    = '`               Average LPPDw',
            fill = 'LPPDw') +
      theme( title        = element_text(angle = 0, hjust = 0.5, size = 10),
             axis.text.y  = element_blank(),
             axis.title.y = element_text(angle = 0,
                                         vjust = 0.5 ),
             axis.title.x = element_blank(),
             axis.text.x  = element_text(angle = 90, 
                                         hjust = 1,
                                         vjust = 0.5) ) +
      scale_fill_viridis( guide = F ) 
    
  }
  
  p5 <- bottom_tile_plot( 'surv', clim_var, mod_perf_df)
  p6 <- bottom_tile_plot( 'grow', clim_var, mod_perf_df) +
    theme( axis.title.x = element_blank(),
           axis.title.y = element_blank() )
  p7 <- bottom_tile_plot( 'fec', clim_var, mod_perf_df) +
    theme( axis.title.x = element_blank(),
           axis.title.y = element_blank() )
  p8 <- bottom_tile_plot( 'log_lambda', clim_var, mod_perf_df) +
    theme( axis.title.x = element_blank(),
           axis.title.y = element_blank(),
           legend.title = element_text(size = 10),
           legend.text = element_text(size = 12) ) +
    scale_fill_viridis() 
  
  
  
  # plot 
  tiff( paste0('results/lppdW_',clim_var,'.tiff'), 
        unit="in", width=10, height=6.3, res=600,compression="lzw" )
  
  grid.arrange(nrow   = 2, ncol = 4, 
               grobs  = list(p1, p2, p3, p4,
                             p5, p6, p7, p8),
               widths  = c(1.8, 0.9, 0.9, 1.4),
               heights = c(2,0.5),
               layout_matrix = rbind(c(1, 2, 3, 4),
                                     c(5,6,7,8))
  )
  
  dev.off()

}

lppdW_tile_plot( 'precip' )
lppdW_tile_plot( 'airt' ) 
