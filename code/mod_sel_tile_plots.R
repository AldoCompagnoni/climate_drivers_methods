rm(list=ls())
library(dplyr)
library(ggplot2)
library(viridis)
library(gridExtra)

# TO DO: read new results from supercomputer runs.

# 1. Compute LOOIC, Delta LOOIC, Z Delta LOOIC, ELPD
# 2. Z Delta LOOIC plots
# 3. Delta LOOIC plot
# 4. ELPD plot 


# # 1. n_dev ------------------------------------------------------
# out_dir  <- 'C:/CODE/climate_drivers_methods/'
# clim_var <- 'precip'
# resp_l   <- c('log_lambda','surv','grow','fec')
# 
# # species list
# spp      <- read.csv("data/all_demog_updt.csv", 
#                    stringsAsFactors = F) %>% 
#               .$SpeciesAuthor %>% 
#               unique 
# 
# # model order no control
# mod_ordr <- c('ctrl1',   'yr1', 'yr2', 'yr3',
#               'gaus',    'expp',
#               'gaus_n',  'expp_n', 'simpl_n',
#               'movb_h',  'movb',
#               'movb_h_n','movb_n')
# 
# # read in files
# sum_f_l   <- grep('mod_summ',
#                   list.files(paste0('results/mod_sel/',
#                                     clim_var)),
#                   value=T) 
# 
# # extract function
# extract_element <- function(x,patt){
#   regmatches(x, 
#            gregexpr(patt, 
#            x) )
# }
#   
# # read all summaries
# f_paths   <- paste0('results/mod_sel/',clim_var,'/',sum_f_l)
# spp_v     <- gsub( paste0('mod_summaries_',resp_l,'_') %>% 
#                     paste(collapse='|') %>% 
#                     paste0('|.csv'),
#                    '',sum_f_l)
# resp_v    <- gsub( 'mod_summaries_', '', sum_f_l) %>%
#                 # delete ANYTHING less than 101 character long that follows:
#                 # either fec or surv or grow or log_lambda
#                 gsub(paste0('(?<=',paste(resp_l,collapse='|'),').{1,100}'), 
#                      '', ., perl=T)
# 
# 
# read_format_summ <- function(x,spp_x,resp_x){
#   read.csv(x) %>% 
#     tibble::add_column( species  = spp_x,  .before = 2) %>% 
#     tibble::add_column( response = resp_x, .before = 2)
# }
# 
# # all model summaries in one data frame
# mod_sum_df <- Map(read_format_summ, 
#                  f_paths, spp_v, resp_v) %>% 
#                 bind_rows 
# 
# mod_div_df <- dplyr::select(mod_sum_df, 
#                             model,species,response,n_diverg) %>% 
#                 mutate( model    = factor(model, levels = mod_ordr),
#                         response = factor( response, 
#                                            levels = c('surv',
#                                                       'grow',
#                                                       'fec',
#                                                       'log_lambda') 
#                                            )
#                         )

# # Store plot
# ggplot(mod_div_df, aes(x=model, y=n_diverg) ) +
#   geom_jitter( aes(color = response) ) +
#   # geom_jitter( aes(x=model) ) +
#   scale_color_viridis( discrete = T ) + 
#   theme( axis.text = element_text( angle = 90) ) +
#   ggsave('results/converg/div_trans_vs_mod.tiff',
#          width=6.3, height=5, compression='lzw')


# 1. Compute LOOIC, Delta LOOIC, Z Delta LOOIC, ELPD ------------------

# data frame containing inputs
input_df  <- expand.grid( resp     = c('log_lambda','surv',
                                       'grow','fec'),
                          clim_var = c('precip','airt'),
                          stringsAsFactors = F ) %>% 
                mutate( resp_clim = paste(resp, clim_var, sep='_') )

loo_l   <- list()

# Store looic and ELPDs
for(ii in 1:nrow(input_df) ){
  
  # out_dir  <- 'C:/CODE/climate_drivers_methods/'
  # clim_var <- 'precip'
  # response <- resp_l[ii]
  
  mod_dir   <- 'E:/work/sApropos/2019.6.11/'
  clim_var  <- input_df$clim_var[ii]
  response  <- input_df$resp_l[ii]
  resp_clim <- input_df$resp_clim[ii]
  
  # species list
  spp       <- read.csv("data/all_demog_updt.csv", 
                        stringsAsFactors = F) %>% 
                 .$SpeciesAuthor %>% 
                 unique
  
  # file list
  sum_f_l   <- grep('mod_summ',
                    list.files(mod_dir),
                    value=T) %>% 
                 grep(resp_clim, ., value=T)
  
  # path list
  path_l    <- paste0(mod_dir, sum_f_l)
  
  # extract species names
  spp_names <- gsub("mod_summaries_", "", sum_f_l ) %>%
                  gsub(paste0("_",resp_clim,".csv"), "", . ) %>% 
                  gsub('array_vr-[0-9]{7}-[0-9]{1,2}_', "", . )
  
  # download files serially
  create_sum_df <- function(spp_n, x){
    
    match_n <- sum( grepl(spp_n, x) )
    
    # proceed only if species name is in file name
    if( match_n == 1 ){
      file_p <- grep(spp_n, x, value=T)
      out_df <- read.csv(file_p, stringsAsFactors = F) %>% 
                  tibble::add_column(species = spp_n,.before=1)
      return( out_df )
    }
    if( match_n > 1 ){
      warning( 'Mistake: more than 1 spp. per clim_var/resp combination')
      return(NULL)
    }
    if( match_n == 0 ) return(NULL) 
    
  }
  
  # model order no control
  mod_ordr <- c('ctrl1',  'yr1','yr2','yr3',
                'gaus',   'expp',
                'gaus_n', 'expp_n', 'simpl_n',
                'mb_h',   'mb',
                'mb_h_n', 'mb_n')
  
  # elpd 
  plot_elpd <- function(sum_df, spp_v){
    
    # get model ranks (still for plotting!)
    model_rank <- function(x){
      sum_df %>%
        subset( species == x) %>%
        arrange( desc(elpd) ) %>% 
        select( model, species, elpd ) %>% 
        mutate( mod_rank_elpd = c(1:nrow(.)) ) %>% 
        mutate( mod_rank_elpd = replace(mod_rank_elpd,
                                        mod_rank_elpd > 3,
                                        NA) ) %>% 
        mutate( mod_rank_elpd = as.character(mod_rank_elpd) ) %>% 
        select( model, species, elpd, mod_rank_elpd )
    }
  
    lapply(spp_v, model_rank) %>% bind_rows
    
  }
  
  # rank the models based on looic
  plot_loo <- function(sum_df, spp_v){
    
    # calculate delta of performance (for plotting)
    delta_perf <- function(x){
      
      # performance of null model
      null_perf <- sum_df %>%
                    subset( species == x) %>% 
                    subset( model == 'ctrl1') %>% 
                    .$looic
      
      # spit out delta WAIC
      sum_df %>%
        subset( species == x) %>% 
        mutate( delta_looic   = looic - null_perf) %>% 
        mutate( delta_looic_z = delta_looic/se_looic )
        
    }
    
    delta_df <- lapply(spp_v, delta_perf) %>% bind_rows
    
    # get "significant" models
    model_sig <- function(x){
      delta_df %>%
        subset( species == x) %>%
        arrange( looic ) %>% 
        mutate( mod_sig = delta_looic_z < -2 ) %>% 
        mutate( mod_sig = as.numeric(mod_sig) %>% as.character) %>% 
        mutate( mod_sig = replace(mod_sig,
                                  mod_sig=='0',
                                  NA) ) %>% 
        select( model, species, looic, #se_looic,
                delta_looic, delta_looic_z, mod_sig )
    }
    
    mod_sig_df <- lapply(spp_v, model_sig) %>% bind_rows
    
    # get model ranks (still for plotting!)
    model_rank <- function(x){
      mod_sig_df %>%
        subset( species == x) %>%
        arrange( looic ) %>% 
        mutate( mod_rank = c(1:nrow(.)) ) %>% 
        mutate( mod_rank = replace(mod_rank,
                                   mod_rank > 3,
                                   NA) ) %>% 
        mutate( mod_rank = as.character(mod_rank) ) %>% 
        select( model, species, looic, #se_looic,
                delta_looic, delta_looic_z, 
                mod_sig, mod_rank )
    }
  
    lapply(spp_v, model_rank) %>% bind_rows
    
  }
  
  # data frame for tile plotss
  elpd_df <- lapply(spp_names,create_sum_df,path_l) %>% 
                bind_rows %>% 
                plot_elpd( spp_names ) 
    
  loo_df  <- lapply(spp_names,create_sum_df,path_l) %>% 
                bind_rows %>% 
                plot_loo( spp_names ) 
  
  tile_df <- full_join(  elpd_df, loo_df ) %>% 
                select(species,model,
                       looic,delta_looic_z, elpd,
                       mod_sig, mod_rank,
                       mod_rank_elpd) %>% 
                mutate( model = factor(model,
                                       levels = mod_ordr ) ) %>% 
                mutate( species    = substr(species,1,25) )
  
  loo_l[[ii]] <- tile_df
  
  # # store tile plots 
  # subset(tile_df, model != 'ctrl1' ) %>%
  #   ggplot( aes(model, species) ) +
  #   geom_tile(aes(fill = delta_looic_z), color = "white") +
  #   geom_point(aes(size  = '0.7',
  #                  shape = mod_sig) ) +
  #   scale_fill_viridis( limits = tile_df$delta_looic_z %>% range ) + 
  #   ggtitle( response ) +
  #   theme(title        = element_text(angle = 0, hjust = 0.5, size = 10),
  #         legend.title = element_text(size = 10),
  #         legend.text = element_text(size = 12),
  #         axis.title = element_blank(),
  #         axis.text.x = element_text(angle = 90, hjust = 1)) + 
  #   ggsave(paste0(out_dir,'results/looic/',
  #                 clim_var,'_',response,'_delta_looic_z.tiff'),
  #          width=6.3,height=6.3,compression='lzw')
  # 
  # # store tile plots 
  # ggplot(tile_df, aes(model, species)) +
  #   geom_tile(aes(fill = delta_looic_z), color = "white") +
  #   geom_point(aes(size  = '0.7',
  #                  shape = mod_rank) ) +
  #   scale_fill_viridis( limits = tile_df$delta_looic_z %>% range ) + 
  #   ggtitle( response ) +
  #   theme(title        = element_text(angle = 0, hjust = 0.5, size = 10),
  #         legend.title = element_text(size = 10),
  #         legend.text = element_text(size = 12),
  #         axis.title = element_blank(),
  #         axis.text.x = element_text(angle = 90, hjust = 1)) + 
  #   ggsave(paste0(out_dir,'results/looic/',
  #                 clim_var,'_',response,'_delta_looic.tiff'),
  #          width=6.3,height=6.3,compression='lzw')

}
  
# provide names to loo tile plots
loo_l <- setNames(loo_l, input_df$resp_clim)


# 2. Delta LOOIC plot -----------------------------------------

# plot delta_z
p_delta_z <- function(ii){
  
  response <- names(loo_l)[ii]

  tile_df <- loo_l[response][[1]] %>% 
    subset( model != 'ctrl1' ) %>%
    mutate( species = substr(species,1,25) ) %>% 
    rename( z = delta_looic_z)  
  
  # if log_lambda, print spp names
  if(ii == 1 | ii == 5 ){
    
    p_out <- tile_df %>% 
      ggplot( aes(model, species) ) +
      geom_tile(aes(fill = z), color = "white") +
      geom_point(aes(size  = '0.7',
                     shape = mod_sig),
                 show.legend = F) +
      scale_fill_viridis( limits = tile_df$z %>% range ) + 
      ggtitle( response ) +
      theme(title        = element_text(angle = 0, hjust = 0.5, size = 10),
            legend.title = element_text(size = 10),
            legend.text  = element_text(size = 12),
            legend.position = 'bottom',
            axis.title   = element_blank(),
            axis.text.x  = element_text(angle = 90, 
                                        hjust = 1,
                                        vjust = 0.5))
    
  }else{
    
    p_out <- tile_df %>% 
      ggplot( aes(model, species) ) +
      geom_tile(aes(fill = z), color = "white") +
      geom_point(aes(size  = '0.7',
                     shape = mod_sig),
                 show.legend = FALSE) +
      scale_fill_viridis( limits = tile_df$z %>% range ) + 
      ggtitle( response ) +
      theme(title        = element_text(angle = 0, hjust = 0.5, size = 10),
            legend.title = element_text(size = 10),
            legend.text  = element_text(size = 12),
            legend.position = 'bottom',
            axis.title   = element_blank(),
            axis.text.y = element_blank(),
            axis.text.x = element_text(angle = 90, 
                                        hjust = 1,
                                        vjust = 0.5))  
    
  }

}

# store precipitation results
p1 <- p_delta_z(1)
p2 <- p_delta_z(2)
p3 <- p_delta_z(3)
p4 <- p_delta_z(4)

# plot results 
tiff( 'results/mod_sel/looic/precip_delta_looic.tiff', 
      unit="in", width=10, 
      height=6.3, res=600,compression="lzw" )

grid.arrange(
     nrow = 1,
     grobs = list(p1, p2, p3, p4),
     widths = c(1.6, 0.9, 0.9, 0.9),
     layout_matrix = rbind(c(1, 2, 3, 4)) )

dev.off()


# store precipitation results
p1 <- p_delta_z(5)
p2 <- p_delta_z(6)
p3 <- p_delta_z(7)
p4 <- p_delta_z(8)

# store plot results 
tiff( 'results/mod_sel/looic/airt_delta_looic.tiff', 
      unit="in", width=10, 
      height=6.3, res=600,compression="lzw" )

grid.arrange(
     nrow = 1,
     grobs = list(p1, p2, p3, p4),
     widths = c(1.6, 0.9, 0.9, 0.9),
     layout_matrix = rbind(c(1, 2, 3, 4)) )

dev.off()


# 3. Delta LOOIC plot  ------------------------------------------

# plot loo plots
p_loo <- function(ii){
  
  response <- names(loo_l)[ii]

  tile_df <- loo_l[response][[1]] %>% 
    mutate( species = substr(species,1,25) ) %>% 
    rename( z = delta_looic_z )
  
  if(ii == 1 | ii == 5){
  
    p_out <- tile_df %>% 
      ggplot( aes(model, species) ) +
      geom_tile(aes(fill = z), color = "white") +
      geom_point(aes(size  = '0.7',
                     shape = mod_rank),
                 show.legend = F) +
      scale_fill_viridis( limits = tile_df$z %>% range ) + 
      ggtitle( response ) +
      theme(title        = element_text(angle = 0, hjust = 0.5, size = 10),
            legend.title = element_text(size = 10),
            legend.text  = element_text(size = 12),
            legend.position = 'bottom',
            axis.title   = element_blank(),
            axis.text.x  = element_text(angle = 90, 
                                        hjust = 1,
                                        vjust = 0.5))
    
  }else{
      
    p_out <- tile_df %>% 
      ggplot( aes(model, species) ) +
      geom_tile(aes(fill = z), color = "white") +
      geom_point(aes(size  = '0.7',
                     shape = mod_rank),
                 show.legend = FALSE) +
      scale_fill_viridis( limits = tile_df$z %>% range ) + 
      ggtitle( response ) +
      theme(title        = element_text(angle = 0, hjust = 0.5, size = 10),
            legend.title = element_text(size = 10),
            legend.text  = element_text(size = 12),
            legend.position = 'bottom',
            axis.title   = element_blank(),
            axis.text.y = element_blank(),
            axis.text.x = element_text(angle = 90, 
                                        hjust = 1,
                                        vjust = 0.5))  
    
  }

}

# Precipitation
p1 <- p_loo(1)
p2 <- p_loo(2)
p3 <- p_loo(3)
p4 <- p_loo(4)

# store plots
tiff( 'results/mod_sel/looic/precip_looic.tiff', 
      unit="in", width=10, 
      height=6.3, res=600,compression="lzw" )

grid.arrange(
     nrow = 1,
     grobs = list(p1, p2, p3, p4),
     widths = c(1.6, 0.9, 0.9, 0.9),
     layout_matrix = rbind(c(1, 2, 3, 4)) )

dev.off()


# Air temperature
p5 <- p_loo(5)
p6 <- p_loo(6)
p7 <- p_loo(7)
p8 <- p_loo(8)

# store plots
tiff( 'results/mod_sel/looic/airt_looic.tiff', 
      unit="in", width=10, 
      height=6.3, res=600,compression="lzw" )

grid.arrange(
     nrow = 1,
     grobs = list(p5, p6, p7, p8),
     widths = c(1.6, 0.9, 0.9, 0.9),
     layout_matrix = rbind(c(1, 2, 3, 4)) )

dev.off()


# 
tiff( 'results/mod_sel/looic/looic.tiff', 
      unit="in", width=10, height=10, res=600,compression="lzw" )

grid.arrange(
     nrow = 2, ncol=4,
     grobs = list(p1, p2, p3, p4,
                  p5, p6, p7, p8),
     widths = c(1.6, 0.9, 0.9, 0.9),
     layout_matrix = rbind(c(1, 2, 3, 4),
                           c(5, 6, 7, 8)) )

dev.off()


# # se_looic by model
# ggplot(visual_df) +
#   geom_point( aes(model, se_looic) ) +
#   theme(title        = element_text(angle = 0, hjust = 0.5, size = 10),
#         legend.title = element_text(size = 10),
#         legend.text = element_text(size = 12),
#         axis.title = element_blank(),
#         axis.text.x = element_text(angle = 90, hjust = 1) ) +
#   ggsave( 'C:/Users/ac22qawo/Desktop/methods_2.21/2.28/se_looic_airt_mod.tiff',
#   # ggsave( 'C:/Users/ac22qawo/Desktop/methods_2.21/2.28/se_looic_prec_mod.tiff',
#           width=6.3, height=6.3, compression='lzw')
# 
# # se_looic by species
# visual_df %>% 
#   mutate( species = substr(species,1,15) ) %>% 
#   ggplot() +
#   geom_point( aes(species,se_looic) ) +
#   theme(title        = element_text(angle = 0, hjust = 0.5, size = 10),
#         legend.title = element_text(size = 10),
#         legend.text = element_text(size = 12),
#         axis.title = element_blank(),
#         axis.text.x = element_text(angle = 90, hjust = 1) ) #+
#   ggsave( 'C:/Users/ac22qawo/Desktop/methods_2.21/2.28/se_looic_airt_spp.tiff',
#   # ggsave( 'C:/Users/ac22qawo/Desktop/methods_2.21/2.28/se_looic_prec_spp.tiff',
#           width=6.3, height=6.3, compression='lzw')

  
# 4. ELPD plot  ------------------------------------------

# plot elpd values
p_elpd <- function(ii){
  
  response <- names(loo_l)[ii]

  tile_df <- loo_l[response][[1]] %>% 
    mutate( species = substr(species,1,25) )
  
  if(ii == 1 | ii == 5){
  
    p_out <- tile_df %>% 
      ggplot( aes(model, species) ) +
      geom_tile(aes(fill = elpd), color = "white") +
      geom_point(aes(size  = '0.7',
                     shape = mod_rank_elpd),
                 show.legend = F) +
      scale_fill_viridis( limits = tile_df$elpd %>% range ) + 
      ggtitle( response ) +
      theme(title        = element_text(angle = 0, hjust = 0.5, size = 10),
            legend.title = element_text(size = 10),
            legend.text  = element_text(size = 12),
            legend.position = 'bottom',
            axis.title   = element_blank(),
            axis.text.x  = element_text(angle = 90, 
                                        hjust = 1,
                                        vjust = 0.5))
    
  }else{
      
    p_out <- tile_df %>% 
      ggplot( aes(model, species) ) +
      geom_tile(aes(fill = elpd), color = "white") +
      geom_point(aes(size  = '0.7',
                     shape = mod_rank_elpd),
                 show.legend = FALSE) +
      scale_fill_viridis( limits = tile_df$elpd %>% range ) + 
      ggtitle( response ) +
      theme(title        = element_text(angle = 0, hjust = 0.5, size = 10),
            legend.title = element_text(size = 10),
            legend.text  = element_text(size = 12),
            legend.position = 'bottom',
            axis.title   = element_blank(),
            axis.text.y = element_blank(),
            axis.text.x = element_text(angle = 90, 
                                       hjust = 1,
                                       vjust = 0.5) )  
    
  }

}

# Precipitation
p1 <- p_elpd(1)
p2 <- p_elpd(2)
p3 <- p_elpd(3)
p4 <- p_elpd(4)

# store plots
tiff( 'results/mod_sel/looic/precip_elpd.tiff', 
      unit="in", width=10, 
      height=6.3, res=600,compression="lzw" )

grid.arrange(
     nrow = 1,
     grobs = list(p1, p2, p3, p4),
     widths = c(1.6, 0.9, 0.9, 0.9),
     layout_matrix = rbind(c(1, 2, 3, 4)) )

dev.off()
  
  
# Air temperature
p5 <- p_elpd(5)
p6 <- p_elpd(6)
p7 <- p_elpd(7)
p8 <- p_elpd(8)

# store plots
tiff( 'results/mod_sel/looic/airt_elpd.tiff', 
      unit="in", width=10, 
      height=6.3, res=600,compression="lzw" )

grid.arrange(
     nrow = 1, 
     grobs = list(p1, p2, p3, p4),
     widths = c(1.6, 0.9, 0.9, 0.9),
     layout_matrix = rbind(c(1, 2, 3, 4)) )

dev.off()
    



# store plots
tiff( 'results/mod_sel/looic/elpd.tiff', 
      unit="in", width=10, 
      height=10, res=600,compression="lzw" )

grid.arrange(
     nrow = 2, ncol=4,
     grobs = list(p1, p2, p3, p4,
                  p5, p6, p7, p8),
     widths = c(1.6, 0.9, 0.9, 0.9),
     layout_matrix = rbind(c(1, 2, 3, 4),
                           c(5, 6, 7, 8)) )

dev.off()
    

  
# waic ------------------------------------------------------------

# rank the models 
plot_waic <- function(sum_df,spp_v){
  
  sel_df <- sum_df %>% 
              # se_waic,looic,se_looic
              select(model,species,waic)

  # calculate delta of performance (for plotting)
  delta_perf <- function(x){
    
    # performance of null model
    null_perf <- sel_df %>%
                  subset( species == x) %>% 
                  subset( model == 'ctrl1') %>% 
                  .$waic
    
    # spit out dealta WAIC
    sel_df %>%
      subset( species == x) %>% 
      mutate( delta_waic = waic - null_perf)
      
  }
  
  delta_df <- lapply(spp_v, delta_perf) %>% bind_rows
  
  # get model ranks (still for plotting!)
  model_rank <- function(x){
    delta_df %>%
      subset( species == x) %>%
      arrange( waic ) %>% 
      mutate( mod_rank = c(1:nrow(.)) ) %>% 
      mutate( mod_rank = replace(mod_rank,
                                 mod_rank > 3,
                                 NA) ) %>% 
      mutate( mod_rank = as.character(mod_rank) ) %>% 
      select( model, species, delta_waic, mod_rank )
  }

  lapply(spp_v, model_rank) %>% bind_rows
  
}

# make WAIC plot
tile_waic <- function(mod_v,spp_names,spp_v,titl){
  
  # data frame for tile plotss
  tile_df <- Map(create_sum_df,sum_f_l,spp_names) %>% 
              bind_rows %>% 
              subset( model %in% mod_v ) %>% 
              plot_waic(spp_v) %>% 
              # right_join(sum_df) %>% 
              subset( species %in% spp_v ) %>% 
              mutate( model   = factor(model,
                                       levels = mod_v ) ) %>% 
              mutate( species = substr(species,1,15) )
  
  # store tile plots 
  ggplot(tile_df, aes(model, species)) +
    geom_tile(aes(fill = delta_waic), color = "white") +
    geom_point(aes(size  = '0.7',
                         shape = mod_rank) ) +
    scale_fill_gradient2( )
    scale_fill_viridis( limits = tile_df$delta_waic %>% range ) + #
    ggtitle('Log Lambda') +
    theme(title        = element_text(angle = 0, hjust = 0.5, size = 10),
          legend.title = element_text(size = 10),
          legend.text = element_text(size = 12),
          axis.title = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1)) +
    ggsave( paste0('C:/Users/ac22qawo/Desktop/methods_2.21/2.28/',titl),
            width=6.3,height=6.3,compression='lzw')

}

# plot it out
tile_waic(mod_labs_no_mb,spp_names,spp_lppd,
          'noMB_tile_waic.tiff')
tile_waic(mod_labs_mb,spp_names,spp_lppd,
          'mb_tile_waic.tiff')
tile_waic(mod_labs_mb_n,spp_names,spp_lppd,
          'mb_n_tile_waic.tiff')



# plot ALL MODELS
setwd('I:/sie/101_data_AC/results/precip')
sum_f_l   <- grep('mod_summ',
                  list.files(),value=T)
spp_names <- gsub('mod_summaries_|.csv','',sum_f_l)
tile_looic(mod_labs_mb_n,spp_names,spp_names,
           'all_prec_ll_looic.tiff.tiff')
tile_looic(mod_labs_no_mb,spp_names,spp_names,
           'noMB_all_prec_ll_looic.tiff.tiff')

setwd('I:/sie/101_data_AC/results/airt')
sum_f_l   <- grep('mod_summ',
                  list.files(),value=T)
spp_names <- gsub('mod_summaries_|.csv','',sum_f_l)
tile_looic(mod_labs_mb_n,spp_names,spp_names,
           'all_airt_ll_looic.tiff.tiff')
tile_looic(mod_labs_no_mb,spp_names,spp_names,
           'noMB_all_airt_ll_looic.tiff.tiff')



# model summaries --------------------------------------
file_l    <- grep('mod_summaries_',list.files(),value=T)
spp_names <- gsub('mod_summaries_|.csv','',file_l)

div_trans <- function(ii){
  data.table::fread(file_l[ii]) %>% 
    select(model,n_diverg,rhat_high) %>% 
    mutate( species = spp_names[ii] )
    # cia %>% 
    # select(n_diverg) %>% 
    # group_by(model) %>% 
    # summarise(div_n = sum(n_diverg) ) %>% 
    # ungroup %>% 
    # mutate( species = spp_names[ii] )
}

div_df <- lapply(1:length(spp_names), div_trans) %>% bind_rows

library(ggplot2)

div_df %>% 
  mutate( model =
            factor(model, levels = c('ctrl1',
                                     'yr1','yr2','yr3',
                                     'gaus','expp',
                                     'movb','movb_c',
                                     'gaus_n', 'expp_n',
                                     'simpl_n') ) ) %>% 
  ggplot() +
  geom_point( aes(x=model,y=n_diverg) ) +
  theme( axis.text.x=element_text(angle = 90)) + 
  ggsave( 'Divergent_trans.tiff',
          height=6.3,width=6.3,compression='lzw')

div_df %>% 
  mutate( model =
            factor(model, levels = c('ctrl1',
                                     'yr1','yr2','yr3',
                                     'gaus','expp',
                                     'movb','movb_c',
                                     'gaus_n', 'expp_n',
                                     'simpl_n') ) ) %>% 
  ggplot() +
  geom_point( aes(x=model,y=rhat_high) ) +
  theme( axis.text.x=element_text(angle = 90)) + 
  ggsave( 'rhat_high.tiff',
          height=6.3,width=6.3,compression='lzw')



# posteriors  ----------------------------------------------------
post_l    <- grep('posterior_',list.files(),value=T)
spp_names <- gsub('posterior_|.csv','',file_l)


ratid_df <- data.table::fread(post_l[9])

ratid_df %>%
  extract() %>% 
  as.data.frame %>%  
  # subset( model == 'movb') %>%
  subset( model == 'movb_c') %>%
  select(beta_1:beta_36) %>% 
  gather(month,beta, beta_1:beta_36) %>% 
  mutate( month = gsub('beta_','',month) ) %>% 
  group_by( month ) %>% 
  summarise( mean = mean(beta),
             maxx = quantile(beta,prob=0.9),
             minn = quantile(beta,prob=0.1) ) %>% 
  ungroup( ) %>% 
  mutate( month = as.numeric(month) ) %>% 
  arrange(month) %>% 
  ggplot() +
  geom_pointrange(aes(x    = month,
                      y    = mean,
                      ymin = minn,
                      ymax = maxx) ) +
  geom_hline(yintercept = 0)
  
fit_mb %>% 
  extract() %>% 
  as.data.frame %>%  
  select(beta.1:beta.36) %>% 
  gather(month,beta, beta.1:beta.36) %>% 
  mutate( month = gsub('beta.','',month) ) %>% 
  group_by( month ) %>% 
  summarise( mean = mean(beta),
             maxx = quantile(beta,prob=0.9),
             minn = quantile(beta,prob=0.1) ) %>% 
  ungroup( ) %>% 
  mutate( month = as.numeric(month) ) %>% 
  arrange(month) %>% 
  ggplot() +
  geom_pointrange(aes(x    = month,
                      y    = mean,
                      ymin = minn,
                      ymax = maxx) ) +
  geom_hline(yintercept = 0) + 
  ggsave('C:/Users/ac22qawo/Desktop/methods_2.21/Ratibida_aldo.tiff',
         height=6.3,width=6.3,compression='lzw')
  



# Paronychia jamesii
tiff('C:/Users/ac22qawo/Desktop/methods_2.21/paronychia_mu_sd.tiff',
     unit="in", width=6.3, height=6.3, res=600,compression="lzw")

paron_df <- data.table::fread(post_l[6])  
  
par(mfrow=c(2,1), mar=c(2,2,1,0.5))
paron_df %>% 
  subset( model == 'gaus') %>% 
  .$sens_mu %>% hist(main="Paronychia jamesii")

paron_df %>% 
  subset( model == 'gaus') %>% 
  .$sens_sd %>% hist(main="Paronychia jamesii")

dev.off()


# Opuntia macrorhiza
tiff('C:/Users/ac22qawo/Desktop/methods_2.21/Opuntia_mu_sd.tiff',
     unit="in", width=6.3, height=6.3, res=600,compression="lzw")

opunt_df <- data.table::fread(post_l[4])  
  
par(mfrow=c(2,1), mar=c(2,2,1,0.5))
opunt_df %>% 
  subset( model == 'gaus') %>% 
  .$sens_mu %>% hist(main="Opuntia sens_mu")

opunt_df %>% 
  subset( model == 'gaus') %>% 
  .$sens_sd %>% hist(main="Opuntia sens_sd")

dev.off()


# Orchis purpurea
tiff('C:/Users/ac22qawo/Desktop/methods_2.21/Orchis_mu_sd.tiff',
     unit="in", width=6.3, height=6.3, res=600,compression="lzw")

orchi_df <- data.table::fread(post_l[5])  
  
par(mfrow=c(2,1), mar=c(2,2,1,0.5))
orchi_df %>% 
  subset( model == 'gaus') %>% 
  .$sens_mu %>% hist(main="Orchis sens_mu")

orchi_df %>% 
  subset( model == 'gaus') %>% 
  .$sens_sd %>% hist(main="Orchis sens_sd")

dev.off()



# Psoralea tenuifolia
tiff('C:/Users/ac22qawo/Desktop/methods_2.21/psoralea_mu_sd.tiff',
     unit="in", width=6.3, height=6.3, res=600,compression="lzw")

# Psoralea tenuifolia
psor_df <- data.table::fread(post_l[7])  

par(mfrow=c(2,1), mar=c(2,2,1,0.5))
psor_df %>% 
  subset( model == 'gaus') %>% 
  .$sens_mu %>% hist(main="Psoralea tenu. sens_mu")

paron_df %>% 
  subset( model == 'gaus') %>% 
  .$sens_sd %>% hist(main="Psoralea tenu. sens_sd")

dev.off()



# Orchis purpurea
orchis_df <- data.table::fread(post_l[5])   

orchis_df %>% 
  subset( model == 'gaus_n') %>% 
  .$sens_mu %>% hist(main="Orchis purpurea")


orchis_df %>% 
  subset( model == 'gaus_n') %>% 
  select(theta_y_1:theta_y_3) %>%
  boxplot(main="Orchis purpurea")
 
# Lomatium
lomat_df <- data.table::fread(post_l[3])   

lomat_df %>% 
  subset( model == 'gaus_n') %>% 
  .$sens_mu %>% hist(main='Lomatium')

lomat_df %>% 
  select(theta_y_1:theta_y_3) %>% 
  boxplot(main='Lomatium')

# Purshia
pur_df   <- data.table::fread(post_l[8])   

pur_df %>% 
  subset( model == 'simpl_n') %>% 
  select(theta_m_1:theta_m_12) %>% 
  boxplot

pur_df %>% 
  subset( model == 'simpl_n') %>% 
  select(theta_y_1:theta_y_3) %>% 
  boxplot



subset(model == 'expp_n') %>% 
  select( beta_1:beta_36) %>% 
  gather( month, beta, beta_1:beta_36) %>% 
  mutate( month = gsub('beta_','',month) ) %>% 
  group_by( month ) %>% 
  summarise( mean = mean(beta),
             maxx = quantile(beta,prob=0.9),
             minn = quantile(beta,prob=0.1) ) %>% 
  ungroup( ) %>% 
  mutate( month = as.numeric(month) ) %>% 
  arrange(month) %>% 
  ggplot() +
  geom_pointrange(aes(x    = month,
                      y    = mean,
                      ymin = minn,
                      ymax = maxx) ) +
  geom_hline(yintercept = 0)


post_ss   <- data.table::fread(file_l[3])
post_ss %>% 
  extract() %>% 
  as.data.frame %>% 
  select( beta.1:beta.36) %>% 
  gather( month, beta, beta.1:beta.36) %>% 
  mutate( month = gsub('beta.','',month) ) %>% 
  group_by( month ) %>% 
  summarise( mean = mean(beta),
             maxx = quantile(beta,prob=0.9),
             minn = quantile(beta,prob=0.1) ) %>% 
  ungroup( ) %>% 
  mutate( month = as.numeric(month) ) %>% 
  arrange(month) %>% 
  ggplot() +
  geom_pointrange(aes(x    = month,
                      y    = mean,
                      ymin = minn,
                      ymax = maxx) ) +
  geom_hline(yintercept = 0) +
  ggsave( "radidiba_mov_beta_c.tiff",
          height=6.3,width=6.3,compression='lzw')
