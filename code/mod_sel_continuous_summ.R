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
  res_folder<- paste0("C:/cloud/Dropbox/sAPROPOS/results/moving_windows/",
                      response,
                      "/",
                      clim_var,pre_chelsa,interval) 
  sum_files <- list.files(res_folder)[grep("mod_summaries_", list.files(res_folder) )] %>% 
                  stringr::str_subset(resp_clim)
  mod_summ  <- lapply(sum_files, function(x) read.csv(paste0(res_folder,"/",x)) ) %>%
                  setNames( gsub("mod_summaries_", "", sum_files ) ) %>%
                  setNames( gsub(paste0(resp_clim,".csv"), "", names(.) ) )
  all_sums  <- Map(function(x,y) tibble::add_column(x, species = y, .before = 1), 
                     mod_summ, names(mod_summ) ) %>% 
                  lapply(function(x) dplyr::select(x,species, model, waic, looic, mse, elpd)) %>% 
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
  
  # spit out 
  left_join(all_sums, rep_n_df) %>% 
    mutate(mean_elpd = elpd / rep_n,
           clim_var  = clim_var,
           response  = response)

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


# Boxplot: prop_dist_to_best_model -------------------------------------
relative_bp <- function(resp, clim_v){
  
  perf_by_spp <- mod_perf_df %>% 
                    subset( response == resp & clim_var == clim_v ) %>% 
                    dplyr::select(species, model, mean_elpd) %>% 
                    mutate( mean_elpd = replace(mean_elpd, mean_elpd == -Inf, NA)) %>% 
                    spread(model, mean_elpd)
  perf_mat    <- dplyr::select(perf_by_spp, -species) %>% t
  
  # get minimum of maximum values
  get_min     <- function(ii) which(perf_mat[,ii] == min(perf_mat[,ii],na.rm=T))
  get_max     <- function(ii) which(perf_mat[,ii] == max(perf_mat[,ii],na.rm=T))
  
  # species sequence
  spp_seq     <- 1:ncol(perf_mat)
  
  # matrix of benchmark ids
  max_ids   <- sapply(spp_seq ,get_max) %>% 
                    cbind( spp_seq )
  min_ids   <- sapply(spp_seq, get_min) %>% 
                    cbind( spp_seq )
  
  # benchmark values
  max_val   <- perf_mat[max_ids]
  min_val   <- perf_mat[min_ids]
  scale_val <- max_val - min_val
  
  # differences matrix 
  diff_mat    <- sweep(perf_mat, 2,  max_val, FUN='-') %>% 
                    sweep(2, scale_val, FUN='/')*-1 
    
  # labels for models to plot on x-axis
  mod_labs <- c('ctrl1','ctrl2','yr1','yr2','yr3','yr_bet',
                'yr_wgt','gaus','expp','gev',
                'expp_n','gev_n','simpl_n')
  
  # differences data frame
  long_df <- diff_mat %>% 
                as.data.frame %>% 
                setNames( perf_by_spp$species ) %>% 
                tibble::add_column(model = row.names(.), .before=1) %>% 
                gather(species,mean_elpd,Actaea_spicata:Thelesperma_megapotamicum) %>% 
                mutate( model = factor(model, levels = mod_labs) )
            
  

  boxplot( mean_elpd ~ model, data = long_df,
         xaxt = "n",  xlab = "", ylab = "",
         main = paste0( resp, '~', clim_v) )
  
  # x axis with ticks but without labels
  axis(1, at=1:13, labels = FALSE)
    
  # Plot x labs at default x position
  text(x =  seq_along(mod_labs), 
       y = par("usr")[3] - 0.05, srt = 45, adj = 1,
       labels = mod_labs, xpd = TRUE)
      
}


tiff("results/prop_dist_to_best_model.tiff",
     unit="in", width=6.3, height=6.3, res=600,compression="lzw")

par(mfrow=c(4,2), mar= c(3, 1.5, 1, 0.2) + 0.1, oma = c(0,3.5,0,0),
    mgp = c(2,0.5,0))
relative_bp('surv','airt')
relative_bp('surv','precip')

relative_bp('grow','airt')
mtext("Proportion of distance from best model (mean LPPD)", 2,line=2.5, at = -0.3,
      cex = 1.3, xpd = NA)
relative_bp('grow','precip')

relative_bp('fec','airt')
relative_bp('fec','precip')

relative_bp('log_lambda','airt')
relative_bp('log_lambda','precip')

dev.off()


# Boxplot: mean LPPD by response / climate variable combination ---------------------------------
tiff("results/mean_lppd_by_resp_climvar.tiff",
     unit="in", width=6.3, height=6.3, res=600,compression="lzw")

# https://stackoverflow.com/questions/12995683/any-way-to-make-plot-points-in-scatterplot-more-transparent-in-r
add_tr <- function(color,trans){
  
  # This function adds transparancy to a color.
  # Define transparancy with an integer between 0 and 255
  # 0 being fully transparant and 255 being fully visable
  # Works with either color and trans a vector of equal length,
  # or one of the two of length 1.
  
  if (length(color)!=length(trans)&!any(c(length(color),length(trans))==1)) stop("Vector lengths not correct")
  if (length(color)==1 & length(trans)>1) color <- rep(color,length(trans))
  if (length(trans)==1 & length(color)>1) trans <- rep(trans,length(color))
  
  num2hex <- function(x)
  {
    hex <- unlist(strsplit("0123456789ABCDEF",split=""))
    return(paste(hex[(x-x%%16)/16+1],hex[x%%16+1],sep=""))
  }
  rgb <- rbind(col2rgb(color),trans)
  res <- paste("#",apply(apply(rgb,2,num2hex),2,paste,collapse=""),sep="")
  return(res)
}


mod_perf_df <- mod_perf_df %>%
                  mutate( mean_elpd = replace(mean_elpd, mean_elpd == -Inf, NA) )

par(mfrow=c(1,1), mar = c(3,3,0.2,0.2), oma = c(0,0,0,0),
    mgp=c(1.5,0.5,0), new=F)
boxplot(mean_elpd ~  clim_var + response, data = mod_perf_df,
        ylab = "Mean LPPD", xaxt='n',cex.lab=1.3, pch=NA)

coord_df <- expand.grid(response = c('fec','grow','log_lambda','surv'),
                        clim_var = c('airt','precip'),
                        stringsAsFactors = F) %>% 
              arrange(response, clim_var) %>% 
              mutate( x = c(1:8) )

points(mean_elpd ~  jitter(x,factor=1.5), data = left_join(mod_perf_df, coord_df),
       new=T, col = add_tr('grey50',50), pch = 16)

axis(1,at=1:8,labels=F)
text(1:8,y=-23.8,c('fec','fec','grow','grow','lambda','lambda','surv','surv'), xpd=NA)
text(1:8,y=-24.8,c('(airt)','(prec)','(airt)','(prec)','(airt)','(prec)','(airt)','(prec)'), xpd=NA)

abline(h=0,lty=2)

dev.off()


# Boxplot: absolute LPPD by species -------------------------------------------------------------
abs_lppd <- function(resp, clim_v){
  
  # absolute measures
  boxp_df <- mod_perf_df %>% 
                subset( response == resp & clim_var == clim_v) %>% 
                mutate( mean_elpd = replace(mean_elpd,
                                            mean_elpd == -Inf,
                                            NA))
  
  if( grepl('grow|surv', resp) ){
  
    boxplot(mean_elpd ~ species, data=boxp_df,
            boxlwd = 0.5, medlwd = 1, outlwd = 0.5, whisklwd = 0.5,
            outcex = 0.5,
            xaxt = 'n', ylim = c(-3,6), cex.main = 1.5,
            main = paste0(resp, ' ~ ', clim_v) )
    # x axis with ticks but without labels
    axis(1, labels = F, tick=F)
    abline(h=0, lty = 2, lwd = 0.5)  
    abline(v=5, lty = 2, lwd = 0.5)  
    abline(v=15, lty = 2, lwd = 0.5)
    abline(v=25, lty = 2, lwd = 0.5)
  
  }else{
    
    boxplot(mean_elpd ~ species, data=boxp_df,
            boxlwd = 0.5, medlwd = 1, outlwd = 0.5, whisklwd = 0.5,
            outcex = 0.5,
            xaxt = 'n', ylim = c(-5,10), cex.main = 1.5,
            main = paste0(resp, ' ~ ', clim_v) )
    # x axis with ticks but without labels
    axis(1, labels = F, tick=F)
    abline(h=0, lty = 2, lwd = 0.5)  
    abline(v=5, lty = 2, lwd = 0.5)  
    abline(v=15, lty = 2, lwd = 0.5)
    abline(v=25, lty = 2, lwd = 0.5)

  }
  
}

# plot species labels
plot_spp_lab <- function(x){

  spp_labs <- unique(x$species) %>% 
                stringr::str_replace('var._gnaphalifolium_2','...')
  
  # Plot x labs at default x position
  text(x =  seq_along(spp_labs), 
       y = -6, srt = 90, adj = 1,
       labels = spp_labs, xpd = NA, cex = 0.7)
}



tiff("results/mean_lppd_by_spp.tiff",
     unit="in", width=6.3, height=6.3, res=600,compression="lzw")

par(mfrow=c(3,2), mar = c(0,2,2,0.2), oma = c(8,2,0,0),
      mgp=c(1.7,0.5,0))

abs_lppd('surv', 'airt')
abs_lppd('surv', 'precip')

abs_lppd('grow', 'airt')
mtext("Mean LPPD", 2,line=2, at = 3,
      cex = 1.3, xpd = NA)
abs_lppd('grow', 'airt')

abs_lppd('fec', 'airt')
plot_spp_lab(mod_perf_df)
abs_lppd('fec', 'airt')
plot_spp_lab(mod_perf_df)

dev.off()



# # elpd heat map
# elpd_heat_map <- function(resp, clim_v){
# 
#   perf_by_spp <- mod_perf_df %>%
#                     subset( response == resp & clim_var == clim_v ) %>%
#                                   dplyr::select(species, model, mean_elpd) %>%
#                                   spread(model, mean_elpd)
#   perf_mat    <- dplyr::select(perf_by_spp, -species) %>% t
#   perf_mat[perf_mat == -Inf] <- NA
# 
#   print(lattice::levelplot(t(perf_mat), col.regions=heat.colors))
# 
# }



# Tileplot: relative (%) distances to best model -----------------------
prop_heat_map <- function(resp, clim_v){

  perf_by_spp <- mod_perf_df %>%
                    subset( response == resp & clim_var == clim_v ) %>%
                    dplyr::select(species, model, rep_yr, rep_n, mean_elpd) %>% 
                    mutate( mean_elpd = replace(mean_elpd, mean_elpd == -Inf, NA)) %>%
                    spread(model, mean_elpd)
  perf_mat    <- dplyr::select(perf_by_spp, -species, -rep_yr, -rep_n) %>% t
  
  # get minimum of maximum values
  get_min     <- function(ii) which(perf_mat[,ii] == min(perf_mat[,ii],na.rm=T))
  get_max     <- function(ii) which(perf_mat[,ii] == max(perf_mat[,ii],na.rm=T))
  
  # species sequence
  spp_seq     <- 1:ncol(perf_mat)
  
  # matrix of benchmark ids
  max_ids   <- sapply(spp_seq ,get_max) %>%
                    cbind( spp_seq )
  min_ids   <- sapply(spp_seq, get_min) %>%
                    cbind( spp_seq )
  
  # benchmark values
  max_val   <- perf_mat[max_ids]
  min_val   <- perf_mat[min_ids]
  scale_val <- max_val - min_val
  
  # differences matrix
  diff_mat    <- sweep(perf_mat, 2,  max_val, FUN='-') %>%
                    sweep(2, scale_val, FUN='/')*-1
  
  # labels for models to plot on x-axis
  mod_labs <- c('ctrl1','ctrl2','yr1','yr2','yr3','yr_bet',
                'yr_wgt','gaus','expp','gev',
                'expp_n','gev_n','simpl_n')
  
  # differences data frame
  long_df <- diff_mat %>%
                as.data.frame %>%
                setNames( perf_by_spp$species ) %>%
                tibble::add_column(model = row.names(.), .before=1) %>%
                gather(species,mean_elpd,Actaea_spicata:Thelesperma_megapotamicum) %>%
                mutate( model = factor(model, levels = mod_labs) ) %>%
                left_join( dplyr::select(perf_by_spp, species, rep_yr, rep_n) ) %>% 
                mutate( species = replace(species,
                                          grepl('Eriogonum',species),
                                          'Eriogonum_longifolium...') ) %>% 
                arrange( rep_yr, rep_n, species, model ) %>% 
                mutate( species = factor(species, levels = unique(species)) )
  
  # plot it out
  ggplot(long_df, aes(model, species, rep_yr)) +
  geom_tile(aes(fill = mean_elpd), color = "white") +
  scale_fill_viridis() + #scale_fill_gradient2(low ="white", high = "steelblue") +
  ylab( "Species" ) +
  xlab( paste0('Model: ',resp," ~ ",clim_v) ) +
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 12),
        plot.title = element_text(size=16),
        axis.title=element_text(size=14,face="bold"),
        axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(fill = "Mean LPPD") +
  ggsave(paste0("results/prop_diff_lppd_",
             resp,"_",clim_v,".tiff"),
         width = 5, height = 5, compression="lzw") #res=600,

  rm( list=c('perf_mat', 'spp_seq', 'max_ids', 'min_ids', 'max_val', 
             'diff_mat', 'long_df') )
  
}

prop_heat_map('surv','airt')
prop_heat_map('surv','precip')
prop_heat_map('grow','airt')
prop_heat_map('grow','precip')
prop_heat_map('fec','airt')
prop_heat_map('fec','precip')
prop_heat_map('log_lambda','airt')
prop_heat_map('log_lambda','precip')


# Tileplot: differences in mean_lppd to best model -------------------------------
mean_lppd_diff_heat_map <- function(resp, clim_v){

  perf_by_spp <- mod_perf_df %>%
                    subset( response == resp & clim_var == clim_v ) %>%
                    dplyr::select(species, model, rep_yr, rep_n, mean_elpd) %>% 
                    mutate( mean_elpd = replace(mean_elpd, mean_elpd == -Inf, NA)) %>%
                    spread(model, mean_elpd)
  perf_mat    <- dplyr::select(perf_by_spp, -species, -rep_yr, -rep_n) %>% t
  
  # get minimum of maximum values
  get_max     <- function(ii) which(perf_mat[,ii] == max(perf_mat[,ii],na.rm=T))
  
  # species sequence
  spp_seq     <- 1:ncol(perf_mat)
  
  # matrix of benchmark ids
  max_ids   <- sapply(spp_seq ,get_max) %>%
                    cbind( spp_seq )
  
  # benchmark values
  max_val   <- perf_mat[max_ids]
  
  # differences matrix
  diff_mat    <- sweep(perf_mat, 2,  max_val, FUN='-') 
  
  # labels for models to plot on x-axis
  mod_labs <- c('ctrl1','ctrl2','yr1','yr2','yr3','yr_bet',
                'yr_wgt','gaus','expp','gev',
                'expp_n','gev_n','simpl_n')
  
  # differences data frame
  long_df <- diff_mat %>%
                as.data.frame %>%
                setNames( perf_by_spp$species ) %>%
                tibble::add_column(model = row.names(.), .before=1) %>%
                gather(species,mean_elpd,Actaea_spicata:Thelesperma_megapotamicum) %>% 
                mutate( model = factor(model, levels = mod_labs) ) %>%
                left_join( dplyr::select(perf_by_spp, species, rep_yr, rep_n) ) %>% 
                mutate( species = replace(species,
                                          grepl('Eriogonum',species),
                                          'Eriogonum_longifolium...') ) %>% 
                arrange( rep_yr, rep_n, species, model ) %>% 
                mutate( species = factor(species, levels = unique(species)) )
  
  # plot it out
  ggplot(long_df, aes(model, species)) +
  geom_tile(aes(fill = mean_elpd), color = "white") +
  scale_fill_viridis() + #scale_fill_gradient2(low ="white", high = "steelblue") +
  ylab( "Species" ) +
  xlab( paste0('Model: ',resp," ~ ",clim_v) ) +
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 12),
        plot.title = element_text(size=16),
        axis.title=element_text(size=14,face="bold"),
        axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(fill = expression(Delta*" Mean LPPD") ) +
  ggsave(paste0("results/mean_lppd_diff_",
             resp,"_",clim_v,".tiff"),
         width = 5, height = 5, compression="lzw") #res=600,

  rm( list=c('perf_mat', 'spp_seq', 'max_ids', 'max_val', 
             'diff_mat', 'long_df') )
  
}

mean_lppd_diff_heat_map('surv','airt')
mean_lppd_diff_heat_map('surv','precip')
mean_lppd_diff_heat_map('grow','airt')
mean_lppd_diff_heat_map('grow','precip')
mean_lppd_diff_heat_map('fec','airt')
mean_lppd_diff_heat_map('fec','precip')
mean_lppd_diff_heat_map('log_lambda','airt')
mean_lppd_diff_heat_map('log_lambda','precip')



# tile plot of lppd differences from best model
lppd_diff_heat_map <- function(resp, clim_v){

  # mod_perf_df <- mod_perf_df %>% 
  #                   mutate( elpd = replace(elpd, 
  #                                          species  == "Orchis_purpurea" & 
  #                                          response == "log_lambda" & 
  #                                          clim_var == "airt" &
  #                                          model    == 'yr_wgt', NA) )
  
  perf_by_spp <- mod_perf_df %>%
                    subset( response == resp & clim_var == clim_v ) %>%
                    dplyr::select(species, model, rep_yr, rep_n, elpd) %>% 
                    mutate( elpd = replace(elpd, elpd == -Inf, NA)) %>%
                    spread(model, elpd)
  perf_mat    <- dplyr::select(perf_by_spp, -species, -rep_yr, -rep_n) %>% t
  
  # get minimum of maximum values
  get_max     <- function(ii) which(perf_mat[,ii] == max(perf_mat[,ii],na.rm=T))
  
  # species sequence
  spp_seq     <- 1:ncol(perf_mat)
  
  # matrix of benchmark ids
  max_ids     <- sapply(spp_seq ,get_max) %>%
                        cbind( spp_seq )
  
  # benchmark values
  max_val     <- perf_mat[max_ids]
  
  # differences matrix
  diff_mat    <- sweep(perf_mat, 2,  max_val, FUN='-') 
  
  # labels for models to plot on x-axis
  mod_labs    <- c('ctrl1','ctrl2','yr1','yr2','yr3','yr_bet',
                   'yr_wgt','gaus','expp','gev',
                   'expp_n','gev_n','simpl_n')
  
  # differences data frame
  long_df <- diff_mat %>%
                as.data.frame %>%
                setNames( perf_by_spp$species ) %>%
                tibble::add_column(model = row.names(.), .before=1) %>%
                gather(species,elpd,Actaea_spicata:Thelesperma_megapotamicum) %>% 
                mutate( model = factor(model, levels = mod_labs) ) %>%
                left_join( dplyr::select(perf_by_spp, species, rep_yr, rep_n) ) %>% 
                mutate( species = replace(species,
                                          grepl('Eriogonum',species),
                                          'Eriogonum_longifolium...') ) %>% 
                arrange( rep_yr, rep_n, species, model ) %>% 
                mutate( species = factor(species, levels = unique(species)) )
  
  # plot it out
  ggplot(long_df, aes(model, species)) +
  geom_tile(aes(fill = elpd), color = "white") +
  # scale_fill_viridis() + #
  scale_fill_viridis( limits = c(-30,0) ) + #
  ylab( "Species" ) +
  xlab( paste0('Model: ',resp," ~ ",clim_v) ) +
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 12),
        plot.title = element_text(size=16),
        axis.title=element_text(size=14,face="bold"),
        axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(fill = expression(Delta*" LPPD") ) +
  ggsave(paste0("results/lppd_diff_",
             resp,"_",clim_v,".tiff"), 
         width = 5, height = 5, compression="lzw") #res=600,

  rm( list=c('perf_mat', 'spp_seq', 'max_ids', 'max_val', 
             'diff_mat', 'long_df') )
  
}

lppd_diff_heat_map('surv','airt')
lppd_diff_heat_map('surv','precip')
lppd_diff_heat_map('grow','airt')
lppd_diff_heat_map('grow','precip')
lppd_diff_heat_map('fec','airt')
lppd_diff_heat_map('fec','precip')
lppd_diff_heat_map('log_lambda','airt')
lppd_diff_heat_map('log_lambda','precip')




# Tileplot: mean lppd -------------------------------------------------------
abs_heat_map <- function(resp, clim_v){

  # limits_vec  <- mod_perf_df %>% 
  #                   mutate( elpd = replace(elpd, elpd == -Inf, NA)) %>% 
  #                   .$elpd %>% range(na.rm=T)
      
  perf_by_spp <- mod_perf_df %>%
                    mutate( elpd = replace(elpd, elpd == -Inf, NA)) %>%
                    subset( response == resp & clim_var == clim_v ) %>%
                    mutate( mean_elpd = elpd / rep_n) %>% 
                    dplyr::select(species, model, rep_p, rep_yr, rep_n, mean_elpd) %>% 
                    arrange( rep_p, rep_yr, species, model ) %>% 
                    mutate( species = replace(species,
                                              grepl('Eriogonum',species),
                                              'Eriogonum_longifolium...') ) %>% 
                    mutate( species = factor(species, levels = unique(species)) )
  
  # plot it out
  ggplot(perf_by_spp, aes(model, species, rep_yr)) +
  geom_tile(aes(fill = mean_elpd), color = "white") +
  scale_fill_viridis( ) +
  ylab( "Species" ) +
  xlab( paste0('Model: ',resp," ~ ",clim_v) ) +
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 12),
        plot.title = element_text(size=16),
        axis.title=element_text(size=14,face="bold"),
        axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(fill = "Mean LPPD") +
  ggsave(paste0("results/abs_lppd_",
             resp,"_",clim_v,".tiff"), 
         width = 5, height = 5, compression="lzw") #res=600,

}

abs_heat_map('surv','airt')
abs_heat_map('surv','precip')
abs_heat_map('grow','airt')
abs_heat_map('grow','precip')
abs_heat_map('fec','airt')
abs_heat_map('fec','precip')
abs_heat_map('log_lambda','airt')
abs_heat_map('log_lambda','precip')



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
  
  find_best_mod <- function(x){
    long_df %>%
      subset( species == x) %>%
      subset( elpd == max(elpd) ) %>%
      select( model, species )
  }
  
  model_rank <- function(x){
    long_df %>%
      subset( species == x) %>%
      arrange( desc(elpd) ) %>% 
      mutate( mod_rank = c(1:nrow(.)) ) %>% 
      mutate( mod_rank = replace(mod_rank,
                                 mod_rank > 3,
                                 NA) ) %>% 
      mutate( mod_rank = as.character(mod_rank) ) %>% 
      # subset( mod_rank == mod_r) %>% 
      select( model, species, mod_rank )
  }
  # find_best_mod <- function(x){
  #   long_df %>%
  #     subset( species == x) %>%
  #     mutate( mod_rank = rank(elpd) ) %>%
  #     select( model, mod_rank, species)
  # }
  
  # pick the best models
  mod_sel_df <- lapply(spp_v, model_rank) %>% 
                    bind_rows %>% 
                    right_join( long_df )
  # mod_2_df <- lapply(spp_v, model_rank, 2) %>% 
  #                   bind_rows %>% 
  #                   mutate( mod_2 = 2 ) %>% 
  #                   right_join( long_df )
  # mod_3_df <- lapply(spp_v, model_rank, 1) %>% 
  #                   bind_rows %>% 
  #                   mutate( mod_3 = 3 ) %>% 
  #                   right_join( long_df )
  
  # mod_sel_df <- Reduce( function(...) full_join(...),
  #                       list(mod_1_df,
  #                            mod_2_df,
  #                            mod_3_df) ) %>% 
  #                 mutate( mod_size = 1 )
  
  # best_mod_df %>% 
  #   select(model, mod_rank, species) %>% 
  #   spread( model, mod_rank ) %>% head
  
  return(mod_sel_df)
    
}


# four tile plot 
four_tile_plot <- function(format_function, fill_var, clim_v, var_lim, 
                           expr_fill, file_title){
  
  # plot it out
  p1 <- ggplot(format_function('surv', clim_v), aes(model, species)) +
        geom_tile(aes_string(fill = fill_var), color = "white") +
        geom_point(aes(size = 0.5,
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
        geom_point(aes(size = 0.5,
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
        geom_point(aes(size = 0.5,
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
        geom_point(aes(size = 0.5,
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


format_function <- form_diff_lppd_df
fill_var    <- 'elpd'
clim_v      <- 'airt'
var_lim     <- c(-30,0)
expr_fill   <- expression(Delta*" LPPD")
file_title  <- 'results/lppd_diff_vr_airt.tiff'


# rank tile plots
four_tile_plot(form_rank_lppd_df, 'rank', 'airt',   c(1,13),  
               'Model rank',    'rank_lppd_vr_airt.tiff')
four_tile_plot(form_rank_lppd_df, 'rank', 'precip', c(1,13),  
               'Model rank',    'rank_lppd_vr_precip.tiff')

# mean lppd
range_m_lppd  <- mod_perf_df$mean_elpd %>% range(na.rm=T)
four_tile_plot(form_mean_lppd_df, 'mean_elpd', 'airt',   c(-15,2), #range_m_lppd,
               "Mean LPPD",       'mean_lppd_vr_airt.tiff')
four_tile_plot(form_mean_lppd_df, 'mean_elpd', 'precip', c(-15,2), #range_m_lppd, 
               "Mean LPPD",        'mean_lppd_vr_precip.tiff')




# mean lppd plots
resp   <- 'fec'
clim_v <- 'precip'


elpd_mod_df %>% 
  group_by(species) %>% 
  summarise( )


mod_perf_df %>% subset( grepl('Trillium',species) & response=='fec' )

# format the data as we want
form_mean_lppd_df <- function(resp, clim_v){
  
  # data for specific response and climate variable
  elpd_mod_df <- mod_perf_df %>%
                    subset( response == resp & clim_var == clim_v ) %>%
                    dplyr::select(species, model, rep_yr, rep_n, mean_elpd) %>% 
                    mutate( mean_elpd = replace(mean_elpd, mean_elpd == -Inf, NA) )
      
  # long df with model ranks!
  mean_elpd_df <- elpd_mod_df %>% 
                    full_join( spp_df ) %>% 
                    mutate( model = factor(model, levels = mod_labs) ) %>%
                    arrange( rep_yr, rep_n, species, model )  %>% 
                    mutate( species = replace(species,
                                              grepl('Eriogonum',species),
                                              'Eriogonum_longifolium...') ) %>% 
                    mutate( species = factor(species, levels = unique(species)) ) 
  
  return(mean_elpd_df)
    
}



# format the data as we want
form_rank_lppd_df <- function(resp, clim_v){
  
  # data for specific response and climate variable
  elpd_mod_df <- mod_perf_df %>%
                    subset( response == resp & clim_var == clim_v ) %>%
                    dplyr::select(species, model, rep_yr, rep_n, elpd) %>% 
                    mutate( elpd = replace(elpd, elpd == -Inf, NA) ) %>% 
                    arrange(species, elpd )
      
  # split data by species (to ricompose later)
  spp_mod_l   <- split(elpd_mod_df, elpd_mod_df$species)
  
  # long df with model ranks!
  mod_rank_df <- lapply(spp_mod_l, function(x) mutate(x, rank = 1:13) ) %>% 
                    # ricompose in one data frame!
                    bind_rows %>% 
                    full_join( spp_df ) %>% 
                    mutate( model = factor(model, levels = mod_labs) ) %>%
                    arrange( rep_yr, rep_n, species, model )  %>% 
                    mutate( species = replace(species,
                                              grepl('Eriogonum',species),
                                              'Eriogonum_longifolium...') ) %>% 
                    mutate( species = factor(species, levels = unique(species)) ) %>% 
                    mutate( rank = replace(rank, is.na(elpd), NA) )
  
  return(mod_rank_df)
    
}

