# produce summary graphs
source("code/format_data.R")
library(dplyr)
library(tidyr)
library(data.table)
library(magrittr)
library(testthat)
library(viridis)
library(egg)
options(stringsAsFactors = F )

today_date    <- gsub("-","_",Sys.time() %>% as.Date %>% as.character)
m_back        <- 36
interval      <- NULL
pre_chelsa    <- NULL # '_pre_chelsa'


# plot standardized betas ------------------------------------------------

# rad
betas_st_df <- read.csv('I:/sie/101_data_AC/betas_st.csv') %>%
                  mutate( species = substr(species,1,20) )

# standardized betas 
ggplot(betas_st_df, aes(x=species, y=b_st) ) +
  geom_point() + 
  facet_grid( ~ model) +
  geom_hline( yintercept = 0 ) +
  theme( axis.text    = element_text( angle = 90,
                                      size  = 5),
         axis.title.x = element_blank() ) +
  ylab( expression('Standardized '*beta) ) +
  ggsave('results/beta_stand.tiff',
         height=3,width=6.3,compression='lzw')
  
# absolute betas 
ggplot(betas_st_df, aes(x=species, y=beta) ) +
  geom_point() + 
  facet_grid( ~ model) +
  geom_hline( yintercept = 0 ) +
  theme( axis.text    = element_text( angle = 90,
                                      size  = 5),
         axis.title.x = element_blank() ) +
  ylab( expression('Standardized '*beta) ) +
  ggsave('results/beta.tiff',
         height=3,width=6.3,compression='lzw')
  
# sd of x_antecedent
ggplot(betas_st_df, aes(x=species, y=sd_x) ) +
  geom_point() + 
  facet_grid( ~ model) +
  theme( axis.text    = element_text( angle = 90,
                                      size  = 5),
         axis.title.x = element_blank() ) +
  ylab( 'sd of x_antecedent' ) +
  ggsave('results/sd_x.tiff',
         height=3,width=6.3,compression='lzw')


# residual standard deviation
ggplot(betas_st_df, aes(x=species, y=sd_y) ) +
  geom_point() + 
  facet_grid( ~ model) +
  geom_hline( yintercept = 0 ) +
  theme( axis.text    = element_text( angle = 90,
                                      size  = 5),
         axis.title.x = element_blank() ) +
  ylab( 'sd of y' ) +
  ggsave('results/sd_y.tiff',
         height=3,width=6.3,compression='lzw')

# plot histograms
beta_h_df <- betas_st_df %>% 
               gather(measure,value, b_st:sd_y) 
         
# calculate moments
beta_mom  <- beta_h_df %>% 
              group_by(measure) %>% 
              summarise( mean = mean(value),
                         med  = median(value) ) %>% 
              ungroup

# histograms, with means and medians
beta_h_df %>% 
  left_join(beta_mom) %>% 
  ggplot() + 
  geom_histogram( aes(value) ) +
  geom_vline( aes(xintercept = mean) ) +
  geom_vline( aes(xintercept = med),
              linetype = 2) +
  facet_grid( ~ measure) +
  ggsave('results/beta_hist.tiff',
         height=3,width=6.3,compression='lzw')

# exploratory
betas_st_df %>% 
  select(beta,sd_x,sd_y) %>% 
  pairs

# result of an ointere
ggplot(betas_st_df) +
  geom_point( aes(x=sd_x,y=sd_y) ) +
  ggsave('results/sd_y_vs_sd_x.tiff',
         height=3,width=6.3,compression='lzw')

betas_st_df %>% 
  subset( sd_x > 0.6)


# Summarize moving windows results by climate variable -------------------------------
par_post <- function(ii){
  
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
  sum_files <- list.files(res_folder)[grep("posterior", list.files(res_folder) )] %>% 
                  stringr::str_subset(resp_clim)
  # read files
  mod_summ  <- lapply(sum_files, function(x) fread(paste0(res_folder,"/",x)) ) %>%
                  setNames( gsub("posterior", "", sum_files ) ) %>%
                  setNames( gsub(paste0(resp_clim,".csv"), "", names(.) ) )
  # all model selection summaries
  all_sums  <- Map(function(x,y) tibble::add_column(x, species = y, .before = 1), 
                   mod_summ, names(mod_summ) ) %>% 
                  # selec ONLY these model selection variables
                  lapply(function(x) x %>% 
                                      dplyr::select(species, 
                                                    model, 
                                                    beta) ) %>% 
                  bind_rows %>% 
                  mutate( clim_var = clim_var )
        
  all_sums

}

# all models results
input_df    <- expand.grid( clim_var = c("precip","airt"),
                            response = "log_lambda", #c("surv","grow","fec",
                            interval = "",
                            stringsAsFactors = F)

# ALL model information
post_df  <- lapply(1:nrow(input_df), par_post) %>%
                  bind_rows %>% 
                  arrange(species) %>% 
                  mutate( species = gsub('^_','',species) )

# Beta means and sd
post_summ <- post_df %>% 
               subset( !(model %in% c('ctrl1','yr_bet')) ) %>% 
               group_by(species, model, clim_var) %>% 
               summarise( beta_m = mean(beta,na.rm=T),
                          beta_05 = quantile(beta,prob=0.05,na.rm=T), 
                          beta_95 = quantile(beta,prob=0.95,na.rm=T),
                          ) %>% 
               ungroup %>% 
               # beta's CI
               mutate( beta_ci = beta_95 - beta_05 ) %>% 
               # trim species names
               mutate( species = strtrim(species, 20) ) %>% 
               # make
               mutate( species  = as.factor(species),
                       model    = as.factor(model),
                       clim_var = as.factor(clim_var) )
  

# plots ----------------------------------------------------

# plot beta means vs. model
ggplot(post_summ, aes(model, beta_m)) + 
  geom_point( aes(color = clim_var) ) +
  scale_color_viridis(discrete=TRUE) +
  theme(axis.text.x  = element_text(angle = 70, hjust = 1) ) +
  geom_hline( yintercept = 0.5,  linetype = 'dashed' ) +
  geom_hline( yintercept = -0.5, linetype = 'dashed' ) +
  ylab("mean beta") + 
  ggsave('results/beta_m_by_mod.tiff',
         width = 6.3, height = 6.3)

# plot beta means vs. species
ggplot(post_summ, aes(species, beta_m)) + 
  geom_point( aes(color = clim_var) ) +
  scale_color_viridis(discrete=TRUE) +
  theme(axis.text.x  = element_text(angle = 70, hjust = 1) ) +
  geom_hline( yintercept = 0.5,  linetype = 'dashed' ) +
  geom_hline( yintercept = -0.5, linetype = 'dashed' ) +
  ggsave('results/betas_m_by_spp.tiff',
         width = 6.3, height = 5)


# plot beta ci width vs. model
ggplot(post_summ, aes(model, beta_ci)) + 
  geom_point( aes(color = clim_var) ) +
  scale_color_viridis(discrete=TRUE) +
  theme(axis.text.x  = element_text(angle = 70, hjust = 1) ) +
  ylab("mean beta") + 
  ggsave('results/beta_ciw_by_mod.tiff',
         width = 6.3, height = 6.3)

# plot beta ci width vs. species
ggplot(post_summ, aes(species, beta_ci)) + 
  geom_point( aes(color = clim_var) ) +
  scale_color_viridis(discrete=TRUE) +
  theme(axis.text.x  = element_text(angle = 70, hjust = 1) ) +
  ylab("mean beta") + 
  ggsave('results/beta_ciw_by_spp.tiff',
         width = 6.3, height = 5)
