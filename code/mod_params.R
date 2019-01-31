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
