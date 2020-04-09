library( tidyverse )
options(stringsAsFactors = F)

# root of slices
dir1 <- 'C:/Users/ac22qawo/sapropos/slices_sim_reg/' 
dir  <- 'C:/Users/ac22qawo/sapropos/slices_sim_glm/' 

# read all simulations at once
glm_df <- lapply( list.files(dir), 
                  function(x) read.csv( paste0(dir,x) ) ) %>% 
             bind_rows

sim_df <- lapply( list.files(dir1), 
                  function(x) read.csv( paste0(dir1,x) ) ) %>% 
             bind_rows %>% 
             mutate( mod = 'gaus' )

all_df <- bind_rows( sim_df, glm_df )

# root mean squared error
rmse        <- function(y, x) (y - x)^2 %>% mean %>% sqrt


# Just look at coefficients ------------------------

plot_coefs <- function( mod_name ){
  sim_df %>% 
    subset( mod == mod_name ) %>% 
    ggplot( ) +
    geom_jitter( aes( true, scrat,
                      color = as.factor(ds) ),
                 width = 0.005,
                 alpha = 0.5 ) +
      theme_minimal() +
      theme( legend.position = 'none') +
      labs( x = 'True coefficients',
            y = 'Recovered coefficients' ) +
      ggsave( paste0('results/simulations/regularization/horseshoe_sims/',mod_name,
                     '_coefs.tiff'),
              width = 6.3, height = 6.3, compression='lzw')
}


all_df %>% 
  ggplot( ) +
  geom_jitter( aes( true, scrat,
                    color = as.factor(ds) ),
               width = 0.005,
               alpha = 0.5 ) +
  facet_wrap( ~ mod ) +
  theme_minimal() +
  theme( legend.position = 'none') +
  labs( x = 'True coefficients',
        y = 'Recovered coefficients' ) +
  ggsave( 'results/simulations/regularization/horseshoe_sims/coefs.tiff',
          width = 6.3, height = 6.3, compression='lzw')


plot_coefs( 'beta' )
plot_coefs( 'gamma' )


# RMSE: GLMNET versus Regularized horseshoe ---------------

# prepare data frame
plot_df    <- sim_df %>% 
                dplyr::select(-scrat_sd) %>% 
                gather( type, beta_e, scrat:glmnet ) %>% 
                mutate( type = replace(type, type == 'scrat', 'Horseshoe') ) %>% 
                mutate( true_sl = round(glob_sc, 2) ) %>% 
                mutate( true_sl = replace( true_sl, 
                                           true_sl == 1,
                                           32 ) ) %>% 
                mutate( true_sl = replace( true_sl, 
                                           true_sl == 0.05,
                                           10 ) ) %>%           
                mutate( true_sl = replace( true_sl, 
                                           true_sl == 0.01,
                                           1 ) )

# True ZEROS
true_zeros <- plot_df %>% 
  subset( true == 0 ) %>%
  group_by( rep, beta_p, sl_df, glob_sc, type ) %>% 
  summarise( rmse = rmse( beta_e, true) ) %>% 
  ungroup %>%
  mutate( coef = 'Zero coefficients' )

# true positive numbers
true_pos <- plot_df %>% 
  subset( true != 0 ) %>%
  group_by( rep, beta_p, sl_df, glob_sc, type ) %>% 
  summarise( rmse = rmse( beta_e, true) ) %>% 
  ungroup %>%
  mutate( coef = 'Positive coefficients' )


bind_rows( true_zeros, true_pos) %>% 
  ggplot() +
  geom_violin( aes( type, rmse,
                    fill = 1) ) +
  facet_wrap( ~ coef ) +
  theme_minimal() +
  theme( legend.position = 'none',
         axis.text.x = element_text( size = 15 ),
         axis.title  = element_text( size = 15 ),
         strip.text  = element_text( size = 15 ) ) +
  labs( x = 'Model type',
        y = 'RMSE of estiamtes vs. true values' )  +
  ggsave( 'results/simulations/regularization/horseshoe_sims/RMSE_glmnet_vs_horse.tiff',
          width = 6.3, height = 6.3, compression='lzw')

# # True ZEROS
# plot_df %>% 
#   subset( true == 0 ) %>%
#   group_by( rep, beta_p, sl_df, glob_sc, type ) %>% 
#   summarise( rmse = rmse( beta_e, true) ) %>% 
#   ungroup %>% 
#   ggplot() +
#   geom_violin( aes( type, rmse,
#                     fill = 1) ) +
#   theme_minimal() +
#   theme( legend.position = 'none' ) +
#   labs( x = 'Model type',
#         y = 'RMSE: estiamtes of true zeros' ) +
#   ggsave( 'results/simulations/regularization/RMSE_zero_glmnet_vs_horse.tiff',
#           width = 6.3, height = 6.3, compression='lzw')
# 
# # True Effects
# plot_df %>% 
#   subset( true != 0 ) %>%
#   group_by( rep, beta_p, sl_df, glob_sc, type ) %>% 
#   summarise( rmse = rmse( beta_e, true) ) %>% 
#   ungroup %>% 
#   ggplot() +
#   geom_violin( aes( type, rmse,
#                     fill = 1) ) +
#   theme_minimal() +
#   theme( legend.position = 'none' ) +
#   labs( x = 'model type',
#         y = 'RMSE: estiamtes of true NONzeros' ) +
#   ggsave( 'results/simulations/regularization/RMSE_nonzero_glmnet_vs_horse.tiff',
#           width = 6.3, height = 6.3, compression='lzw')
# 
# # ALL effects
# plot_df %>% 
#   group_by( rep, beta_p, sl_df, glob_sc, type ) %>% 
#   summarise( rmse = rmse( beta_e, true) ) %>% 
#   ungroup %>% 
#   ggplot() +
#   geom_violin( aes( type, rmse,
#                     fill = 1) ) +
#   theme_minimal() +
#   theme( legend.position = 'none' ) +
#   labs( x = 'model type',
#         y = 'RMSE: estiamtes of ALL betas' ) +
#   ggsave( 'results/simulations/regularization/RMSE_all_glmnet_vs_horse.tiff',
#           width = 6.3, height = 6.3, compression='lzw')



# PRIOR effect on Bayesian estimation --------------------

# prior values for the "zero" coefficients
prior_zero <- plot_df %>% 
  subset( true == 0 ) %>% 
  group_by( beta_p, sl_df, true_sl, type, rep ) %>% 
  summarise( rmse = rmse( beta_e, true) ) %>% 
  ungroup %>% 
  mutate( sl_gl = paste0(sl_df, '-', true_sl) ) %>% 
  subset( type == 'Horseshoe' ) %>%
  mutate( coef = 'Zero coefficients' )

# prior values for the "positive" coefficients
prior_pos <- plot_df %>% 
  subset( true != 0 ) %>% 
  group_by( beta_p, sl_df, true_sl, type, rep ) %>% 
  summarise( rmse = rmse( beta_e, true) ) %>% 
  ungroup %>% 
  mutate( sl_gl = paste0(sl_df, '-', true_sl) ) %>% 
  subset( type == 'Horseshoe' ) %>%
  mutate( coef = 'Positive coefficients' )


# overall 
bind_rows( prior_zero, prior_pos ) %>% 
  ggplot() +
  geom_violin( aes( sl_gl, rmse,
                    fill = 1) ) +
  theme_minimal() +
  facet_wrap( ~ coef ) +
  theme( legend.position = 'none',
         axis.text.x = element_text( size = 15,
                                     angle = 75 ),
         axis.title  = element_text( size = 15 ),
         strip.text  = element_text( size = 15 )) +
  labs( x = 'Prior: slab deg.freed.-num. of true slopes',
        y = 'RMSE: estiamtes of true zero betas' ) +
  ggsave( 'results/simulations/regularization/horseshoe_sims/RMSE_gscale_slab_df.tiff',
          width = 6.3, height = 6.3, compression='lzw')

# 
# # zero values by slab DF
# plot_df %>% 
#   subset( true == 0 ) %>% 
#   group_by( beta_p, sl_df, true_sl, type, rep ) %>% 
#   summarise( rmse = rmse( beta_e, true) ) %>% 
#   ungroup %>% 
#   mutate( sl_gl = paste0(sl_df, '-', true_sl) ) %>% 
#   subset( type == 'Horseshoe' ) %>%
#   ggplot() +
#   geom_violin( aes( sl_gl, rmse,
#                     fill = 1) ) +
#   theme_minimal() +
#   theme( axis.text.x = element_text( angle = 70),
#          legend.position = 'none' ) +
#   labs( x = 'Prior: slab deg.freed.-num. of true slopes',
#         y = 'RMSE: estiamtes of true zero betas' ) +
#   ggsave( 'results/simulations/regularization/RMSE_zeros_slab_df.tiff',
#           width = 6.3, height = 6.3, compression='lzw')
# 
# # nonzero values by slab DF
# plot_df %>% 
#   subset( true != 0 ) %>%
#   group_by( beta_p, sl_df, true_sl, type, rep ) %>% 
#   summarise( rmse = rmse( beta_e, true) ) %>% 
#   ungroup %>% 
#   mutate( sl_gl = paste0(sl_df, '-', true_sl) ) %>% 
#   subset( type == 'Horseshoe' ) %>%
#   ggplot() +
#   geom_violin( aes( sl_gl, rmse,
#                     fill = 1) ) +
#   theme_minimal() +
#   theme( axis.text.x = element_text( angle = 70),
#          legend.position = 'none' ) +
#   labs( x = 'Prior: slab deg.freed.-num. of true slopes',
#         y = 'RMSE: estiamtes of true NON-zero betas' ) +
#   ggsave( 'results/simulations/regularization/RMSE_nonzeros_slab_df.tiff',
#           width = 6.3, height = 6.3, compression='lzw')
# 
# 
# # all values by slab DF
# plot_df %>% 
#   subset( true != 0 ) %>%
#   group_by( beta_p, sl_df, true_sl, type, rep ) %>% 
#   summarise( rmse = rmse( beta_e, true) ) %>% 
#   ungroup %>% 
#   mutate( sl_gl = paste0(sl_df, '-', true_sl) ) %>% 
#   ggplot() +
#   geom_violin( aes( sl_gl, rmse,
#                     fill = 1) ) +
#   theme_minimal() +
#   theme( axis.text.x = element_text( angle = 70),
#          legend.position = 'none' ) +
#   labs( x = 'Prior: slab deg.freed.-num. of true slopes',
#         y = 'RMSE: estiamtes of true NON-zero betas' ) +
#   ggsave( 'results/simulations/regularization/RMSE_all_priors.tiff',
#           width = 6.3, height = 6.3, compression='lzw')
# 
