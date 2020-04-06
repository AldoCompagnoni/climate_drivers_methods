# 
setwd('E:/work/sApropos/2019.6.11')

35/40

cases_df <- expand.grid( clim_var = c('precip','airt'),
                         resp     = c('log_lambda', 'surv', 
                                      'grow', 'fec'),
                         file     = c('crossval', 'diagnostics',
                                      'mod_summaries', 'posterior'),
                         spp      = 1:39,
                         stringsAsFactors = F ) %>% 
                mutate( clim_resp = paste0(resp,'_',clim_var) )


# csv data files
csv_l <- grep('.csv', list.files(), value=T)
wrk_l <- grep('workshace', list.files(), value=T)
        
  
grep('surv_airt.csv',csv_l,value=T)
  
# create df of results files 
check_avail_files <-function(ii){
  
  # check whether you find a file or not
  out <- grep( cases_df$clim_resp[ii], csv_l, value=T) %>% 
          grep( cases_df$file[ii],.,value=T ) %>% 
          grep( paste0('array_vr-[0-9]{1,7}-',
                       cases_df$spp[ii],'_'),
                .,
                value=T )
  
  if( length(out) == 1 ){
    return( cases_df[ii,] )
  } else{
    return( NULL )
  }
  
}

# available files (or not)
avail_df <- lapply(1:nrow(cases_df), check_avail_files) %>% bind_rows

# 
anti_join( cases_df, avail_df) %>% 
  subset( !(spp %in% c(35:39)) )


grep('log_lambda_precip.csv', csv_l, value=T) %>% 
  grep('posterior', ., value=T) %>% 
  length
grep('log_lambda_airt.csv', csv_l, value=T) %>% 
  grep('posterior', ., value=T) %>% 
  length

grep('grow_precip.csv', csv_l, value=T) %>% 
  grep('posterior', ., value=T) %>% 
  length
grep('grow_airt.csv', csv_l, value=T) %>% 
  grep('posterior', ., value=T) %>% 
  length

grep('surv_precip.csv', csv_l, value=T) %>% 
  grep('posterior', ., value=T) %>% 
  length
grep('surv_airt.csv', csv_l, value=T) %>% 
  grep('posterior', ., value=T) %>% 
  length

grep('fec_precip.csv', csv_l, value=T) %>% 
  grep('posterior', ., value=T) %>% 
  length
grep('fec_airt.csv', csv_l, value=T) %>% 
  grep('posterior', ., value=T) %>% 
  length



gsub('array_vr-{7}-')