library(dplyr)

# get random species NOT from hays
spp_keep <- c('Actaea_spicata', 
              'Eryngium_alpinum',
              'Cirsium_pitcheri_8')

# Read CHELSA data from my own sApropos folder
# Randomly choose Actaea_spicata at Site A
airt <- read.csv('C:/cloud/Dropbox/sApropos/airt_chelsa_hays.csv') %>% 
          subset( species    == "Actaea_spicata" &
                  population == 'Site A' )
prec <- read.csv('C:/cloud/Dropbox/sApropos/precip_chelsa_hays.csv') %>% 
          subset( species    == "Actaea_spicata" &
                  population == 'Site A' )

# write data 
write.csv(airt , 'C:/CODE/climate_drivers_methods/data/demo_airt.csv', row.names=F)
write.csv(prec , 'C:/CODE/climate_drivers_methods/data/demo_prec.csv', row.names=F)
