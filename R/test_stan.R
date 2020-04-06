install.packages("rstan", repos = "https://cloud.r-project.org/", dependencies = TRUE)
library(rstan)
setwd('C:/CODE/climate_drivers_methods')

# data
schools_dat <- list(J = 8, 
                    y = c(28,  8, -3,  7, -1,  1, 18, 12),
                    sigma = c(15, 10, 16, 11,  9, 11, 10, 18))

fit <- stan(file = 'code/stan/8schools.stan', data = schools_dat)

