# trying to understand what moment matching does in beta
library(ggthemes)
library(ggplot2)

# beta moment matching ---------------------------------

a_param <- function( mu, sig ){
  
  # mu * ( ((mu*(1 - mu)) / sig) - 1 )
  # (((1 - mu) / sig) - (1 / mu)) * mu ^ 2
  # (mu * ( sig + mu^2 - mu) ) / sig
  (mu * (sig + mu^2 - mu)) / sig
  
}

mu_v <- runif(10000, 0, 1)
si_v <- runif(10000, 0.001, 0.25)

y    <- Map( a_param, mu_v, si_v ) %>% unlist

df   <- data.frame( y = y,
                    mu = mu_v,
                    si = si_v,
                    z = 0 ) %>% 
          mutate( z = replace(z, y<0, 1) ) %>% 
          mutate( z = as.factor(z) )

df %>% 
  ggplot() +
  geom_point(aes(mu, si,
                color = z) ) +
  scale_color_discrete()

mu_v[which(y<0)[1]]
si_v[which(y<0)[1]]


# gamma moment matching ---------------------------------


gamma_a_par <- function( mu, sig ){

  mu / sig
  
}

mu_v <- runif(10000, 0, 100)
si_v <- runif(10000, 0.001, 1)

y    <- Map( gamma_a_par, mu_v, si_v ) %>% unlist
min( y )

df   <- data.frame( y = y,
                    mu = mu_v,
                    si = si_v,
                    z = 0 ) %>% 
  mutate( z = replace(z, y<0, 1) ) %>% 
  mutate( z = as.factor(z) )

df %>% 
  ggplot() +
  geom_point(aes(mu, si,
                 color = z) ) +
  scale_color_discrete()
