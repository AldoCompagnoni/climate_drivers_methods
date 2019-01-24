library(rstan)
rstan_options(auto_write = TRUE)

# 
writeLines(readLines("code/eight_schools_cp.stan"))

# customizations?
c_light <- c("#DCBCBC")
c_light_highlight <- c("#C79999")
c_mid <- c("#B97C7C")
c_mid_highlight <- c("#A25050")
c_dark <- c("#8F2727")
c_dark_highlight <- c("#7C0000")

# read data and fit model
input_data  <- read_rdump("code/eight_schools.data.R")
fit_cp      <- stan(file='code/eight_schools_cp.stan', data=input_data,
                    iter=1200, warmup=500, chains=1, 
                    seed=483892929, refresh=1200)

# Examine output --------------------------------------------
params_cp <- as.data.frame(extract(fit_cp, permuted=FALSE))
names(params_cp) <- gsub("chain:1.", "", names(params_cp), fixed = TRUE)
names(params_cp) <- gsub("[", ".", names(params_cp), fixed = TRUE)
names(params_cp) <- gsub("]", "", names(params_cp), fixed = TRUE)
params_cp$iter <- 1:700

par(mar = c(4, 4, 0.5, 0.5))
plot(params_cp$iter, log(params_cp$tau), col=c_dark, pch=16, cex=0.8,
     xlab="Iteration", ylab="log(tau)", ylim=c(-6, 4))

# calculate running mean 
running_means <- sapply(params_cp$iter, function(n) mean(log(params_cp$tau)[1:n]))
par(mar = c(4, 4, 0.5, 0.5))
plot(params_cp$iter, running_means, col=c_dark, pch=16, cex=0.8, ylim=c(0, 2),
    xlab="Iteration", ylab="MCMC mean of log(tau)")
# different from true mean
abline(h=0.7657852, col="grey", lty="dashed", lwd=3)

# Count divergent transitions
divergent <- get_sampler_params(fit_cp, inc_warmup=FALSE)[[1]][,'divergent__']
sum(divergent)

