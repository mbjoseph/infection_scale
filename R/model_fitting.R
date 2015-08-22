# fitting a simple random effect model to the rib abundance data
source('R/d_clean_joint.R')
source('R/helpers.R')
library(rstan)

stan_d <- list(n = nrow(pd), 
               nsite = max(pd$numsite), 
               nregion = max(pd$region), 
               site = pd$numsite, 
               nyear = length(unique(pd$year)),
               year = pd$year - min(pd$year) + 1,
               region = pd$region, 
               y = pd$rib, 
               nspec = length(unique(pd$speciescode)),
               h_spec = as.numeric(pd$speciescode),
               local_richness = c(scale(pd$local_rich)), 
               region_richness = c(scale(pd$reg_rich)), 
               isd = c(scale(log(pd$ISD + .1))))

watch <- c("eta_site", "eta_region", 
           "eta_species", "eta_year",
           "a0", 
           "sigma_site", "sigma_region", 
           "sigma_species", 
           "sigma_indiv",
           "beta_local", "beta_region", "beta_isd", "beta_intxn",
           "log_lik")

m_init <- stan('R/isd.stan', data=stan_d, chains=1, iter=1, pars=watch)
m_fit <- stan(fit=m_init, data=stan_d, pars=watch, chains=3, 
              iter=500, cores=3)

traceplot(m_fit, pars=watch[-c(1, length(watch))], inc_warmup=F)
waic(m_fit)

library(ggmcmc)
ggd <- ggs(m_fit)
ggs_caterpillar(ggd, 'eta_year')
ggs_caterpillar(ggd, 'eta_species')
ggs_caterpillar(ggd, 'eta_region')
ggs_caterpillar(ggd, 'eta_site')

post <- rstan::extract(m_fit)
nspec <- stan_d$nspec

par(mfrow=c(2, 2))
lims <- range(c(post$beta_local, post$beta_region, post$beta_isd))
br <- seq(lims[1], lims[2], length.out=40)
hist(post$beta_local, breaks=br, 
     main="Effect of local richness")
abline(v=0, lty=2, col="red")
hist(post$beta_region, breaks=br, 
     main="Effect of regional richness")
abline(v=0, lty=2, col="red")
hist(post$beta_isd, breaks=br, 
     main="Effect of infected snail density")
abline(v=0, lty=2, col="red")
hist(post$beta_isd2, breaks=br, 
     main="Effect of sq(infected snail density)")
abline(v=0, lty=2, col="red")

# plot relationship between isd and mu_infection
n <- 50
xvals <- seq(min(stan_d$isd), max(stan_d$isd), length.out=n)
calc_eff <- function(x, post){
  post$beta_isd * x + post$beta_isd2 * x^2
}
mu_eff <- array(dim=c(n, length(post$a0)))
for (i in 1:n){
  mu_eff[i, ] <- calc_eff(xvals[i], post)
}

med_eff <- apply(mu_eff, 1, median)
hdi_eff <- apply(mu_eff, 1, HDI)

par(mfrow=c(1, 1))
plot(xvals, med_eff, type='l', ylim=range(mu_eff))
lines(xvals, hdi_eff[1, ], lty=2)
lines(xvals, hdi_eff[2, ], lty=2)
rug(stan_d$isd)





# plot species specific intercepts
par(mfrow=c(2, 3), mar=c(5, 4, 4, 2) + 0.1)
for (i in 1:nspec){
  hist(post$eta_species[, i], breaks=40, 
       main=paste("Species intercept", levels(pd$speciesname)[i]),
       xlab="Value", ylab="Posterior density")
  hdi_local <- HDI(post$beta_spec[, i, 1])
  abline(v=hdi_local, lty=2, col="red")
  abline(v=0, lty=2, lwd=2)
}

# plot responses to local richness
par(mfrow=c(2, 3), mar=c(5, 4, 4, 2) + 0.1)
for (i in 1:nspec){
  hist(post$beta_spec[, i, 2], breaks=40, 
       main=paste("Local richness effect", levels(pd1$HOSTSPECIES)[i]),       xlab="Value", ylab="Posterior density")
  hdi_local <- HDI(post$beta_spec[, i, 2])
  abline(v=hdi_local, lty=2, col="red")
  abline(v=0, lty=2, lwd=2)
}


# plot responses to regional richness
par(mfrow=c(2, 3), mar=c(5, 4, 4, 2) + 0.1)
for (i in 1:nspec){
  hist(post$beta_spec[, i, 3], breaks=40, 
       main=paste("Regional richness effect", levels(pd1$HOSTSPECIES)[i]),
       xlab="Value", ylab="Posterior density")
  hdi_local <- HDI(post$beta_spec[, i, 3])
  abline(v=hdi_local, lty=2, col="red")
  abline(v=0, lty=2, lwd=2)
}

waic_fit <- waic(m_fit)

str(post)
med_mu <- exp(apply(post$log_mu, 2, median))
plot(med_mu, stan_d$y)
abline(0, 1, lty=2)

plot(stan_d$y, stan_d$y - med_mu)
abline(h=0, lty=2)


# variance decomposition
par(mfrow=c(2, 3))
br <- seq(0, 13, length.out=40)
xvals <- seq(0, 13, length.out=100)
hist(post$sigma_region, breaks=br, main="Among region sd", freq=F)
lines(x=xvals, y=dcauchy(xvals, 0, 3))
hist(post$sigma_site, breaks=br, main="Among site sd", freq=F)
lines(x=xvals, y=dcauchy(xvals, 0, 3))
hist(post$sigma_spec[, 1], breaks=br, main="Among host species sd", freq=F)
lines(x=xvals, y=dcauchy(xvals, 0, 3))
hist(post$sigma_spec[, 2], breaks=br, main="Among host species sd", freq=F)
lines(x=xvals, y=dcauchy(xvals, 0, 3))
hist(post$sigma_spec[, 3], breaks=br, main="Among host species sd", freq=F)
lines(x=xvals, y=dcauchy(xvals, 0, 3))
hist(post$sigma_indiv, breaks=br, main="Among individual sd", freq=F)
lines(x=xvals, y=dcauchy(xvals, 0, 3))
par(mfrow=c(1, 1))






