# fitting a simple random effect model to the rib abundance data
source('R/d_clean_joint.R')
source('R/helpers.R')
library(rstan)

stan_d <- list(n = nrow(pd1), 
               nsite = length(unique(pd1$SiteCode)), 
               nregion = max(coord_d$nbhood), 
               site = pd1$numsite, 
               region = region, 
               y = pd1$RIB, 
               nspec = length(unique(pd1$HOSTSPECIES)),
               h_spec = as.numeric(pd1$HOSTSPECIES),
               local_richness = host_data$LARVAMPRICHNOTAGR, 
               region_richness = reg_rich)

watch <- c("sigma_indiv", "sigma_site", "sigma_region", "alpha_region", 
           "alpha_site", "alpha_0", "beta_local", "beta_region", 
           "log_mu", "log_lik", "alpha_spec")

m_init <- stan('R/mod.stan', data=stan_d, chains=1, iter=1, pars=watch)
source('http://mc-stan.org/rstan/stan.R')
m_fit <- stan(fit=m_init, data=stan_d, pars=watch, chains=3, 
              iter=1200, warmup=500)
rm(stan)
library(rstan)

traceplot(m_fit, 
          pars = c("sigma_indiv", "sigma_site", "sigma_region", "alpha_region", 
                   "alpha_0", "beta_local", "beta_region", "alpha_spec"),
          inc_warmup=T)


post <- rstan::extract(m_fit)
par(mfrow=c(1, 2), mar=c(5, 4, 4, 2) + 0.1)
br <- seq(-10, 10, .2)
hist(post$beta_local, breaks=br, main="Local richness effect", 
     xlab="Value", ylab="Posterior density")
hdi_local <- HDI(post$beta_local)
abline(v=hdi_local, lty=2, col="red")

hist(post$beta_region, breaks=br, main="Regional richness effect", 
     xlab="Value", ylab="Posterior density")
hdi_region <- HDI(post$beta_region)
abline(v=hdi_region, lty=2, col="red")
par(mfrow=c(1, 1))


waic_fit <- waic(m_fit)

str(post)
med_mu <- exp(apply(post$log_mu, 2, median))
plot(med_mu, stan_d$y)
abline(0, 1, lty=2)

plot(stan_d$y, stan_d$y - med_mu)
abline(h=0, lty=2)

