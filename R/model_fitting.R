# fitting a simple random effect model to the rib abundance data
library(rstan)
library(ggmcmc)
stan_d <- list(n = nrow(pd1), 
               nsite = length(unique(pd1$SiteCode)), 
               nregion = max(coord_d$nbhood), 
               site = pd1$numsite, 
               region = region, 
               y = pd1$RIB, 
               local_richness = host_data$LARVAMPRICHNOTAGR, 
               region_richness = reg_rich)

watch <- c("sigma_indiv", "sigma_site", "sigma_region", "alpha_region", 
           "alpha_site", "alpha_0", "phi", "beta_local", "beta_region")

m_init <- stan('R/mod.stan', data=stan_d, chains=1, iter=1, pars=watch)
m_fit <- stan(fit=m_init, data=stan_d, pars=watch, chains=3, iter=1000)

traceplot(m_fit)


ggd <- ggs(m_fit)
ggs_caterpillar(ggd, "alpha_site")
