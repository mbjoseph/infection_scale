# Prepping data for the non-aggregated single-year analysis
source('R/single_year/d_clean_joint.R')
stan_d <- list(
               # amphibian Rib infection data
               n_a = nrow(pd1), 
               nsite = length(levels(pd1$SiteCode)), 
               y_a = pd1$RIB, 
               site_a = as.numeric(pd1$SiteCode), 
               rich = site_sum$local_rich,
               # snail Rib infection data
               n_s = nrow(snail_agg), 
               site_s = as.numeric(snail_agg$fsite), 
               inf_s = snail_agg$n_infected, 
               samp_s = snail_agg$n_sampled, 
               # snail density data
               n_sd = nrow(snail_den), 
               y_sd = snail_den$count, 
               site_sd = as.numeric(snail_den$Site))

library(rstan)
watch <- c('mu_a', 'sd_a', 'mu_s', 'sd_s', 'mu_sd', 'sd_sd', 
           'alpha_a','alpha_s', 'alpha_sd', 
           'isd', 'beta_isd', 'beta_rich', 'beta_intxn', 
           'phi_a', 'phi_sd')

m_init <- stan('R/single_year/no_agg.stan', 
               chains=1, iter=10, data=stan_d)
m_fit <- stan(fit=m_init, pars=watch,
              chains=2, cores=2, iter=2000, data=stan_d)

rstan::traceplot(m_fit, inc_warmup=F)

library(ggmcmc)
ggd <- ggs(m_fit)
ggs_caterpillar(ggd, 'alpha_s\\.')
ggs_caterpillar(ggd, 'alpha_sd')
ggs_caterpillar(ggd, 'isd\\.')
ggs_caterpillar(ggd, 'alpha_a\\.')

# plot response surface
post <- rstan::extract(m_fit)
isdm <- apply(post$isd, 2, median)
rich_vals <- min(stan_d$rich):max(stan_d$rich)
d <- expand.grid(richness = rich_vals, 
                 isd = seq(min(isdm), max(isdm), length.out=40))
d$mu <- NA
d$hdi_lo <- NA
d$hdi_hi <- NA

mu_pred <- array(dim=c(nrow(d), length(post$mu_a)))
for (i in 1:nrow(d)){
  for (j in 1:length(post$mu_a)){
    mu_pred[i, j] <- post$mu_a[j] + post$beta_rich[j] * d$richness[i] + 
                      post$beta_isd[j] * d$isd[i] + 
                      post$beta_intxn[j] * d$isd[i] *d$richness[i]
  }
  hdi <- HDI(mu_pred[i, ])
  d$hdi_lo[i] <- hdi[1]
  d$hdi_hi[i] <- hdi[2]
  d$mu[i] <- median(mu_pred[i, ])
}

site_sum$isd <- isdm
site_sum$richness <- site_sum$local_rich
site_sum$logisd <- log(site_sum$isd)

isd_d <- melt(post$isd, varnames=c('iter', 'site'), value.name='isd')
isd_d$richness <- site_sum$local_rich[match(isd_d$site, as.numeric(factor(site_sum$site)))]
isd_d$logisd <- log(isd_d$isd)

d$logisd <- log(d$isd)

ggplot(d, aes(x=isd)) + 
  facet_wrap(~richness) + 
  geom_ribbon(aes(ymin = hdi_lo, ymax=hdi_hi)) + 
  geom_line(aes(y=mu), color='red') +
  geom_rug(data=site_sum)

ggplot(isd_d, aes(x=logisd, group=site)) + geom_density() +
  facet_wrap(~richness) 
  
