# Prepping data for the non-aggregated single-year analysis
source('R/single_year/d_clean_joint.R')
source('R/helpers.R')

rib_pres <- FALSE
if (rib_pres){
  keep <- keep[!keep %in% unique(rib_abs$SiteCode)]
  gd <- droplevels(subset(gd, site %in% keep))
  snail_infection <- droplevels(subset(snail_infection, site %in% keep))
  snail_infection$fsite <- factor(snail_infection$site)
  pd1 <- droplevels(subset(pd1, SiteCode %in% keep))
  coord_d <- droplevels(subset(coord_d, sitename %in% keep))
  snail_den <- droplevels(subset(snail_den, Site %in% keep))
  
  library(reshape2)
  msnail <- melt(snail_infection, 
                 id.vars = c('fsite', 'DissectCode'), 
                 measure.vars = c('ribprez'))
  
  snail_agg <- ddply(msnail, ~ fsite, 
                     summarize,
                     n_infected = sum(value), 
                     n_sampled = length(value), 
                     rib_in_snails = any(value > 0))
  
  # calculate observed local richness
  site_sum <- ddply(subset(gd, Species != 'HELI'), 
                    ~ site, summarize, 
                    local_rich = length(unique(Species[Count > 0])))
  
  all.equal(levels(factor(gd$site)), levels(coord_d$fsite))
  all.equal(levels(factor(pd1$SiteCode)), levels(coord_d$fsite))
  all.equal(levels(factor(snail_infection$site)), levels(coord_d$fsite))
  all.equal(levels(factor(pd1$SiteCode)), levels(coord_d$fsite))
  all.equal(levels(factor(pd1$SiteCode)), levels(snail_den$Site))
}


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
              chains=3, cores=3, iter=2000, data=stan_d)

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
                 isd = seq(min(isdm), max(isdm), length.out=20))
d$mu <- NA
d$hdi_lo <- NA
d$hdi_hi <- NA
d$q25 <- NA
d$q75 <- NA

mu_pred <- array(dim=c(nrow(d), length(post$mu_a)))
for (i in 1:nrow(d)){
  mu_pred[i, ] <- post$mu_a + post$beta_rich * d$richness[i] + 
    post$beta_isd * d$isd[i] + 
    post$beta_intxn * d$isd[i] *d$richness[i]
  hdi <- HDI(mu_pred[i, ])
  d$hdi_lo[i] <- hdi[1]
  d$hdi_hi[i] <- hdi[2]
  d$q25[i] <- quantile(mu_pred[i, ], probs = .25)
  d$q75[i] <- quantile(mu_pred[i, ], probs = .75)
  d$mu[i] <- median(mu_pred[i, ])
}

site_sum$isd <- isdm
site_sum$richness <- site_sum$local_rich
site_sum$logisd <- log(site_sum$isd)

isd_d <- melt(post$isd, varnames=c('iter', 'site'), value.name='isd')
isd_d$richness <- site_sum$local_rich[match(isd_d$site, 
                                            as.numeric(factor(site_sum$site)))]
isd_d$logisd <- log(isd_d$isd)

d$logisd <- log(d$isd)

d$panel_lab <- paste(d$richness, 'larval amphibian species')
site_sum$panel_lab <- paste(site_sum$local_rich, 'larval amphibian species')

ggplot(d, aes(x=isd)) + 
  facet_wrap(~panel_lab) + 
  geom_ribbon(aes(ymin = hdi_lo, ymax=hdi_hi), alpha=.5) + 
  geom_ribbon(aes(ymin = q25, ymax=q75), alpha=.5) + 
  geom_line(aes(y=mu, group=richness), size=2) +
  geom_rug(data=site_sum) + 
  theme_bw()

ggplot(isd_d, aes(x=logisd, group=site)) + geom_density() +
  facet_wrap(~richness) 
  
