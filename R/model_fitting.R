# fitting a simple random effect model to the rib abundance data
source('R/d_clean.R')
source('R/helpers.R')
library(rstan)
library(maps)
library(scales)

## Define spatial neighbors -----------------------
agg_dist <- 2
nbs <- dnearneigh(s_coord, d1 = 0, d2 = agg_dist, longlat = T)
par(mfrow=c(1, 1))
plot(nbs, coords=s_coord)
(neighborhoods <- n.comp.nb(nbs))
num_neigh <- neighborhoods$nc
coord_d$nbhood <- neighborhoods$comp.id
plot(coord_d$Lon, coord_d$Lat, col=coord_d$nbhood)
map('county', add=T, col=alpha(1, .3))
map('state', add=T)

all(levels(pd$fsite) == levels(coord_d$site))
pd$numsite <- as.numeric(pd$fsite)
pd$region <- coord_d$nbhood[pd$numsite]
region <- coord_d$nbhood[order(coord_d$site)]

# Calculate apparent regional host richness ----------------------------
host_cols <- c("amca", "bubo", "psre", "raca", "radr", "tagr", "tato")
nyear <- length(unique(pd$year))

all(levels(snails$fsite) == levels(pd$fsite))
snails$num_site <- as.numeric(snails$fsite)
snails$region <- region[snails$num_site]

snails$reg_rich <- NA
for (i in 1:max(snails$region)){
  assmtyears <- unique(snails$AssmtYear_1[snails$region == i])
  if (length(assmtyears) == 0){
    print(paste('No assessment years for region', i))
    next
  }
  for (j in 1:length(assmtyears)){
    indx <- which(snails$region == i & 
                    snails$AssmtYear_1 == assmtyears[j])
    if (length(indx) == 0){
      print(paste('length of indx = 0 for i, j', i, j, sep=' '))
    }
    subd <- snails[indx, host_cols]
    seen_at_all <- apply(subd, 2, FUN=function(x) sum(x) > 0)
    snails$reg_rich[indx] <- sum(seen_at_all)
  }
}

# Add regional richness data from snails to pd ----------------------
snails$reg_year <- paste(snails$region, snails$AssmtYear_1, sep='.')
pd$helisoma_density <- snails$`HELI_ SW_1`[match(pd$fsite, snails$fsite)]
pd$heli_rib_prevalence <- snails$HELIRIB_3[match(pd$fsite, snails$fsite)]
pd$ISD <- pd$helisoma_density * pd$heli_rib_prevalence
plot(jitter(snails$reg_rich), jitter(snails$local_rich))
pd$reg_year <- paste(pd$region, pd$year, sep='.')
pd$reg_rich <- snails$reg_rich[match(pd$reg_year, snails$reg_year)]
plot(jitter(pd$reg_rich), jitter(pd$local_rich))

pd <- pd[!is.na(pd$reg_rich), ]

# Summarize Rib data by region & year ------------------------------
library(plyr)
summary_d <- ddply(pd, c("region", 'year'), 
                   summarise, 
                   rib_pres = as.numeric(any(rib > 0) | any(ISD > 0)),
                   reg_rich = unique(reg_rich),
                   mean_ISD = mean(ISD),
                   sd_ISD = sd(ISD),
                   nsampled = length(site), 
                   mean_rib = mean(rib))
summary_d
plot(summary_d$nsampled, summary_d$reg_rich)
cc <- complete.cases(summary_d[c('mean_ISD', 'sd_ISD')])
cor(summary_d$mean_ISD[cc], summary_d$sd_ISD[cc])
summary_d$reg_year <- paste(summary_d$region, summary_d$year, sep='.')
pd$mean_ISD <- summary_d$mean_ISD[match(pd$reg_year, summary_d$reg_year)]

# Subset data to locations where Rib detected --------------------------
rib_occ <- ddply(pd, c('fsite', 'year'), summarize, 
                 rib_occ=as.numeric(any(rib > 0) | any(ISD > 0)))
rib_occ$fsiteyear <- paste(rib_occ$fsite, rib_occ$year, sep='.')
pd$fsiteyear <- paste(pd$fsite, pd$year, sep='.')
pd$rib_occ <- rib_occ$rib_occ[match(pd$fsiteyear, rib_occ$fsiteyear)]

pd$fsiteyear <- as.factor(pd$fsiteyear)
snails$fsiteyear <- factor(paste(snails$fsite, snails$AssmtYear_1, sep='.'))
all(levels(pd$fsiteyear) == levels(snails$fsiteyear))
all(levels(snails$fsiteyear) == levels(factor(rib_occ$fsiteyear)))
snails$rib_occ <- rib_occ$rib_occ[match(snails$fsiteyear, rib_occ$fsiteyear)]
pd <- subset(pd, rib_occ == 1)
snails <- subset(snails, rib_occ == 1)
pd <- droplevels(pd)
snails <- droplevels(snails)

# Extract snail density data ---------------------------------------------
# treat helisoma density in sweeps and seines as a multimethod indicator

names(snails)
with(snails, {
  plot(log(1+`HELI_ SW_1`), log(1+`HELI_ SE_1`))
})

library(reshape2)
snail_density_df <- melt(snails, 
                         id.vars=c('AssmtCode_1', 'fsite', 'AssmtYear_1', 
                                   'region'), 
                         measure.vars = c('HELI_ SW_1', 'HELI_ SE_1'))

snail_density_df <- snails[!duplicated(snails$AssmtCode_1), 
                           c('fsite', 'AssmtYear_1', 
                             'region', 'HELI_ SW_1', 'HELI_ SE_1')]
snail_density_df$heli_sweep <- c(scale(log(1 + snail_density_df$`HELI_ SW_1`)))
snail_density_df$heli_seine <- c(scale(log(1 + snail_density_df$`HELI_ SE_1`)))

with(snail_density_df, {
  plot(heli_seine, heli_sweep)
})

snails$heli_infected <- round(snails$HELDISS3 * snails$HELIRIB_3)

# Extract snail infection data
snail_inf_df <- snails[!duplicated(snails$AssmtCode_1), 
                         c('fsite', 'AssmtYear_1', 'region', 
                           "HELDISS3", "heli_infected")]

# final verification that factor levels are comparable across datasets
all(levels(pd$fsite) == levels(snails$fsite))
all(levels(snail_density_df$fsite) == levels(snail_inf_df$fsite))
all(levels(pd$fsite) == levels(snail_inf_df$fsite))

# create vector of length nsite that describes regional identity
reg <- pd[!duplicated(pd$fsite), 'region']

# Bundle data for stan ---------------------------------
stan_d <- list(n = nrow(pd), 
               nsite = length(levels(pd$fsite)), 
               nregion = max(pd$region), 
               site = as.numeric(pd$fsite), 
               nyear = length(unique(pd$year)),
               year = pd$year - min(pd$year) + 1,
               region = reg, 
               y = pd$rib, 
               richness = pd$reg_rich, 
               nsnail = nrow(snail_inf_df),
               y_snail = snail_inf_df$heli_infected, 
               k_snail = snail_inf_df$HELDISS3, 
               site_snail = as.numeric(snail_inf_df$fsite), 
               n_snail_density = nrow(snail_density_df), 
               snail_density1 = snail_density_df$heli_sweep, 
               snail_density2 = snail_density_df$heli_seine, 
               snail_den_site = as.numeric(snail_density_df$fsite))

watch <- c("eta_site", 
           "eta_year",
           "a0", 
           "sigma_site", 
           "sigma_indiv",
           "beta_rich", 
           "mu_p",
           "p_site",
           "sigma_p", 
           "p_region", 
           "sigma_p_region", 
           'd_site', 
           'sigma_snail_density', 
           'beta_y', 
           'd_region', 
           'sigma_d_region', 
           'beta_isd', 
           'beta_intxn', 
           'mu_den')

m_init <- stan('R/isd.stan', data=stan_d, chains=1, iter=1, pars=watch)
m_fit <- stan(fit=m_init, data=stan_d, pars=watch, chains=3, 
              iter=1000, cores=3, warmup = 400)

traceplot(m_fit, pars=watch, inc_warmup=T)
waic(m_fit)

library(ggmcmc)
ggd <- ggs(m_fit)
ggs_caterpillar(ggd, 'eta_year')
ggs_caterpillar(ggd, 'eta_site')
ggs_caterpillar(ggd, 'p_site')
ggs_caterpillar(ggd, 'p_region')
ggs_caterpillar(ggd, 'd_site')
ggs_caterpillar(ggd, 'd_region')

post <- rstan::extract(m_fit)

par(mfrow=c(3, 1))
v <- unlist(post[c('beta_rich', 'beta_isd', 'beta_intxn')])
br <- seq(min(v), max(v), length.out=60)
hist(post$beta_rich, breaks=br, 
     main=paste('Fixed effs: agg. dist. =', agg_dist, "km"))
abline(v=HDI(post$beta_rich), lty=2)
hist(post$beta_isd, breaks=br, main="")
abline(v=HDI(post$beta_isd), lty=2)
hist(post$beta_intxn, breaks=br, main="")
abline(v=HDI(post$beta_intxn), lty=2)



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






