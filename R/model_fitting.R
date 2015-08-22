# fitting a simple random effect model to the rib abundance data
source('R/d_clean.R')
source('R/helpers.R')
library(rstan)
library(maps)
library(scales)

## Define spatial neighbors -----------------------
agg_dist <- 5
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

# Bundle data for stan ---------------------------------------------
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
               site_snail = as.numeric(snail_inf_df$fsite))

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
           "sigma_p_region")

m_init <- stan('R/isd.stan', data=stan_d, chains=1, iter=1, pars=watch)
m_fit <- stan(fit=m_init, data=stan_d, pars=watch, chains=3, 
              iter=500, cores=3)

traceplot(m_fit, pars=watch, inc_warmup=T)
waic(m_fit)

library(ggmcmc)
ggd <- ggs(m_fit)
ggs_caterpillar(ggd, 'eta_year')
ggs_caterpillar(ggd, 'eta_site')
ggs_caterpillar(ggd, 'p_site')
ggs_caterpillar(ggd, 'p_region')

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






