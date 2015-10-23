# Cleaning the parasite infection data, and
# defining spatial neighbors
source('R/d_clean.R')
library(effects)
# d2 is km distance
nseq <- 10
dseq <- seq(.2, 1.5, length.out=nseq)
logistic_res <- array(dim=c(nseq, 3, 3))

pois_res <- array(dim=c(nseq, 4, 3))
rho_isd <- rep(NA, nseq)
effplots <- list()
for (s in 1:nseq){
  nbs <- dnearneigh(s_coord, d1=0, d2=dseq[s], longlat=T)
  par(mfrow=c(1, 1))
  #plot(nbs, coords=s_coord)
  (neighborhoods <- n.comp.nb(nbs))
  num_neigh <- neighborhoods$nc
  coord_d$nbhood <- neighborhoods$comp.id
  #plot(coord_d$Lon, coord_d$Lat, col=coord_d$nbhood)
  
  all(levels(pd$fsite) == levels(coord_d$site))
  pd$numsite <- as.numeric(pd$fsite)
  pd$region <- coord_d$nbhood[pd$numsite]
  
  region <- coord_d$nbhood[order(coord_d$site)]
  
  # calculate regional host richness based on regional definition
  host_cols <- names(snails)[c(57:63)]
  nyear <- length(unique(pd$year))
  
  all(levels(snails$fsite) == levels(pd$fsite))
  snails$num_site <- as.numeric(snails$fsite)
  snails$region <- region[snails$num_site]
  
  # calculate regional richness
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
  
  snails$reg_year <- paste(snails$region, snails$AssmtYear_1, sep='.')
  
  # add regional richness data to pd
  pd$helisoma_density <- snails$`HELI_ SW_1`[match(pd$fsite, snails$fsite)]
  pd$heli_rib_prevalence <- snails$HELIRIB_3[match(pd$fsite, snails$fsite)]
  pd$ISD <- pd$helisoma_density * pd$heli_rib_prevalence
  #plot(jitter(snails$reg_rich), jitter(snails$local_rich))
  pd$reg_year <- paste(pd$region, pd$year, sep='.')
  pd$reg_rich <- snails$reg_rich[match(pd$reg_year, snails$reg_year)]
  #plot(jitter(pd$reg_rich), jitter(pd$local_rich))

  pd <- pd[!is.na(pd$reg_rich), ]
  
  # is Rib present in all regions? 
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
  cc <- complete.cases(summary_d[c('mean_ISD', 'sd_ISD')])
  rho_isd[s] <- cor(summary_d$mean_ISD[cc], summary_d$sd_ISD[cc])
  summary_d$reg_year <- paste(summary_d$region, summary_d$year, sep='.')
  pd$mean_ISD <- summary_d$mean_ISD[match(pd$reg_year, summary_d$reg_year)]
  # only want to fit model to data at sites where rib occurs
  rib_occ <- ddply(pd, c('fsite', 'year'), summarize, 
                   rib_occ=as.numeric(any(rib > 0) | any(ISD > 0)))
  
  # fit logistic regression for presence/absence of rib in region
  if (sd(summary_d$rib_pres) != 0){
    mod <- glm(rib_pres ~ reg_rich + nsampled, family='binomial', 
               data=summary_d)
    logistic_res[s, , 1] <- coef(mod)
    if (any(summary(mod)$coefficients[, 'z value'] == 0)){
      print('glm coefs cannot generate CIs')
    }else{
      try(logistic_res[s, , 2:3] <- confint(mod))
    }
  }
  
  # fit abundance model to look at intxn between richness and ISD
  summary_d$log_rib <- log(summary_d$mean_rib)
  #pd$rib_pres <- summary_d$rib_pres[match(pd$reg_year, summary_d$reg_year)]
  

  rib_occ$fsiteyear <- paste(rib_occ$fsite, rib_occ$year, sep='.')
  pd$fsiteyear <- paste(pd$fsite, pd$year, sep='.')
  pd$rib_occ <- rib_occ$rib_occ[match(pd$fsiteyear, rib_occ$fsiteyear)]

  
  mod2 <- glmer(rib ~ reg_rich*mean_ISD + #year +
                  (1|fsite),# + (1|Dissection), 
             data=subset(pd, rib_occ==1), family='poisson')
  effplots[[s]] <- allEffects(mod2)
  pois_res[s, , 1] <- fixef(mod2)
  pois_res[s, , 2:3] <- confint(mod2, method='Wald')[-1, ]
  
  if (s == nseq){
    colnames(logistic_res) <- names(coef(mod))
    colnames(pois_res) <- names(fixef(mod2))
  }
}

plot_res <- function(res, param, dseq, ...){
  plot(dseq, res[, param, 1], type='b', 
       ylim=range(res[, param, ], na.rm=TRUE), 
       xlab='Aggregation distance', ...)
  lines(dseq, res[, param, 2], lty=2)
  lines(dseq, res[, param, 3], lty=2)
  abline(h=0, lty=3)
}

dimnames(logistic_res)
par(mfrow=c(1, 2))
plot_res(logistic_res[1:9, , ], 'reg_rich', dseq[1:9], 
         ylab='Effect of regional richness on occurrence')
plot_res(logistic_res[1:9, , ], 'nsampled', dseq[1:9], 
         ylab='Effect of number of sites in region on occurrence')

par(mfrow=c(1, 3))
#plot_res(pois_res, 'reg_rich', dseq, ylab='Effect of regional richness')
#plot_res(pois_res, 'mean_ISD', dseq, ylab='Effect of mean infected snail density')
#plot_res(pois_res, "reg_rich:mean_ISD", dseq, ylab='Intxn: regional richness and mean ISD')

par(mfrow=c(1, 1))
par(mar=c(5, 4.4, 4, 2) + 0.1)
#plot(dseq, rho_isd, xlab="Distance of aggregation", 
#     ylab=expression(
#       paste("Correlation between ", hat(mu[ISD]), " and ", hat(sigma[ISD]))
#       )
#     )

mod3 <- glmer(rib ~ local_rich * ISD + 
                (1|fsite) + (1|year), 
              data=subset(pd, rib_occ==1), family='poisson')


mod4 <- glmer(rib ~ local_rich + ISD + 
                (1|fsite) + (1|year) + (1|Dissection), 
              data=subset(pd, rib_occ==1), family='poisson')

#summary(mod3)
