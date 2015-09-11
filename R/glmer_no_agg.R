# GLMM to evaluate relationship bt ISD and richness on PSRE infection
rm(list=ls())
library(readxl)
library(ggplot2)
library(dplyr)
library(lme4)

# load data
pd <- read.csv("data/necropsy20122014_all_rib.csv")
pd <- subset(pd, speciescode == "PSRE")
pd <- pd[complete.cases(pd),]
pd <- droplevels(pd)
snails <- read_excel('data/CA_Master_20122014v7_snails.xlsx')
snails <- snails[complete.cases(snails[c('HELIRIB_3', "HELI_ SW_1")]), ]
snails <- droplevels(snails)
snails$local_rich <- NA
snails$amca <- NA
snails$bubo <- NA
snails$psre <- NA
snails$raca <- NA
snails$radr <- NA
snails$tagr <- NA
snails$tato <- NA

for (i in 1:nrow(snails)){
  snails$amca[i] <- any(c(snails$AMCA_SW_1[i], snails$AMCA_SE_1[i]) > 0)
  snails$bubo[i] <- any(c(snails$BUBO_SW_1[i], snails$BUBO_SE_1[i]) > 0)
  snails$psre[i] <- any(c(snails$PSRE_SW_1[i], snails$PSRE_SE_1[i]) > 0)
  snails$raca[i] <- any(c(snails$RACA_SW_1[i], snails$RACA_SE_1[i]) > 0)
  snails$radr[i] <- any(c(snails$RADR_SW_1[i], snails$RADR_SE_1[i]) > 0)
  snails$tagr[i] <- any(c(snails$TAGR_SW_1[i], snails$TAGR_SE_1[i]) > 0)
  snails$tato[i] <- any(c(snails$TATO_SW_1[i], snails$TATO_SE_1[i]) > 0)
}

snails$local_rich <- apply(snails[, c('amca', 'bubo', 'psre', 'raca', 
                                      'radr', 'tagr', 'tato')], 1, sum)
snails <- snails[!(snails$local_rich==0), ] # remove sites with richness = 0

names(snails)
names(pd)

levels(pd$siteyear)[levels(pd$siteyear) == 'HIDDEN_2014'] <- 'Hidden_2014'
levels(pd$siteyear)[levels(pd$siteyear) == 'FROG_2014'] <- 'Frog_2014'

pd$local_rich <- snails$local_rich[match(pd$siteyear, snails$siteyear)]
pd$site <- as.character(pd$sitecode)

# fix some site names
pd$sitename <- as.character(pd$sitename)
pd$sitename[pd$sitename == 'GRANT LAKE'] <- 'Grant Lake'
pd$sitename[pd$sitename == 'MCCREERY LAKE'] <- 'McCreery Lake'
snails$SiteName[snails$SiteName == 'Mccreery Lake'] <- "McCreery Lake"
pd$sitename[pd$sitename == 'SF25'] <- 'CA-SF25'
pd$sitename[pd$sitename == 'SF85a'] <- 'CA-SF85a'
pd$sitename[pd$sitename == 'Ca-SF31'] <- 'CA-SF31'
pd$sitename[pd$sitename == 'CA-SF40'] <- 'SF-40'
pd$sitename[pd$sitename == 'CA-SF41'] <- 'SF-41'
pd$sitename[pd$sitename == 'MITTENS'] <- 'Mittens'
pd$sitename[pd$sitename == 'OWLHENGE'] <- 'Owlhenge'
pd$sitename[pd$sitename == 'BARN'] <- 'Barn'
snails$SiteName[snails$SiteName == 'BARN'] <- 'Barn'
pd$sitename[pd$sitename == 'CA-SF85a'] <- 'CA-SF85'

setdiff(pd$sitename, snails$SiteName)
intersect(pd$sitename, snails$SiteName)

pd <- subset(pd, sitename %in% snails$SiteName)
pd <- droplevels(pd)
snails <- subset(snails, SiteName %in% pd$sitename)
snails <- droplevels(snails)

pd$heli_den <- snails$`HELI_ SW_1`[match(pd$siteyear, snails$siteyear)]
pd$heli_prev <- snails$HELIRIB_3[match(pd$siteyear, snails$siteyear)]
pd$isd <- pd$heli_den * pd$heli_prev

pd <- pd[!is.na(pd$local_rich), ]

ggplot(pd, aes(x=isd, y=rib)) + 
  facet_wrap(~local_rich) +
  geom_jitter(shape=1) + 
  xlab('Infected snail density') + 
  ylab('Rib infections')

mod <- glm(rib ~ isd * local_rich, family='poisson', data=pd)
library(effects)
plot(allEffects(mod))

library(lme4)
mod2 <- glmer(rib ~ isd * local_rich + 
                (1|sitename) + (1|year) + (1|Dissection), 
              family='poisson', data=pd)
plot(allEffects(mod2))

# maybe this is driven by outliers in infected snail density
pdsub <- subset(pd, isd < 1)
mod3 <- glm(rib ~ isd * local_rich, family='poisson', data=pdsub)
plot(allEffects(mod3))

mod4 <- glmer(rib ~ isd * local_rich + 
                (1|sitename) + (1|year) + (1|Dissection), 
              family='poisson', data=pdsub)
plot(allEffects(mod4))

# maybe it's a 2013 thing
mod5 <- glmer(rib ~ isd * local_rich + 
                (1|sitename) + (1|Dissection), 
              family='poisson', data=subset(pd, year==2013))
# comes up as not significant
