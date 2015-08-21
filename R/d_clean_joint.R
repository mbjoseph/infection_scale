# Cleaning the parasite infection data, and
# defining spatial neighbors

rm(list=ls())
library(readxl)
library(ggplot2)
library(dplyr)

# load data
pd <- read.csv("data/necropsy20122014_all_rib.csv")
pd <- pd[complete.cases(pd),]
pd <- droplevels(pd)

# load snail infection data
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

# insert siteXyear richness values into pd data.frame
pd$local_rich <- NA
for (i in 1:length(levels(pd$siteyear))){
  siteyear <- levels(pd$siteyear)[i]
  indx1 <- which(pd$siteyear == siteyear)
  indx2 <- which(snails$siteyear == siteyear)
  if (length(indx1) == 0) print(paste("siteyear not in pd: ", siteyear))
  if (length(indx2) == 0){
    print(paste("siteyear not in snails: ", siteyear))
    indx2 <- grep(siteyear, snails$siteyear,
                  ignore.case=TRUE, value=TRUE)
    indx2 <- which(snails$siteyear == indx2)
    print(paste("new indx2:", indx2))
    if (length(indx2) == 0){
      print("Removing site-year from data, because we lack sweep/seine data")
      pd <- pd[-indx1, ] # remove siteyear from parasite dissection data
      next # move on
    }
  }
  pd$local_rich[indx1] <- rep(snails$local_rich[indx2], length(indx1))
  if (any(is.na(pd$local_rich[indx1]))) print(paste(siteyear, "has NA"))
  if (snails$local_rich[indx2] == 0) print(paste(siteyear, "has 0 richness"))
}

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

setdiff(pd$sitename, snails$SiteName)
intersect(pd$sitename, snails$SiteName)

pd <- subset(pd, sitename %in% snails$SiteName)
pd <- droplevels(pd)
snails <- subset(snails, SiteName %in% pd$sitename)
snails <- droplevels(snails)


# classify sites into regions based on spatial neighborhoods
snails$fsite <- as.factor(snails$SiteName)
pd$fsite <- as.factor(pd$sitename)
cbind(levels(pd$fsite), levels(snails$fsite))

coord_d <- subset(snails, !duplicated(snails$fsite), 
                  c('fsite', 'Longitude_1', 'Latitude_1'))
coord_d <- coord_d[order(coord_d$fsite), ]
names(coord_d) <- c('site', 'Lon', 'Lat')
plot(coord_d$Lon, coord_d$Lat)

# define neighborhoods
library(spdep)
class(coord_d) <- 'data.frame'
s_coord <- coordinates(coord_d[, c("Lon", "Lat")])

# d2 is km distance
nseq <- 10
dseq <- seq(0, 2.2, length.out=nseq)
logistic_res <- array(dim=c(nseq, 5, 3))

for (s in 1:nseq){
  nbs <- dnearneigh(s_coord, d1=0, d2=dseq[s], longlat=T)
  par(mfrow=c(1, 1))
  plot(nbs, coords=s_coord)
  (neighborhoods <- n.comp.nb(nbs))
  num_neigh <- neighborhoods$nc
  coord_d$nbhood <- neighborhoods$comp.id
  plot(coord_d$Lon, coord_d$Lat, col=coord_d$nbhood)
  
  all(levels(pd$fsite) == levels(coord_d$site))
  pd$numsite <- as.numeric(pd$fsite)
  pd$region <- coord_d$nbhood[pd$numsite]
  
  region <- coord_d$nbhood[order(coord_d$site)]
  
  # calculate regional host richness based on regional definition
  host_cols <- names(snails)[c(57:63)]
  nyear <- length(unique(pd$year))
  
  # assume regional richness is constant across years, otherwise run into issues 
  # with uneven sampling, and lack of detection of species
  reg_rich <- array(NA, dim=c(num_neigh))
  n_sampled <- array(dim=c(num_neigh))
  n_in_region <- rep(NA, num_neigh)
  
  all(levels(snails$fsite) == levels(pd$fsite))
  snails$num_site <- as.numeric(snails$fsite)
  snails$region <- region[snails$num_site]
  
  # calculate regional richness
  snails$reg_rich <- NA
  for (i in 1:max(snails$region)){
    indx <- which(snails$region == i)
    subd <- snails[indx, host_cols]
    seen_at_all <- apply(subd, 2, FUN=function(x) sum(x) > 0)
    snails$reg_rich[indx] <- sum(seen_at_all)
  }

  pd$reg_rich <- snails$reg_rich[match(pd$region, snails$region)]
  
  # add regional richness data to pd
  pd$helisoma_density <- snails$`HELI_ SW_1`[match(pd$fsite, snails$fsite)]
  pd$heli_rib_prevalence <- snails$HELIRIB_3[match(pd$fsite, snails$fsite)]
  pd$ISD <- pd$helisoma_density * pd$heli_rib_prevalence
  plot(jitter(pd$reg_rich), jitter(pd$local_rich))
  str(pd)
  
  # is Rib present in all regions? 
  library(plyr)
  summary_d <- ddply(pd, c("region", "year"), 
        summarise, 
        rib_pres = as.numeric(any(rib > 0)),
        reg_rich = unique(reg_rich),
        mean_ISD = mean(ISD),
        nsampled = length(site))
  summary_d
  
  # fit logistic regression for presence/absence of rib in region
  mod <- glm(rib_pres ~ reg_rich * mean_ISD + year, family='binomial', 
             data=summary_d)
  logistic_res[s, , 1] <- coef(mod)
  try(logistic_res[s, , 2:3] <- confint(mod))
  if (s == nseq){
    colnames(logistic_res) <- names(coef(mod))
  }
}

plot(dseq, logistic_res[, 'reg_rich', 1], type='b', 
     ylim=range(logistic_res[, 'reg_rich', ]), 
     ylab='Effect of richness on Rib occurrence', 
     xlab='Distance of aggregation')
lines(dseq, logistic_res[, 'reg_rich', 2], lty=2, col='grey')
lines(dseq, logistic_res[, 'reg_rich', 3], lty=2, col='grey')
abline(h=0, lty=3)



library(lme4)
mod <- glmer(rib ~ reg_rich * ISD + 
               (1|sitecode) + (1|speciesname) + (1|Dissection) + (1|year), 
             family='poisson', data=pd)
summary(mod)
confint(mod, method='Wald')

eff <- allEffects(mod)
plot(eff)
#mod <- glmer(rib ~ (1|Dissection) + (1|siteyear) + 
#               (1|speciescode) + ISD + local_rich + reg_rich, family='poisson', 
#             data=pd)
#summary(mod)
