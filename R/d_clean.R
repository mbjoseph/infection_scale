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

# load snail infection data
snails <- read_excel('data/CA_Master_20122014v7_snails.xlsx')
snails <- snails[complete.cases(snails[c('HELIRIB_3', "HELI_ SW_1", 
                                         "HELDISS3")]), ]
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

setdiff(pd$sitename, snails$SiteName)
intersect(pd$sitename, snails$SiteName)

pd <- subset(pd, sitename %in% snails$SiteName)
pd <- droplevels(pd)
snails <- subset(snails, SiteName %in% pd$sitename)
snails <- droplevels(snails)


setdiff(snails$siteyear, pd$siteyear)
intersect(snails$siteyear, pd$siteyear)

snails <- subset(snails, siteyear %in% pd$siteyear)
pd <- subset(pd, siteyear %in% snails$siteyear)
snails <- droplevels(snails)
pd <- droplevels(pd)

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