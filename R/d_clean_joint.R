# Cleaning the parasite infection data, and
# defining spatial neighbors

rm(list=ls())
library(readxl)
library(ggplot2)
library(dplyr)

pd <- read.csv("data/necropsy20122014_all_rib.csv")

#macroparasites <- c("RIB", "ECHINO", "GLOB", "ALARIA")
#microparasites <- c("RVLOAD", "BDZE")

#pd1 <- pd[c("AssmtCode_1", "Year", "SiteCode", "HOSTSPECIES", 
#            "STAGE", "SVL", "RIB")]
#pd1 <- subset(pd1, HOSTSPECIES == "PSRE")
pd <- pd[complete.cases(pd), ] # remove any NA values
pd <- droplevels(pd)

# load snail infection data
snails <- read_excel('data/CA_Master_20122014v7_snails.xlsx')
nrow(snails)
#snails <- complete.cases(snails)

# extract covariate information on local and regional host richness
names(snails)
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

#snails <- subset(snails, !is.na(local_rich))
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

setdiff(pd$sitename, snails$SiteName)
intersect(pd$sitename, snails$SiteName)

pd <- subset(pd, sitename %in% snails$SiteName)
pd <- droplevels(pd)
snails <- subset(snails, SiteName %in% pd$sitename)
snails <- droplevels(snails)



# classify sites into regions based on spatial neighborhoods
coord_d <- read.csv("data/powell.csv")
coord_d$sitename <- as.character(coord_d$ID)

matching_sites <- pd$sitecode %in% coord_d$sitename

# print sites that don't match
unique(pd$sitecode[!matching_sites])

pd$sitecode[pd$sitecode == "Barn"] <- "BARN"
pd$sitecode[pd$sitecode == "EagleGCP"] <- "EAGLEGCP"
pd$sitecode[pd$sitecode == "Frog"] <- "FROG"
pd$sitecode[pd$sitecode == "Hidden"] <- "HIDDEN"
pd$sitecode[pd$sitecode == "Ca-SF79"] <- "CA-SF79"


coord_d$sitename[coord_d$sitename=="WestWing"] <- "WESTWING"
coord_d$sitename[coord_d$sitename=="Barn"] <- "BARN"
coord_d$sitename[coord_d$sitename=="Frog"] <- "FROG"
coord_d$sitename[coord_d$sitename=="HeronGCP"] <- "HERONGCP"
coord_d$sitename[coord_d$sitename=="Hidden"] <- "HIDDEN"
coord_d$sitename[coord_d$sitename=="Ca-SF31"] <- "CA-SF31"
coord_d$sitename[coord_d$sitename=="ca-SF85a"] <- "CA-SF85A"
coord_d$sitename[coord_d$sitename=="Ca-SF101"] <- "CA-SF101"
coord_d$sitename[coord_d$sitename=="Ca-SF25"] <- "CA-SF25"
coord_d$sitename[coord_d$sitename=="Ca-SF79"] <- "CA-SF79"
coord_d$sitename[coord_d$sitename=="EagleGCP"] <- "EAGLEGCP"
coord_d$sitename[coord_d$sitename=="Beaver"] <- "BEAVER"

matching_sites <- pd$sitecode[pd$sitecode %in% coord_d$sitename]
double_matches <- snails$SiteCode_1[snails$SiteCode_1 %in% matching_sites]
length(unique(double_matches))
# unique(pd$sitecode[!matching_sites])
matches <- unique(double_matches)

coord_d$fsite <- as.factor(coord_d$sitename)

coord_d <- subset(coord_d, sitename %in% matches)

coord_d <- droplevels(coord_d)
pd <- subset(pd, site %in% matches)
pd <- droplevels(pd)
snails <- subset(snails, SiteCode_1 %in% matches)
snails <- droplevels(snails)

pd$fsite <- as.factor(pd$site)

cbind(levels(coord_d$fsite), levels(pd$fsite))

coord_d$site <- as.numeric(coord_d$fsite)

elev <- coord_d$Elevation[order(coord_d$site[!duplicated(coord_d$site)])]
elev <- scale(elev)

coord_d <- subset(coord_d, !duplicated(coord_d$sitename))

par(mfrow=c(1, 2))
plot(coord_d$Lon, coord_d$Lat)
#points(coord_d$Lon, coord_d$Lat, col="red", pch=2)

# define neighborhoods
library(spdep)
s_coord <- coordinates(coord_d[, c("Lon", "Lat")])
# d2 is km distance
nbs <- dnearneigh(s_coord, d1=0, d2=7, longlat=T)
plot(nbs, coords=s_coord)
(neighborhoods <- n.comp.nb(nbs))
num_neigh <- neighborhoods$nc
coord_d$nbhood <- neighborhoods$comp.id

all(levels(pd$fsite) == levels(coord_d$fsite))

pd$numsite <- as.numeric(pd$fsite)

region <- coord_d$nbhood[order(coord_d$site)]

# calculate regional host richness based on regional definition
host_cols <- names(snails)[c(57:63)]

nyear <- length(unique(pd$year))
# assume regional richness is constant across years, otherwise run into issues 
# with uneven sampling, and lack of detection of species
reg_rich <- array(NA, dim=c(num_neigh))
n_sampled <- array(dim=c(num_neigh))
n_in_region <- rep(NA, num_neigh)

snails$num_site <- as.numeric(as.factor(snails$SiteName))

for (i in 1:num_neigh){
  sites_in_region <- which(region == i)
  n_in_region[i] <- length(sites_in_region)
  subd <- snails[(snails$num_site %in% sites_in_region), host_cols]
  n_sampled[i] <- nrow(subd)
  reg_rich[i] <- sum(colSums(subd) > 0)
}
plot(jitter(c(n_sampled), .5), c(reg_rich))
#snails <- snails[order(snails$num_site), ]

pd$region <- region[pd$numsite]

# add regional richness data to pd
pd$reg_rich <- NA
pd$ISD <- NA
snails$ISD <- snails$`HELI_ SW_1` * snails$HELIRIB_3
for (i in 1:nrow(pd)){
  pd$reg_rich[i] <- reg_rich[pd$region[i]]
  snail_row <- which(snails$num_site == pd$numsite[i] & 
                       snails$AssmtYear_1 == pd$year[i])
  if (length(snail_row) == 0){
    next
  } else {
    pd$ISD[i] <- snails$ISD[snail_row]
    if (length(snail_row) != 1){
      print(paste("i =", i, ", length of snail_row:", length(snail_row)))
    }
  }
}
plot(jitter(pd$reg_rich), jitter(pd$local_rich))
str(pd)
pd <- pd[complete.cases(pd), ]
pd <- droplevels(pd)
pd$numsite <- as.numeric(pd$fsite)

# calculate number of snails infected
snails$num_infected <- snails$HELDISS3 * snails$HELIRIB_3

# is Rib present in all regions? 
library(plyr)
ddply(pd, c("region", "year"), 
      summarise, 
      rib_pres = any(rib > 0), 
      nsites = length(site))

library(lme4)
mod <- glmer(rib ~ (1|Dissection) + (1|siteyear) + 
               (1|speciescode) + ISD + local_rich + reg_rich, family='poisson', 
             data=pd)
summary(mod)
