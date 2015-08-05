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
snails$SiteName[snails$SiteName == "Barn"] <- "BARN"
pd$site[pd$site == "Barn"] <- "BARN"
snails$SiteName[snails$SiteName == "Yerba Buena"] <- "YBBA"
snails$SiteName[snails$SiteName == "Frog"] <- "FROG"
pd$site[pd$site == "Frog"] <- "FROG"
snails$SiteName[snails$SiteName == "Hidden"] <- "HIDDEN"
pd$site[pd$site == "Hidden"] <- "HIDDEN"
snails$SiteName[snails$SiteName=="WestWing"] <- "WESTWING"
snails$SiteName[snails$SiteName=="Barn"] <- "BARN"
snails$SiteName[snails$SiteName=="Frog"] <- "FROG"
snails$SiteName[snails$SiteName=="HeronGCP"] <- "HERONGCP"
snails$SiteName[snails$SiteName=="Hidden"] <- "HIDDEN"
snails$SiteName[snails$SiteName=="Ca-SF31"] <- "CA-SF31"
snails$SiteName[snails$SiteName=="ca-SF85a"] <- "CA-SF85A"
snails$SiteName[snails$SiteName=="Ca-SF101"] <- "CA-SF101"
snails$SiteName[snails$SiteName=="Ca-SF25"] <- "CA-SF25"
snails$SiteName[snails$SiteName=="Ca-SF79"] <- "CA-SF79"
snails$SiteName[snails$SiteName=="EagleGCP"] <- "EAGLEGCP"
snails$SiteName[snails$SiteName=="Beaver"] <- "BEAVER"

setdiff(pd$site, snails$SiteName)
intersect(pd$site, snails$SiteName)

pd <- subset(pd, site %in% snails$SiteName)
pd <- droplevels(pd)
snails <- subset(snails, SiteName %in% pd$site)
snails <- droplevels(snails)

snails$num_site <- as.numeric(as.factor(snails$SiteName))

# classify sites into regions based on spatial neighborhoods
coord_d <- read.csv("data/powell.csv")
coord_d$sitename <- as.character(coord_d$ID)

matching_sites <- pd$sitecode %in% coord_d$sitename

# print sites that don't match
unique(pd$sitecode[!matching_sites])


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

matching_sites <- pd$site %in% coord_d$sitename
unique(pd$site[!matching_sites])

coord_d$fsite <- as.factor(coord_d$sitename)

coord_d <- subset(coord_d, sitename %in% unique(pd$site[matching_sites]))

coord_d <- droplevels(coord_d)
pd <- subset(pd, site %in% unique(pd$site[matching_sites]))
pd <- droplevels(pd)

pd$fsite <- as.factor(pd$site)

cbind(levels(coord_d$fsite), levels(pd$fsite))

coord_d$site <- as.numeric(coord_d$fsite)

elev <- coord_d$Elevation[order(coord_d$site[!duplicated(coord_d$site)])]
elev <- scale(elev)

coord_d <- subset(coord_d, !duplicated(coord_d$sitename))

plot(coord_d$Lon, coord_d$Lat)
points(coord_d$Lon, coord_d$Lat, col="red", pch=2)

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
reg_rich <- array(NA, dim=c(num_neigh, nyear))
n_sampled <- array(dim=c(num_neigh, nyear))
n_in_region <- rep(NA, num_neigh)
for (i in 1:num_neigh){
  for (j in 1:nyear){
    sites_in_region <- which(region == i)
    n_in_region[i] <- length(sites_in_region)
    subd <- snails[(snails$num_site %in% sites_in_region) & 
                     (snails$AssmtYear_1 == (2011 + j)), host_cols]
    n_sampled[i, j] <- nrow(subd)
    reg_rich[i, j] <- sum(colSums(subd) > 0)
  }
}
plot(jitter(c(n_sampled), .5), c(reg_rich))
#snails <- snails[order(snails$num_site), ]

pd$region <- region[pd$numsite]

# add regional richness data to pd
pd$reg_rich <- NA
pd$ISD <- NA
snails$ISD <- snails$`HELI_ SW_1` * snails$HELIRIB_3
for (i in 1:nrow(pd)){
  pd$reg_rich[i] <- reg_rich[pd$region[i], pd$year[i] - 2011]
  snail_row <- which(snails$num_site == pd$numsite[i] & 
                       snails$AssmtYear_1 == pd$year[i])
  if (length(snail_row) == 0){
    next
  } else {
    pd$ISD[i] <- snails$ISD[snail_row]
  }
}
plot(jitter(pd$reg_rich), jitter(pd$local_rich))
str(pd)
pd <- pd[complete.cases(pd), ]
pd <- droplevels(pd)
pd$numsite <- as.numeric(pd$fsite)

# calculate number of snails infected
snails$num_infected <- snails$HELDISS3 * snails$HELIRIB_3
