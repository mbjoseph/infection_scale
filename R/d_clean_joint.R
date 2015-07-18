# Cleaning the parasite infection data, and
# defining spatial neighbors

rm(list=ls())
library(readxl)
library(ggplot2)
library(dplyr)

pd <- read.csv("data/coinf.csv")

macroparasites <- c("RIB", "ECHINO", "GLOB", "ALARIA")
microparasites <- c("RVLOAD", "BDZE")

pd1 <- pd[c("AssmtCode_1", "Year", "SiteCode", "HOSTSPECIES", 
            "STAGE", "SVL", "RIB")]
pd1 <- subset(pd1, HOSTSPECIES == "PSRE")
pd1 <- pd1[complete.cases(pd1), ] # remove any NA values
pd1 <- droplevels(pd1)

# extract covariate information on local and regional host richness
host_data <- read.csv("data/ca_site2.csv")
host_data <- subset(host_data, AssmtYear_1 == 2013)
host_data <- droplevels(host_data)
unique(host_data$SiteCode_1)
host_data$site <- as.character(host_data$SiteCode_1)
host_data$site[host_data$site=="WestWing"] <- "WESTWING"
host_data$site[host_data$site=="Barn"] <- "BARN"
host_data$site[host_data$site=="Frog"] <- "FROG"
host_data$site[host_data$site=="HeronGCP"] <- "HERONGCP"
host_data$site[host_data$site=="Hidden"] <- "HIDDEN"
host_data$site[host_data$site=="Ca-SF31"] <- "CA-SF31"
host_data$site[host_data$site=="ca-SF85a"] <- "CA-SF85A"
host_data$site[host_data$site=="Ca-SF101"] <- "CA-SF101"
host_data$site[host_data$site=="Ca-SF25"] <- "CA-SF25"
host_data$site[host_data$site=="Ca-SF79"] <- "CA-SF79"
host_data$site[host_data$site=="EagleGCP"] <- "EAGLEGCP"

# figure out which sites we have parasite and host data for
pd1$site <- as.character(pd1$SiteCode)
setdiff(pd1$site, host_data$site)

pd1 <- subset(pd1, site %in% host_data$site)
pd1 <- droplevels(pd1)
host_data <- subset(host_data, site %in% pd1$site)
host_data <- droplevels(host_data)
host_data$num_site <- as.numeric(host_data$SiteCode_1)

# classify sites into regions based on spatial neighborhoods
coord_d <- read.csv("data/powell.csv")
coord_d$sitename <- as.character(coord_d$ID)

matching_sites <- pd1$SiteCode %in% coord_d$sitename

# print sites that don't match
unique(pd1$SiteCode[!matching_sites])


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

matching_sites <- pd1$SiteCode %in% coord_d$sitename
unique(pd1$SiteCode[!matching_sites])


coord_d$fsite <- as.factor(coord_d$sitename)

coord_d <- subset(coord_d, sitename %in% unique(pd1$SiteCode[matching_sites]))

coord_d <- droplevels(coord_d)
pd1 <- subset(pd1, SiteCode %in% unique(pd1$SiteCode[matching_sites]))
pd1 <- droplevels(pd1)

coord_d$site <- as.numeric(coord_d$fsite)

elev <- coord_d$Elevation[order(coord_d$site[!duplicated(coord_d$site)])]
elev <- scale(elev)

coord_d <- subset(coord_d, !duplicated(coord_d$sitename))

plot(coord_d$Lon, coord_d$Lat)

# define neighborhoods
library(spdep)
s_coord <- coordinates(coord_d[, c("Lon", "Lat")])
# 2 km distance
nbs <- dnearneigh(s_coord, d1=0, d2=5, longlat=T)
plot(nbs, coords=s_coord)
(neighborhoods <- n.comp.nb(nbs))
num_neigh <- neighborhoods$nc
coord_d$nbhood <- neighborhoods$comp.id

all(levels(pd1$SiteCode) == levels(coord_d$fsite))

pd1$numsite <- as.numeric(pd1$SiteCode)

region <- coord_d$nbhood[order(coord_d$site)]

# calculate regional host richness based on regional definition
host_cols <- names(host_data)[c(49:53, 55)]

reg_rich <- rep(NA, num_neigh)
n_in_region <- rep(NA, num_neigh)
for (i in 1:num_neigh){
  sites_in_region <- which(region == i)
  n_in_region[i] <- length(sites_in_region)
  subd <- host_data[host_data$num_site %in% sites_in_region, host_cols]
  reg_rich[i] <- sum(colSums(subd) > 0)
}
plot(n_in_region, reg_rich)
host_data <- host_data[order(host_data$num_site), ]
