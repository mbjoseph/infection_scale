# Exploring larval amphib abundance data
rm(list=ls())
library(readxl)
library(ggplot2)
library(dplyr)

d <- read_excel("R/single_year/Re due Maxv2.xlsx")
str(d)
d$Sppcode[which(d$Sppcode == "BUFO")] <- "BUBO"
pull_site <- function(x){
  split <- strsplit(x, split="_")
  split[[1]][1]
}

d$site <- NULL
for (i in 1:nrow(d)){
  d$site[i] <- pull_site(d$`Survery code`[i])
}

plot(log(sort(d$Count) + 1))
names(d)[names(d) == "Depth (ft) "] <- "depth"
names(d)[names(d) == "Distance (ft)"] <- "distance"
names(d)[names(d) == "Max Depth"] <- "maxdepth"

nz <- d$Area > 0 & d$Perimeter > 0
plot(log(d$Area[nz]), log(d$Perimeter[nz]))

# load snail data
snail_d <- read.csv("R/single_year/snail_d.csv")

# load parasite data
pd <- read.csv("R/single_year/coinf.csv")

# unpack dates, sampling occasions, type and year from seine codes
d$date <- rep(NA, nrow(d))
d$samp <- rep(NA, nrow(d))
d$type <- rep(NA, nrow(d))
d$year <- rep(NA, nrow(d))
for (i in 1:nrow(d)){
  split <- unlist(strsplit(d$SeineCode[i], split="_"))
  d$date[i] <- split[2]
  d$samp[i] <- split[3]
  d$type[i] <- unlist(strsplit(split[3], "-"))[1]
  d$year[i] <- substr(d$date[i], 1, 4)
}
d$year <- as.numeric(d$year)

# fix some inconsistent names for seine codes
d$type[d$type == "ns01"] <- "NS01"
d$type[d$type == "NS1"] <- "NS01"
d$type[d$type == "NSO1"] <- "NS01"
d$type[d$type == "ns02"] <- "NS02"

# do the same to get dates and sampling occasions in snail data
snail_d$date <- rep(NA, nrow(snail_d))
for (i in 1:nrow(snail_d)){
  split <- unlist(strsplit(as.character(snail_d$SeineCode[i]), split="_"))
  snail_d$date[i] <- split[2]
  typ <- unlist(strsplit(split[3], "-"))[1]
  snail_d$samp[i] <- split[3]
}

sweep_d <- subset(d, !is.na(SeineCode) &
                  SizeStage != "ADULT" 
                  & !Sppcode %in% c("AMSP", "AMMA", "PSMA") 
                  & year == 2013)

#snail_d <- subset(snail_d, Sppcode == 'HELI')

# fix site name inconsistencies
sweep_d$site[sweep_d$site == 'Ca-mud67'] <- 'CA-MUD67'
sweep_d$site[sweep_d$site == 'Ca-SF79'] <- 'CA-SF79'
sweep_d$site[sweep_d$site == 'CA-SF85a'] <- 'CA-SF85'
sweep_d$site[sweep_d$site == 'Ca-Pig'] <- 'CA-PIG'
sweep_d$site[sweep_d$site == 'NTalkGCP'] <- 'NTALKGCP'
sweep_d$site[sweep_d$site == 'Ronjrgcp'] <- 'RONJRGCP'
sweep_d$site[sweep_d$site == 'RonJrGCP'] <- 'RONJRGCP'
snail_d$Site.Code[snail_d$Site.Code == 'RonJrGCP'] <- 'RONJRGCP'
snail_d$Site.Code[snail_d$Site.Code == 'Ca-Pig'] <- 'CA-PIG'
unique(sweep_d$site[!(sweep_d$site %in% snail_d$Site.Code)])

library(dplyr)
# combine amphib sweep data with snail sweep data
combined_sweep_d <- data.frame(
  SeineCode = c(sweep_d$SeineCode, as.character(snail_d$SeineCode)),
  Sppcode = c(sweep_d$Sppcode, as.character(snail_d$Sppcode)), 
  Count = c(sweep_d$Count, snail_d$Count), 
  SizeStage = c(sweep_d$SizeStage, as.character(snail_d$SizeStage)), 
  Site = c(sweep_d$site, as.character(snail_d$Site.Code)), 
  date = c(sweep_d$date, snail_d$date)
)

combined_sweep_d$date <- as.character(combined_sweep_d$date)

# set any NAs in counts to zero
library(tidyr)
sd <- spread(combined_sweep_d, Sppcode, Count)
for (i in unique(combined_sweep_d$Sppcode)){
  subd <- unlist(sd[i])
  subd[is.na(subd)] <- 0
  sd[i] <- subd
}

gd <- gather(sd, Species, Count, 
             which(names(sd) %in% levels(combined_sweep_d$Sppcode)))
gd$site <- as.character(gd$Site)

# keep only sites with host AND parasite data
# fix some names
gd$site[gd$site=="WestWing"] <- "WESTWING"
gd$site[gd$site=="Barn"] <- "BARN"
gd$site[gd$site=="Frog"] <- "FROG"
gd$site[gd$site=="HeronGCP"] <- "HERONGCP"
gd$site[gd$site=="Hidden"] <- "HIDDEN"


macroparasites <- c("RIB", "ECHINO", "GLOB", "ALARIA")
microparasites <- c("RVLOAD", "BDZE")

pd1 <- pd[c("AssmtCode_1", "Year", "SiteCode", "HOSTSPECIES", 
            "STAGE", "SVL",
            macroparasites, microparasites)]
pd1 <- pd1[complete.cases(pd1), ] # remove any NA values
pd1 <- subset(pd1, HOSTSPECIES == 'PSRE')
pd1 <- droplevels(pd1)

#na_areas <- which(is.na(gd$Area))
#gd <- gd[-na_areas, ]
gd <- droplevels(gd)
gd$site[gd$site == 'Ca-SF25'] <- 'CA-SF25'
gd$site[gd$site == 'CA-SF85'] <- 'CA-SF85A'
gd$site[gd$site == 'Ca-SF101'] <- 'CA-SF101'
gd$site[gd$site == 'NTALKGCP'] <- 'NTalkGCP'
gd$site[gd$site == 'RONJRGCP'] <- 'RonJrGCP'

keeper_sites <- levels(pd1$SiteCode)[levels(pd1$SiteCode) %in% unique(gd$site)]

sp_rm <- c('GYRA', 'LYCO', 'LYMN', 'PHSP', 'RARI')
gd <- subset(gd, !(Species %in% sp_rm))

gd <- subset(gd, site %in% keeper_sites)
pd1 <- subset(pd1, SiteCode %in% keeper_sites)
gd <- droplevels(gd)
pd1 <- droplevels(pd1)

ggplot(gd, aes(x=log(1 +Count))) + 
  geom_histogram() + 
  facet_wrap(~Species, scales="free")

par(mfrow=c(1, 1))
nsweeps <- c(c(table(gd$site)) / length(unique(gd$Species)))
hist(nsweeps,  breaks=0:max(nsweeps))


# extract snail infection data ---------------------------------------
snail_infection <- read_excel('R/single_year/Snail_Dissection_2013_updated5.xlsx', 
                              sheet=1)
snail_infection$Spp_Code[snail_infection$Spp_Code == 'Heli'] <- "HELI"
snail_infection <- subset(snail_infection, Spp_Code == 'HELI')
str(snail_infection)
# note: GLOB is stylet (styprez)
length(unique(snail_infection$SiteCode))
#sub_snails <- subset(snail_infection, SiteCode %in% keeper_sites)
unique(snail_infection$SiteCode[snail_infection$SiteCode %in% gd$Site])
unique(gd$Site[!gd$site %in% snail_infection$SiteCode])

# fix sitenames 
snail_infection$site <- NA
for (i in 1:length(unique(snail_infection$SiteCode))){
  temp_site <- unique(snail_infection$SiteCode)[i]
  indx <- which(snail_infection$SiteCode == temp_site)
  match <- temp_site %in% levels(factor(gd$site))
  if (match){
    snail_infection$site[indx] <- temp_site
    next
  }
  # try removing CA prefix
  trimmed_site <- strsplit(temp_site, split = '\\-')[[1]][2]
  match_site <- grep(trimmed_site, levels(factor(gd$site)),
                     ignore.case=TRUE, value=TRUE)
  if (deparse(match_site) == 'character(0)'){
    print(paste('no match:', trimmed_site))
    next
  } else {
    snail_infection$site[indx] <- match_site
  }
}

sum(unique(snail_infection$site) %in% levels(factor(gd$site)))
levels(gd$Site)[levels(gd$Site) %in% snail_infection$site]

keep <- levels(gd$Site)[levels(gd$Site) %in% snail_infection$site]

# get spatial coords 
coord_d <- read.csv("R/single_year/powell.csv")
coord_d$sitename <- as.character(coord_d$ID)
coord_d$sitename[coord_d$sitename=="WestWing"] <- "WESTWING"
coord_d$sitename[coord_d$sitename=="Barn"] <- "BARN"
coord_d$sitename[coord_d$sitename=="Frog"] <- "FROG"
coord_d$sitename[coord_d$sitename=="HeronGCP"] <- "HERONGCP"
coord_d$sitename[coord_d$sitename=="Hidden"] <- "HIDDEN"
coord_d$sitename[coord_d$sitename=="Ca-SF31"] <- "CA-SF31"
coord_d$sitename[coord_d$sitename=="ca-SF85a"] <- "CA-SF85A"
coord_d$sitename[coord_d$sitename=="Ca-SF25"] <- "CA-SF25"
coord_d$sitename[coord_d$sitename=="Ca-SF79"] <- "CA-SF79"
coord_d$fsite <- as.factor(coord_d$sitename)

keep <- keep[keep %in% coord_d$sitename]

coord_d <- subset(coord_d, sitename %in% keep, drop=TRUE)
coord_d <- droplevels(coord_d)
coord_d <- coord_d[!duplicated(coord_d$sitename), ]

gd <- droplevels(subset(gd, site %in% keep))
snail_infection <- droplevels(subset(snail_infection, site %in% keep))
snail_infection$fsite <- factor(snail_infection$site)
pd1 <- droplevels(subset(pd1, SiteCode %in% keep))

library(reshape2)
msnail <- melt(snail_infection, 
               id.vars = c('fsite', 'DissectCode'), 
               measure.vars = c('ribprez'))

snail_agg <- ddply(msnail, ~ fsite, 
                   summarize,
                   n_infected = sum(value), 
                   n_sampled = length(value))

# calculate observed local richness
site_sum <- ddply(subset(gd, Species != 'HELI'), 
                  ~ site, summarize, 
                  local_rich = length(unique(Species[Count > 0])))

all.equal(levels(factor(gd$site)), levels(coord_d$fsite))
all.equal(levels(factor(pd1$SiteCode)), levels(coord_d$fsite))
all.equal(levels(factor(snail_infection$site)), levels(coord_d$fsite))
all.equal(levels(factor(pd1$SiteCode)), levels(coord_d$fsite))

coord_d$site <- as.numeric(coord_d$fsite)

plot(coord_d$Lon, coord_d$Lat)

snail_den <- ddply(subset(gd, Species == 'HELI'), 
                   ~ Site + SeineCode, summarize,
                   count = sum(Count))
