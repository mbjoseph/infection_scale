# extracting snail infection prevalence and snail density data
rm(list=ls())
library(readxl)
library(ggplot2)
library(dplyr)
library(lme4)
library(tidyr)

# load snail infection data
snail_inf <- read_excel('data/CA_Master_20122014v7_snails.xlsx')
snail_inf <- snail_inf[complete.cases(snail_inf[
  c('HELIRIB_3', "HELI_ SW_1", 'HELDISS3')]), ]

# load snail density data
snail_den <- read.csv('data/snail_d.csv')

library(data.table)
sdt <- data.table(snail_den)
sdt <- sdt[, sum(Count), by=c('SeineCode', 'Sppcode')]

snail_den <- snail_den %>%
  expand(SeineCode, Sppcode) %>%
  data.table() %>%
  merge(sdt, by=c('SeineCode', 'Sppcode'), all.x=TRUE) %>%
  subset(Sppcode == 'HELI') %>%
  as.data.frame()

names(snail_den)[names(snail_den) == 'V1'] <- 'count'
snail_den$count[is.na(snail_den$count)] <- 0

snail_den$site <- snail_den$SeineCode %>%
  as.character() %>%
  strsplit('\\_') %>%
  lapply(FUN = `[`, 1) %>%
  unlist()

# fix some fucked up site names
snail_den$site[snail_den$site == 'Ca-Pig'] <- 'CA-PIG'
snail_den$site[snail_den$site == 'Ca-SF25'] <- 'CA-SF25'
snail_den$site[snail_den$site == 'Ca-SF79'] <- 'CA-SF79'
snail_den$site[snail_den$site == 'Ca-SF101'] <- 'CA-SF101'
snail_den$site[snail_den$site == 'HERONGCP'] <- 'HeronGCP'
snail_den$site[snail_den$site == 'NTALKGCP'] <- 'NTalkGCP'
snail_den$site[snail_den$site == 'RONJRGCP'] <- 'RonJrGCP'

snail_den$date <- snail_den$SeineCode %>%
  as.character() %>%
  strsplit('\\_') %>%
  lapply(FUN = `[`, 2) %>%
  unlist()

snail_den$year <- snail_den$date %>%
  substr(1, 4)

snail_den$assmtcode <- paste(snail_den$site, snail_den$date, sep="_")
snail_den$siteyear <- paste(snail_den$site, snail_den$year, sep="_")

# fill in missing seine records
snail_den$seine <- snail_den$SeineCode%>%
  as.character() %>%
  strsplit('\\_') %>%
  lapply(FUN = `[`, 3) %>%
  unlist()

substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

snail_den$num <- snail_den$SeineCode%>%
  as.character() %>% 
  substrRight(2) %>%
  as.numeric()

# for each assessment code, fill in any gaps in the sweep records
assmtcodes <- unique(snail_den$assmtcode)

for (i in 1:length(assmtcodes)){
  code <- assmtcodes[i]
  indx <- which(snail_den$assmtcode == code)
  missing_sweeps <- !(1:10 %in% snail_den$num[indx])
  if (any(missing_sweeps)){
    n <- sum(missing_sweeps)
    num <- which(missing_sweeps)
    seine_ids <- ifelse(num < 10, paste0('NS01-0', num), 
                        paste0('NS01-', num))
    siteyears <- rep(snail_den$siteyear[indx[1]], n)
    assmt <- rep(code, n)
    year <- rep(snail_den$year[indx[1]], n)
    date <- rep(snail_den$date[indx[1]], n)
    site <- rep(snail_den$site[indx[1]], n)
    count <- rep(0, n)
    spp <- rep('HELI', n)
    SeineCodes <- paste(site, date, seine_ids, sep="_")
    snail_den <- rbind(snail_den, 
                       data.frame(SeineCode = SeineCodes, Sppcode = spp, 
                                  count=count, site=site, date=date, 
                                  year=year, assmtcode=assmt, 
                                  siteyear=siteyears, seine=seine_ids, 
                                  num=num))
  }
}

snail_den <- snail_den[order(snail_den$siteyear, snail_den$date, snail_den$num), ]

# determine overlap between siteyears in the snail density and infection data
setdiff(snail_den$siteyear, snail_inf$siteyear)
keep <- intersect(snail_den$siteyear, snail_inf$siteyear)

# subset data to include only site years common to both data sets
snail_den <- subset(snail_den, siteyear %in% keep)
snail_inf <- subset(snail_inf, siteyear %in% keep)

