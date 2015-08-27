# evaluating the site X year autocorrelation as a function of value at time t
# fitting a simple random effect model to the rib abundance data
source('R/d_clean.R')
source('R/helpers.R')
library(rstan)
library(maps)
library(scales)
library(plyr)
str(pd)
summ <- ddply(pd, ~ site + year, summarize, 
              mean_rib = mean(rib))
plot(x=NULL, y=NULL, 
     xlim=range(log(1 + summ$mean_rib)), ylim=range(log(1 + summ$mean_rib)), 
     xlab='Mean Rib in year t', ylab='Mean Rib in year t + 1')
abline(0, 1, lty=2)
for (i in 1:length(unique(summ$site))){
  sub <- subset(summ, site == unique(summ$site)[i])
  unique_years <- unique(sub$year)
  if (length(unique_years) > 1){
    sorted_t <- sort(unique_years)
    for (j in 2:length(unique_years)){
      points(x=log(1 + sub$mean_rib[sub$year == unique_years[j - 1]]), 
             y=log(1 + sub$mean_rib[sub$year == unique_years[j]]))
    }
  }
}