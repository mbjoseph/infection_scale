# Cleaning the parasite infection data, and
source('R/d_clean.R')

# unaggregated effect of local richness on rib pr(presence)
summary_d <- pd %>%
  group_by(siteyear) %>%
  summarize(rib_pres = as.numeric(any(rib > 0)), 
                   local_rich = unique(local_rich), 
                   site = unique(site), 
                   year = unique(year))
plot(jitter(summary_d$local_rich), 
     jitter(summary_d$rib_pres), 
     xlab='Local amphibian richness', 
     ylab='Any PSRE infected')

pres_mod <- glm(rib_pres ~ local_rich, family='binomial', data=summary_d)
summary(pres_mod)

