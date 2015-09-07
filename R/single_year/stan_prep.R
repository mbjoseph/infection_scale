# preparing data for stan
# Fitting aspatial joint host-symbiont model
source("data/d_clean_joint.R")

library(rstan)

seen <- rep(NA, N_h)
#seen_net <- rep(NA, N_h)

for (i in 1:H){
  for (j in 1:nsite){
    indices <- which(spec_h == i & site_h == j)
    dates <- unique(julian[indices])
    for (k in 1:length(dates)){
      ind2 <- which(spec_h == i & site_h == j & julian == dates[k])
      ever_seen <- any(y_h[ind2] > 0)
      seen[ind2] <- as.numeric(ever_seen)
      #ind3 <- which(spec_h == i & site_h == j & julian == dates[k] & seine==0)
      #seen_net[ind2] <- as.numeric(any(y_h[ind3] > 0))
    }
  }
}

surv_d <- data.frame(spec_h, site_h, y_h, julian, seen)
surv_d <- surv_d[with(surv_d, order(spec_h, site_h, julian, -y_h)), ]

seen_macro <- rep(NA, N_macro)
for (i in 1:Mac){
  for (j in 1:nsite){
    indices <- which(S_macro == i & site_macro == j)
    dates <- unique(macro_d$AssmtCode_1[indices])
    for (k in 1:length(dates)){
      ind2 <- which(S_macro == i & site_macro == j & macro_d$AssmtCode_1 == dates[k])
      ever_seen <- any(y_macro[ind2] > 0)
      seen_macro[ind2] <- as.numeric(ever_seen)
    }
  }
}

seen_micro <- rep(NA, N_micro)
for (i in 1:Mic){
  for (j in 1:nsite){
    indices <- which(S_micro == i & site_micro == j)
    dates <- unique(micro_d$AssmtCode_1[indices])
    for (k in 1:length(dates)){
      ind2 <- which(S_micro == i & site_micro == j & micro_d$AssmtCode_1 == dates[k])
      ever_seen <- any(y_micro[ind2] > 0)
      seen_micro[ind2] <- as.numeric(ever_seen)
    }
  }
}


stan_d <- list(nsite=nsite,
               H=H,
               n_amphib = 5, 
               n_snail = 4,
               S=S,
               N_h=N_h,
               site_h=site_h,
               spec_h=spec_h,
               y_h=y_h,
               net_id = as.numeric(factor(sumd$SeineCode)),
               n_net = length(unique(sumd$SeineCode)),
               seen=seen,
               Mac=Mac,
               N_macro=N_macro,
               site_macro=site_macro,
               H_macro=H_macro,
               S_macro=S_macro,
               hsamp=hsamp,
               tad = tad,
               svl = c(svl),
               macro_host_id=macro_host_id,
               y_macro=y_macro,
               seen_macro = seen_macro,
               Mic=Mic,
               N_micro=N_micro,
               site_micro=site_micro,
               H_micro=H_micro,
               S_micro=S_micro,
               micro_host_id=micro_host_id,
               y_micro=y_micro,
               seen_micro = seen_micro,
               julian_h=c(julian_h),
               julian_macro = c(julian_macro),
               julian_micro = c(julian_micro),
               elev = c(elev),
               area= c(area),
               snails_dissected = nrow(snail_agg), 
               y_snail = snail_agg$n_infected, 
               k_snail = snail_agg$n_sampled)
