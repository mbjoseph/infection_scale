# fitting simple snail density/infection model
source('R/eiv-extravaganza/d_clean.R')
source('R/helpers.R')
stan_d <- list(
  n_sweep = nrow(snail_den), 
  n_site = length(levels(snail_den$fsite)),
  site_sw = as.numeric(snail_den$fsite), 
  snail_sw = snail_den$count, 
  n_dissect = nrow(snail_inf), 
  site_dissect = as.numeric(snail_inf$fsite), 
  n_shed = snail_inf$HELDISS3, 
  n_inf = as.integer(snail_inf$n_infected), 
  seen = seen$seen
)

library(rstan)
watch <- c('b0_den', 'b0_inf', 'phi', 'sd_den', 'sd_inf', 
           'alpha_sw', 'alpha_inf')
m_init <- stan('R/eiv-extravaganza/mod.stan', data=stan_d, chains=1, iter=10)
m_fit <- stan(fit=m_init, data=stan_d, cores=3, pars=watch, chains=3)
traceplot(m_fit, inc_warmup=TRUE)
par(mar=c(5, 4, 4, 2) + 0.1)
library(ggmcmc)
ggd <- ggs(m_fit)
ggs_caterpillar(ggd, 'alpha_sw')
ggs_caterpillar(ggd, 'alpha_inf')

post <- extract(m_fit)
alpha_inf <- apply(post$alpha_inf, 2, 
      function(x){
        hdi <- HDI(x)
        c(hdi[1], median(x), hdi[2])
      })
Palpha_inf <- data.frame(t(plogis(alpha_inf)))
alpha_inf <- data.frame(t(alpha_inf))

names(alpha_inf) <- c('lo', 'med', 'hi')
names(Palpha_inf) <- c('lo', 'med', 'hi')
ord <- order(Palpha_inf$med)
par(mfrow=c(1, 1))
plot(x=Palpha_inf$med[ord], y=1:nrow(Palpha_inf), xlim=c(0, 1))
segments(x0=Palpha_inf$lo[ord], x1=Palpha_inf$hi[ord], 
         y0=1:nrow(Palpha_inf), y1=1:nrow(Palpha_inf))

plot(stan_d$n_shed, y=Palpha_inf$med, ylim=c(0, 1), 
     xlab='Number of snails shed', 
     ylab='Posterior probability of infection')
segments(y0=Palpha_inf$lo, y1=Palpha_inf$hi, 
         x0=stan_d$n_shed, x1=stan_d$n_shed)

alpha_sw <- apply(post$alpha_sw, 2, 
                   function(x){
                     hdi <- HDI(x)
                     c(hdi[1], median(x), hdi[2])
                   })
alpha_sw <- data.frame(t(alpha_sw))
names(alpha_sw) <- c('lo', 'med', 'hi')

plot(alpha_sw$med, alpha_inf$med, ylim=c(-12, 2), xlim=c(-4, 5))
segments(y0=alpha_inf$lo, y1=alpha_inf$hi, 
         x0=alpha_sw$med, x1=alpha_sw$med)
segments(x0=alpha_sw$lo, x1=alpha_sw$hi, 
         y0=alpha_inf$med, y1=alpha_inf$med)
