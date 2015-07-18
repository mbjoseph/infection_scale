data {
  int<lower=0> n;
  int<lower=0> nsite;
  int<lower=0> nregion;
  int<lower=0> y[n];
  int<lower=1, upper=nsite> site[n];
  int<lower=1, upper=nregion> region[nsite];
  // covariates
  vector<lower=0>[nsite] local_richness;
  vector<lower=0>[nregion] region_richness;
}
parameters {
  real alpha_0;
  vector[n] alpha_indivR;
  vector[nsite] alpha_siteR;
  vector[nregion] alpha_regionR;
  real<lower=0> sigma_region;
  real<lower=0> sigma_site;
  real<lower=0> sigma_indiv;
  real<lower=0> phi;
  real beta_local;
  real beta_region;
}
transformed parameters {
  vector[n] alpha_indiv;
  vector[nsite] alpha_site;
  vector[nregion] alpha_region;
  vector[n] log_mu;
  
  alpha_indiv <- alpha_indivR * sigma_indiv;
  alpha_site <- beta_local * local_richness  + alpha_siteR * sigma_site;
  alpha_region <- beta_region * region_richness + alpha_regionR * sigma_region;
  
  for (i in 1:n){
    log_mu[i] <- alpha_region[region[site[i]]] 
                  + alpha_site[site[i]];
  }
  log_mu <- alpha_0 + log_mu + alpha_indiv;
}
model {
  alpha_0 ~ normal(0, 1);
  beta_local ~ normal(0, 1.5);
  beta_region ~ normal(0, 1.5);
  sigma_region ~ cauchy(0, 2);
  sigma_site ~ cauchy(0, 2);
  sigma_indiv ~ cauchy(0, 2);
  alpha_indivR ~ normal(0, 1);
  alpha_siteR ~ normal(0, 1);
  alpha_regionR  ~ normal(0, 1);
  phi ~ cauchy(0, 5);
  
  y ~ neg_binomial_2_log(log_mu, phi);
}
generated quantities {
  vector[n] log_lik;
  
  for (i in 1:n){
    log_lik[i] <- neg_binomial_2_log_log(y[i], log_mu[i], phi);
  }
}