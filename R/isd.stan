data {
  int<lower=0> n;
  int<lower=0> nsite;
  int<lower=0> nregion;
  int<lower=0> y[n];
  int<lower=1, upper=nsite> site[n];
  int<lower=1, upper=nregion> region[nsite];
  int<lower=1> nyear;
  int<lower=1, upper=nyear> year[n];
  // covariates
  vector[n] richness;
  // snail data
  int nsnail;
  int<lower=0> y_snail[nsnail]; // number of infected snails
  int<lower=0> k_snail[nsnail]; // number of snails tested
  int<lower=1, upper=nsite> site_snail[nsnail]; // site belonging of above
}

parameters {
  real a0;
  vector[nsite] eta_siteR;
  real<lower=0> sigma_site;
  vector[n] eta_indivR;
  real<lower=0> sigma_indiv;
  vector[nyear] eta_year;
  real beta_rich;
  real mu_p;
  vector[nsite] p_siteR;
  real<lower=0> sigma_p;
  vector[nregion] p_regionR;
  real<lower=0> sigma_p_region;
}
transformed parameters {
  vector[nsite] eta_site;
  vector[n] eta_indiv;
  vector[n] log_mu;
  vector[nsnail] logit_p;
  vector[nsite] mu_p_site;
  vector[nregion] p_region;
  vector[nsite] p_site;
  
  eta_site <- eta_siteR * sigma_site;
  eta_indiv <- eta_indivR * sigma_indiv;
  p_region <- p_regionR * sigma_p_region;


{ // temp scope
  vector[n] site_eff;
  vector[n] year_eff;
  for (i in 1:n){
    site_eff[i] <- eta_site[site[i]];
    year_eff[i] <- eta_year[year[i]];
  }
  log_mu <- a0 + site_eff + year_eff + eta_indiv;
}

  for (i in 1:nsite){
    mu_p_site[i] <- p_region[region[site_snail[i]]];
  }
    
  p_site <- mu_p_site + p_siteR * sigma_p;
  
  for (i in 1:nsnail){
    logit_p[i] <- mu_p + p_site[site_snail[i]];
  }
    
}

model {
  // amphib Rib parameters
  a0 ~ normal(0, 3);
  beta_rich ~ normal(0, 3);
  sigma_site ~ normal(0, 3);
  eta_siteR ~ normal(0, 1);
  eta_year ~ normal(0, 1);
  eta_indivR ~ normal(0, 1);

  // snail prevalence parameters
  sigma_p ~ normal(0, 3);
  sigma_p_region ~ normal(0, 2);
  p_regionR ~ normal(0, 1);
  p_siteR ~ normal(0, 1);
  
  // likelihood
  y ~ poisson_log(log_mu);
  y_snail ~ binomial_logit(k_snail, logit_p);
}
