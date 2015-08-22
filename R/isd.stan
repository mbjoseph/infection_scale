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
  // snail infection data
  int nsnail;
  int<lower=0> y_snail[nsnail]; // number of infected snails
  int<lower=0> k_snail[nsnail]; // number of snails tested
  int<lower=1, upper=nsite> site_snail[nsnail]; // site belonging of above
  
  // snail density data
  int<lower=0> n_snail_density;
  vector[n_snail_density] snail_density1;
  vector[n_snail_density] snail_density2;
  int<lower=1, upper=nsite> snail_den_site[n_snail_density];
}

parameters {
  real a0;
  vector[nsite] eta_siteR;
  real<lower=0> sigma_site;
  vector[n] eta_indivR;
  real<lower=0> sigma_indiv;
  vector[nyear] eta_year;
  real beta_rich;
  real<lower=0> beta_isd;
  real beta_intxn;
  
  // snail infection params
  real mu_p;
  vector[nsite] p_siteR;
  real<lower=0> sigma_p;
  vector[nregion] p_regionR;
  real<lower=0> sigma_p_region;
  
  // snail density params
  vector<lower=0>[2] sigma_snail_density;
  real<lower=0> beta_y; // loading for second y variable
  vector[nsite] d_siteR;
  real<lower=0> sigma_d;
  vector[2] mu_den;
  vector[nregion] d_regionR;
  real<lower=0> sigma_d_region;
}
transformed parameters {
  // frog infection
  vector[nsite] eta_site;
  vector[n] eta_indiv;
  vector[n] log_mu;
  
  // snail infection
  vector[nsnail] logit_p;
  vector[nsite] mu_p_site;
  vector[nregion] p_region;
  vector[nsite] p_site;
  
  // density
  vector[nsite] d_site;
  vector[n_snail_density] mu1;
  vector[n_snail_density] mu2;
  vector[nsite] mu_d_site;
  vector[nregion] d_region;
  
  // regional estimate of mean infected snail density
  vector[nregion] regional_isd;
  
  eta_site <- eta_siteR * sigma_site;
  eta_indiv <- eta_indivR * sigma_indiv;
  p_region <- p_regionR * sigma_p_region;

  for (i in 1:nsite){
    mu_p_site[i] <- p_region[region[site_snail[i]]];
  }
    
  p_site <- mu_p_site + p_siteR * sigma_p;
  
  for (i in 1:nsnail){
    logit_p[i] <- mu_p + p_site[site_snail[i]];
  }
    
  // snail density
  d_region <- d_regionR * sigma_d_region;
  
  for (i in 1:nsite) mu_d_site[i] <- d_region[region[snail_den_site[i]]];
  
  d_site <- mu_d_site + d_siteR * sigma_d;
  for (i in 1:n_snail_density){
    mu1[i] <- mu_den[1] + d_site[snail_den_site[i]];
    mu2[i] <- mu_den[2] + beta_y * d_site[snail_den_site[i]];
  }
  
  for (i in 1:nregion)
    regional_isd[i] <- exp(d_region[i]) * inv_logit(p_region[i]);
  
  { // temp scope for Rib expectations
  vector[n] site_eff;
  vector[n] year_eff;
  vector[n] reg_isd;
  for (i in 1:n){
    site_eff[i] <- eta_site[site[i]];
    year_eff[i] <- eta_year[year[i]];
    reg_isd[i] <- regional_isd[region[site[i]]];
  }
  log_mu <- a0 + site_eff + year_eff + eta_indiv
              + beta_rich * richness 
              + beta_isd * reg_isd
              + beta_intxn * (richness .* reg_isd);
  }
}

model {
  // amphib Rib parameters
  a0 ~ normal(0, 5);
  beta_rich ~ normal(0, 5);
  beta_isd ~ normal(0, 5);
  beta_intxn ~ normal(0, 5);
  sigma_site ~ normal(0, 3);
  eta_siteR ~ normal(0, 1);
  eta_year ~ normal(0, 1);
  eta_indivR ~ normal(0, 1);

  // snail prevalence parameters
  sigma_p ~ normal(0, 3);
  sigma_p_region ~ normal(0, 2);
  p_regionR ~ normal(0, 1);
  p_siteR ~ normal(0, 1);
  
  // snail density parameters
  d_siteR ~ normal(0, 1);
  d_regionR ~ normal(0, 1);
  sigma_d ~ normal(0, 3);
  sigma_d_region ~ normal(0, 2);
  beta_y ~ normal(1, .2);
  sigma_snail_density ~ normal(0, 1);
  mu_den ~ normal(0, 1);
  
  // likelihood
  y ~ poisson_log(log_mu);
  y_snail ~ binomial_logit(k_snail, logit_p);
  snail_density1 ~ normal(mu1, sigma_snail_density[1]);
  snail_density2 ~ normal(mu2, sigma_snail_density[2]);
}
