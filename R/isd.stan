data {
  int<lower=0> n;
  int<lower=0> nsite;
  int<lower=0> nregion;
  int<lower=0> y[n];
  int<lower=1, upper=nsite> site[n];
  int<lower=1> nspec;
  int<lower=1, upper=nspec> h_spec[n];
  int<lower=1, upper=nregion> region[n];
  int<lower=1> nyear;
  int<lower=1, upper=nyear> year[n];
  // covariates
  vector[n] local_richness;
  vector[n] region_richness;
  // snail data
  vector[n] isd; // infected snail density
}

transformed data {
  vector[n] isd2;
  
  isd2 <- isd .* isd;
}

parameters {
  real a0;
  vector[nsite] eta_siteR;
  real<lower=0> sigma_site;
  //vector[nregion] eta_regionR;
  //real<lower=0> sigma_region;
  vector[nspec] eta_speciesR;
  real<lower=0> sigma_species;
  //vector[nyear] eta_yearR;
  //real<lower=0> sigma_year;
  vector[n] eta_indivR;
  real<lower=0> sigma_indiv;
  real beta_local;
  real beta_region;
  real beta_isd;
  real beta_isd2;
}
transformed parameters {
  vector[nsite] eta_site;
  //vector[nregion] eta_region;
  vector[nspec] eta_species;
  //vector[nyear] eta_year;
  vector[n] eta_indiv;
  vector[n] log_mu;
  
  eta_site <- eta_siteR * sigma_site;
  //eta_region <- eta_regionR * sigma_region;
  eta_species <- eta_speciesR * sigma_species;
  //eta_year <- eta_yearR * sigma_year;
  eta_indiv <- eta_indivR * sigma_indiv;

{ // temp scope
  vector[n] site_eff;
  vector[n] spec_eff;
  //vector[n] reg_eff;
  //vector[n] year_eff;
  for (i in 1:n){
    site_eff[i] <- eta_site[site[i]];
    spec_eff[i] <- eta_species[h_spec[i]];
    //reg_eff[i] <- eta_region[region[i]];
    //year_eff[i] <- eta_year[year[i]];
  }
  log_mu <- site_eff //+ reg_eff
              + spec_eff// + year_eff 
              + a0
              + beta_local * local_richness 
              + beta_isd * isd + beta_isd2 + isd2
              + beta_region * region_richness + eta_indiv;
}
}

model {
  a0 ~ normal(0, 1);
  beta_local ~ normal(0, 1);
  beta_region ~ normal(0, 1);
  beta_isd ~ normal(0, 1);
  beta_isd2 ~ normal(0, 1);
  sigma_site ~ normal(0, 1);
  //sigma_region ~ normal(0, 1);
  sigma_species ~ normal(0, 1);
  //sigma_year ~ normal(0, 1);
  eta_siteR ~ normal(0, 1);
  //eta_regionR ~ normal(0, 1);
  eta_speciesR ~ normal(0, 1);
  //eta_yearR ~ normal(0, 1);
  eta_indivR ~ normal(0, 1);
  y ~ poisson_log(log_mu);
}

generated quantities {
  vector[n] log_lik;
  
  for (i in 1:n){
    log_lik[i] <- poisson_log_log(y[i], log_mu[i]);
  }
}