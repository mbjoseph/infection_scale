data {
  int<lower=0> n;
  int<lower=0> nsite;
  int<lower=0> nregion;
  int<lower=0> y[n];
  int<lower=1, upper=nsite> site[n];
  int<lower=1> nspec;
  int<lower=1, upper=nspec> h_spec[n];
  int<lower=1, upper=nregion> region[nsite];
  // covariates
  vector[nsite] local_richness;
  vector[nregion] region_richness;
  // snail data
  int<lower=0> n_snails;
  int<lower=0, upper=1> y_snails[n_snails];
  int<lower=1, upper=nsite> site_snails[n_snails];
}
transformed data {
  vector[2] zeros;
  
  for (i in 1:2){
    zeros[i] <- 0;
  }
}
parameters {
  corr_matrix[2] Rho;
  real alpha_site;
  real beta_site;
  real beta_rho_site;
  vector[2] mu;
  vector<lower=0>[2] sigma_site;
  matrix[nsite, 2] eta_site;
  vector<lower=0>[2] sigma_region;
  matrix[nregion, 2] eta_region;
  real alpha_region;
  real beta_region;
  real beta_rho_region;
  real<lower=0> sigma_species;
  vector[nspec] eta_species;
}
transformed parameters {
  vector[n_snails] logit_p;
  vector[n] log_mu;
  real sig_prod;
  vector[2] var_vec;
  real sig_reg_prod;
  vector[2] var_reg_vec;

  sig_prod <- sigma_site[1] * sigma_site[2];
  var_vec <- sigma_site .* sigma_site;
  
  sig_reg_prod <- sigma_region[1] * sigma_region[2];
  var_reg_vec <- sigma_region .* sigma_region;
  
  for (i in 1:n){
    log_mu[i] <- eta_site[site[i], 1] + eta_species[h_spec[i]] 
    + eta_region[region[site[i]], 1];
  }

  for (i in 1:n_snails){
    logit_p[i] <- eta_site[site_snails[i], 2] 
    + eta_region[region[site_snails[i]], 2];
  }

}
model {
  alpha_site ~ normal(0, 3);
  beta_site ~ normal(0, 3);
  beta_rho_site ~ normal(0, 3);
  sigma_site ~ normal(0, 3);
  alpha_region ~ normal(0, 3);
  beta_region ~ normal(0, 3);
  beta_rho_region ~ normal(0, 3);
  sigma_region ~ normal(0, 3);
  sigma_species ~ normal(0, 3);
  mu ~ normal(0, 1.5);
  eta_species ~ normal(0, sigma_species);
  { matrix[2, 2] Sigma;
    real x;
    real rho;
    for (i in 1:nsite){
    Sigma <- diag_matrix(var_vec);
    x <- alpha_site + beta_rho_site * local_richness[i];
    rho <- 2 * exp(x) / (1 + exp(x)) - 1;
    Sigma[1, 2] <- rho * sig_prod;
    Sigma[2, 1] <- Sigma[1, 2];
    eta_site[i] ~ multi_normal(mu + beta_site * local_richness[i], Sigma);
  }
  }
  
  // regional ranefs
  { matrix[2, 2] Sigma;
    real x;
    real rho;
    for (i in 1:nregion){
    Sigma <- diag_matrix(var_reg_vec);
    x <- alpha_region + beta_rho_region * region_richness[i];
    rho <- 2 * exp(x) / (1 + exp(x)) - 1;
    Sigma[1, 2] <- rho * sig_reg_prod;
    Sigma[2, 1] <- Sigma[1, 2];
    eta_region[i] ~ multi_normal(zeros + beta_region * region_richness[i], Sigma);
  }
  }
  
  y_snails ~ bernoulli_logit(logit_p);
  y ~ poisson_log(log_mu);
}
