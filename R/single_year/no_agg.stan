data {
  // amphib rib infection
  int<lower=1> n_a; // # amphib infections
  int<lower=1> nsite; // # sites
  int<lower=0> y_a[n_a]; // rib intensities
  int<lower=1, upper=nsite> site_a[n_a]; // site indices 
  vector[nsite] rich;
  // snail infection
  int<lower=1> n_s; // # snail infections
  int<lower=1, upper=nsite> site_s[n_s];
  int<lower=0> inf_s[n_s]; // snail infections
  int<lower=1> samp_s[n_s]; // snails sampled
  // snail density
  int<lower=1> n_sd; // # snail density observations
  int<lower=0> y_sd[n_sd];
  int<lower=1, upper=nsite> site_sd[n_sd];
}

parameters {
  // amphibian infections
  vector[nsite] alpha_aR;
  real mu_a;
  real<lower=0> sd_a;
  real beta_isd;
  real beta_rich;
  real beta_intxn;
  real<lower=0> phi_a;
  
  // snail infections
  vector[nsite] alpha_sR; // site level ranefs
  real<lower=0> sd_s;
  real mu_s;
  
  // snail density
  vector[nsite] alpha_sdR;
  real mu_sd;
  real<lower=0> sd_sd;
  real<lower=0> phi_sd;
}

transformed parameters {
  vector[nsite] alpha_a;
  vector[n_a] alpha_a_vec;
  vector[nsite] alpha_s;
  vector[n_s] alpha_s_vec;
  vector[nsite] alpha_sd;
  vector[n_sd] alpha_sd_vec;
  vector[nsite] isd;
  
  
  alpha_s <- mu_s + alpha_sR * sd_s;
  for (i in 1:n_s){
    alpha_s_vec[i] <- alpha_s[site_s[i]];
  }
  
  alpha_sd <- mu_sd + alpha_sdR * sd_sd;
  for (i in 1:n_sd){
    alpha_sd_vec[i] <- alpha_sd[site_sd[i]];
  }
  
  for (i in 1:nsite){
      isd[i] <- inv_logit(alpha_s[i]);
  }
  
  isd <- isd .* exp(alpha_sd);
  
  alpha_a <- mu_a + alpha_aR * sd_a
              + beta_isd * isd + beta_rich * rich + beta_intxn * (isd .* rich);
  for (i in 1:n_a){
    alpha_a_vec[i] <- alpha_a[site_a[i]];
  }
}

model {
  // priors
  mu_a ~ normal(0, 1);
  sd_a ~ normal(0, 1);
  alpha_aR ~ normal(0, 1);
  beta_isd ~ normal(0, 3);
  beta_rich ~ normal(0, 3);
  beta_intxn ~ normal(0, 3);
  
  mu_s ~ normal(0, 1);
  sd_s ~ normal(0, 1);
  alpha_sR ~ normal(0, 1);
  
  mu_sd ~ normal(0, 1);
  alpha_sdR ~ normal(0, 1);
  sd_sd ~ normal(0, 1);
  
  phi_a ~ normal(0, 5);
  phi_sd ~ normal(0, 5);
  
  // likelihood
  y_a ~ neg_binomial_2_log(alpha_a_vec, phi_a);
  inf_s ~ binomial_logit(samp_s, alpha_s_vec);
  y_sd ~ neg_binomial_2_log(alpha_sd_vec, phi_sd);
}
