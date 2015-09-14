data {
  int n_sweep;
  int n_site;
  int<lower=1, upper=n_site> site_sw[n_sweep];
  int<lower=0> snail_sw[n_sweep];
  int n_dissect;
  int site_dissect[n_dissect];
  int n_shed[n_dissect];
  int n_inf[n_dissect];
  int<lower=0, upper=1> seen[n_site];
}

parameters {
  real b0_den;
  real b0_inf;
  real<lower=0> phi;
  real<lower=0> sd_den;
  real<lower=0> sd_inf;
  vector[n_site] alpha_swR;
  vector[n_site] alpha_infR;
}

transformed parameters {
  vector[n_site] alpha_sw;
  vector[n_site] alpha_inf;
  vector[n_sweep] mu_sw;
  vector[n_dissect] mu_diss;
  
  alpha_sw <- alpha_swR * sd_den + b0_den;
  alpha_inf <- alpha_infR * sd_inf + b0_inf;
  
  for (i in 1:n_sweep){
    mu_sw[i] <- alpha_sw[site_sw[i]];
  }

  for (i in 1:n_dissect){
    mu_diss[i] <- alpha_inf[site_dissect[i]];
  }
}

model {
  b0_den ~ normal(-2, 1.5);
  b0_inf ~ normal(0, 1.5);
  phi ~ normal(0, 1);
  alpha_swR ~ normal(0, 1);
  alpha_infR ~ normal(0, 1);

  snail_sw ~ neg_binomial_2_log(mu_sw, phi);
  n_inf ~ binomial_logit(n_shed, mu_diss);
}