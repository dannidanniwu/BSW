data {
  int<lower=0> num_data;            // number of data points
  real k[num_data];                 // The time or feature for the GP
  vector[num_data] y;               // response variable
  int<lower=0, upper=1> A[num_data];// binary variable
  int<lower=1> num_sites;           // number of sites
  int site[num_data];               // site id for each observation
}

parameters {
  real<lower=0> sigma;             // noise standard deviation
  real beta_A_raw;                 // coefficient for A
  real beta_0;                     // intercept
  vector[num_sites] beta_A_site_raw;
  real<lower=0> sigma_lk;
  
  real<lower=0> tau;
  real<lower=0> alpha;            // GP length scale for f_overall
  real<lower=0> rho;              // GP scale for f_overall
  
  real<lower=0> alpha_site;       // GP length scale for f_site
  real<lower=0> rho_site;         // GP scale for f_site
  
  vector[num_data] f_overall;     // Overall GP function values
  vector[num_data] f_site[num_sites]; // GP function values for each site
}

transformed parameters {
  vector[num_data] Y_hat;          // predicted values
  
  real beta_A = 0 + 5 * beta_A_raw;
  vector[num_sites] beta_A_site = beta_A + tau * beta_A_site_raw; // Calculate the actual parameter using the raw parameter

  for (i in 1:num_data) {
    Y_hat[i] = beta_0 + beta_A_site[site[i]] * A[i] + f_site[site[i]][i];
  }
}

model {
  matrix[num_data, num_data] L_K_overall;
  matrix[num_data, num_data] K_overall;
  
  matrix[num_data, num_data] L_K_site;
  matrix[num_data, num_data] K_site;
  
  real sq_sigma_lk = square(sigma_lk);

  // GP covariance matrix for overall effect
  K_overall = cov_exp_quad(k, alpha, rho);
  for (n in 1:num_data) {
    K_overall[n, n] = K_overall[n, n] + sq_sigma_lk;  // adding jitter for numerical stability
  }
  
  L_K_overall = cholesky_decompose(K_overall);

  // GP covariance matrix for site-specific effects
  K_site = cov_exp_quad(k, alpha_site, rho_site);
  for (n in 1:num_data) {
    K_site[n, n] = K_site[n, n] + sq_sigma_lk;  // adding jitter for numerical stability
  }

  L_K_site = cholesky_decompose(K_site);

  tau ~ normal(0, 1);
  
  sigma ~ student_t(3, 0, 2.5);
  beta_A_raw ~ normal(0, 1);
  beta_0 ~ normal(0, 1);
  beta_A_site_raw ~ normal(0, 1);
  
  alpha ~ normal(0, 1);
  rho ~ inv_gamma(5, 5);
  
  alpha_site ~ normal(0, 1);
  rho_site ~ inv_gamma(5, 5);
  
  f_overall ~ multi_normal_cholesky(rep_vector(0, num_data), L_K_overall);
  sigma_lk ~ std_normal();

  // Hierarchical GPs for sites
  for (s in 1:num_sites) {
    f_site[s] ~ multi_normal_cholesky(f_overall, L_K_site);
  }

  y ~ normal(Y_hat, sigma);
}
