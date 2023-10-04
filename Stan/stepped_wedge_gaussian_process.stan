//https://mc-stan.org/docs/stan-users-guide/fit-gp.html
data {
  int<lower=0> num_data;            // number of data points
  vector[num_data] y;               // response variable
  int<lower=0, upper=1> A[num_data];// binary variable
  int<lower=1> num_sites;           // number of sites
  int site[num_data];               // site id for each observation
  int<lower=0> num_times;          // number of unique time points
  real unique_k[num_times];        // unique time values for the GP
  int time_idx[num_data];          // index mapping each data point to its time point in unique_k
 
}

parameters {
  real<lower=0> sigma;             // noise standard deviation
  real beta_A;                 // coefficient for A
  real beta_0;                     // intercept
  vector[num_sites] beta_A_site_raw;
  real<lower=0> sigma_lk;
  vector[num_times] z_site[num_sites];  // New parameter for non-centered parametrization

  real<lower=0> tau;
  real<lower=0> alpha;            // GP length scale for f_overall
  real<lower=0> rho;              // GP scale for f_overall
  
  real<lower=0> alpha_site;       // GP length scale for f_site
  real<lower=0> rho_site;         // GP scale for f_site
  
  vector[num_times] z_overall;  // New parameter for non-centered parametrization of f_overall

}

transformed parameters {
  vector[num_data] Y_hat;          // predicted values
  vector[num_times] f_site[num_sites]; // GP function values for each site
  matrix[num_times, num_times] L_K_site;
  matrix[num_times, num_times] K_site;
  matrix[num_times, num_times] L_K_overall;
  real sq_sigma_lk = square(sigma_lk);
  vector[num_sites] beta_A_site = beta_A + tau * beta_A_site_raw; // Calculate the actual parameter using the raw parameter
  vector[num_times] f_overall;
  matrix[num_times, num_times] K_overall;

  // GP covariance matrix for overall effect
  K_overall = cov_exp_quad(unique_k, alpha, rho);
  for (n in 1:num_times) {
    K_overall[n, n] = K_overall[n, n] + sq_sigma_lk;  // adding jitter for numerical stability
  }
  
  L_K_overall = cholesky_decompose(K_overall);
  
  f_overall = L_K_overall * z_overall;
  K_site = gp_exp_quad_cov(unique_k, alpha_site, rho_site);
  
  for (n in 1:num_times) {
    K_site[n, n] = K_site[n, n] + sq_sigma_lk;  // adding jitter for numerical stability
  }
  
  L_K_site = cholesky_decompose(K_site);
  for (s in 1:num_sites) {
    f_site[s] = f_overall + L_K_site * z_site[s];
  }
  
  

 for (i in 1:num_data) {
    Y_hat[i] = beta_0 + beta_A_site[site[i]] * A[i] + f_site[site[i]][time_idx[i]];
  }
}

model {

  z_overall ~ std_normal();


  for (s in 1:num_sites) {
    z_site[s] ~ std_normal();  // New prior for non-centered parametrization
  }


  tau ~ std_normal();
  
  sigma ~ student_t(3, 0, 2.5);
  beta_A ~ normal(0, 5);
  beta_0 ~ std_normal();
  beta_A_site_raw ~ std_normal();
  
  alpha ~ std_normal();
  rho ~ inv_gamma(5, 5);
  
  alpha_site ~ std_normal();
  rho_site ~ inv_gamma(5, 5);
  
  sigma_lk ~ std_normal();

  y ~ normal(Y_hat, sigma);
}
