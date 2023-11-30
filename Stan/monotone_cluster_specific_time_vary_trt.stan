data {
  int<lower=0> num_t_ex;   // number of steps in the effect curve
  int<lower=0> num_data;  // number of observations
  int num_basis;                    // number of basis
  matrix[num_basis, num_data] B;    // B-spline basis matrix
  int<lower=0, upper=1> A[num_data];// binary variable
  vector[num_data] y;               // response variable
  int<lower=0, upper=num_t_ex> t_ex[num_data];
  int<lower=1> num_sites;           // number of sites
  int site[num_data];               // site id for each observation
  vector[num_t_ex] c; // fixed constants for the Dirichlet prior
}

parameters {
  real delta;                     // trt effect
  vector[num_basis] mu_a_raw;         // global mean for spline coefficients
  real<lower=0.01, upper=100> omega; // omega, constrained between 0.01 and 100
  real<lower=1e-3> sigma_mu_a; // scale for the a coefficients
  simplex[num_t_ex] alpha_raw;  // raw alpha for the standard Dirichlet prior;the elements are constrained to be non-negative and sum to 1
  real  beta_0;
  real<lower=1e-3> sigma_a; 
  matrix[num_sites, num_basis] a_site_raw;
  real<lower=1e-3> sigma;
  real<lower=1e-3> sigma_u;
  vector[num_sites] u_raw;
}

transformed parameters {
  vector[num_data] Y_hat;
  matrix[num_sites, num_basis] a_site;
  vector[num_basis] mu_a;
  vector<lower=0>[num_t_ex] alpha; // the actual alpha parameters for the model
  vector[num_sites] u;
  
  u= u_raw*sigma_u;
  mu_a[1] = mu_a_raw[1];
   //print(" mu_a[1]: ",  mu_a[1]);
  for (i in 2:num_basis){
    mu_a[i] =  mu_a[i-1] + mu_a_raw[i] * sigma_mu_a;
    //print(" mu_a[", i, "]: ",  mu_a[i]);
  }
  
  for (j in 1:num_basis) {
    a_site[:,j] =  mu_a[j]+ sigma_a* a_site_raw[:,j];
    //print(" sigma_a: ",  sigma_a);
    //print("a_site_raw[", j, "]: ",  a_site_raw[:,j]);  
    //print("a_site[", j, "]: ",  a_site[:,j]);  
  }
  
  
  // Scale the standard Dirichlet-distributed vector by omega
  for (j in 1:(num_t_ex)) {
    alpha[j] = alpha_raw[j] * omega * c[j];
  }

  // Ensure alpha sums to 1
  alpha = alpha / sum(alpha);
  
    
  for (i in 1:num_data) {
    real mu_n = 0; // Initialize mu_n for each data point

    // Only add alpha[t] to mu_n if t > 0
    if (t_ex[i] > 0) {
      for (t in 1:t_ex[i]) {
        mu_n += alpha[t]; // Add the step function effect
      }
    }
   
    // Compute the predicted value Y_hat for each data point
    Y_hat[i] = beta_0 + delta * mu_n * A[i] * exp(u[site[i]]) + dot_product(a_site[site[i]], B[:,i]);
  }
}

model {
  // Priors
  delta ~ normal(0, 5);
  omega ~ uniform(0.01, 100); // Uniform prior for omega
  alpha_raw ~ dirichlet(rep_vector(1, num_t_ex)); // Standard Dirichlet prior for raw alpha
  to_vector(a_site_raw) ~ std_normal();
  sigma ~ student_t(3, 0, 2.5); 
  sigma_u ~ normal(0,0.2);
  u_raw ~ std_normal();
  mu_a_raw ~ std_normal();
  sigma_mu_a ~ std_normal();
  beta_0 ~ std_normal();
  sigma_a ~ student_t(3, 0, 2.5);
   y ~ normal(Y_hat, sigma);
}
generated quantities {
  vector[num_t_ex] phi;
  real phi_mean;
  real mu_n_posterior = 0; 
  vector[num_data] log_lik;  // num_obs should be the total number of observations
  for (n in 1:num_data) {
    log_lik[n] = normal_lpdf(y[n] | Y_hat[n], sigma);
  }
  
  for (i in 1:num_t_ex) {
    mu_n_posterior += alpha[i];
    phi[i] = mu_n_posterior;
  }
   phi_mean = mean(delta * phi);
}