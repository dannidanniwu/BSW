//Baysian spline model
//for modeling overall time varying treatment effect and cluster specific time trend
//penalization only through random walk prior
data {
  int<lower=0> num_data;            // number of data points
  int num_basis;                    // number of basis
  int<lower=0> num_t_ex;            //number od the maximum exposure time points
  vector[num_data] y;               // response variable
  matrix[num_basis, num_data] B;    // B-spline basis matrix
  matrix[num_basis, num_data] B_t;    // B-spline basis matrix
  matrix[num_basis, num_t_ex] B_test_t;    // B-spline basis matrix fo unique exposure time, start from 1
  int<lower=0, upper=1> A[num_data];// binary variable
  int<lower=1> num_sites;           // number of sites
  int site[num_data];               // site id for each observation
}

parameters {
  real<lower=1e-3> sigma;
  real<lower=1e-3> sigma_beta_A;
  real<lower=1e-3> sigma_a;          // variance for site-level spline coefficients
  vector[num_basis] mu_a_raw;         // global mean for spline coefficients
  vector[num_basis] beta_A_raw;
  real  beta_0;
  real<lower=1e-3> lambda[2];
  real<lower=1e-3> sigma_mu_a; // scale for the a coefficients
  matrix[num_sites, num_basis] a_site_raw;
}

transformed parameters {
  vector[num_data] Y_hat;
  matrix[num_sites, num_basis] a_site;
  vector[num_basis] mu_a;         // global mean for spline coefficients
  vector[num_basis] beta_A;
  
  beta_A[1] = beta_A_raw[1];
  for (i in 2:num_basis){
    beta_A[i] =  beta_A[i-1] + beta_A_raw[i] * sigma_beta_A;
  }

  mu_a[1] = mu_a_raw[1];
  for (i in 2:num_basis){
    mu_a[i] =  mu_a[i-1] + mu_a_raw[i] * sigma_mu_a;
  }
  
  for (j in 1:num_basis) {
    a_site[:,j] =  mu_a[j]+ sigma_a* a_site_raw[:,j];
  }
    
  for (i in 1:num_data) {
    Y_hat[i] = beta_0 + dot_product(beta_A, B_t[:,i]) * A[i] + dot_product(a_site[site[i]], B[:,i]);//B[:,i] will give you a vector that consists of all elements from the ith column of that matrix.
  }
}

model {
  sigma ~ student_t(3, 0, 2.5);
  beta_A_raw ~ std_normal();
  sigma_beta_A ~ std_normal();
  to_vector(a_site_raw) ~ std_normal();
  mu_a_raw ~ std_normal();
  sigma_a ~ student_t(3, 0, 2.5);
  beta_0 ~ std_normal();
  lambda ~ student_t(3, 0, 2.5);//large lambda encourage smooth
  sigma_mu_a ~ std_normal();
  
  y ~ normal(Y_hat, sigma);
}
generated quantities {
  vector[num_t_ex] phi;
  real phi_mean;
  real phi_max;
  vector[num_data] log_lik;  // num_obs should be the total number of observations 
   
  for (n in 1:num_data) { 
    log_lik[n] = normal_lpdf(y[n] | Y_hat[n], sigma); 
  } 
  
  for (i in 1:num_t_ex) {
    phi[i] = dot_product(beta_A, B_test_t[:,i]);
  }
   phi_mean = mean(phi);
   phi_max = phi[num_t_ex];

}