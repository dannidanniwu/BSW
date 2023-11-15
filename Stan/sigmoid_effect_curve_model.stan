data {
  int<lower=0> num_data;            // number of data points
  int num_basis;                    // number of basis
  int<lower=0> num_t_ex;            //number od the maximum exposure time points
  vector[num_data] y;               // response variable
  matrix[num_basis, num_data] B;    // B-spline basis matrix
  int<lower=0, upper=num_t_ex> t_ex[num_data];
  int<lower=0, upper=1> A[num_data];// binary variable
  int<lower=1> num_sites;           // number of sites
  int site[num_data];               // site id for each observation
}

parameters {
  real  beta_0;
  real delta;                     // trt effect
  real<lower=1e-3> sigma;
  real<lower=1e-3> sigma_a;          // variance for site-level spline coefficients
  vector[num_basis] mu_a_raw;         // global mean for spline coefficients
  real beta_A;
  real  k;
  real  eta;
  real<lower=1e-3> lambda;
  real<lower=1e-3> sigma_mu_a; // scale for the a coefficients
  matrix[num_sites, num_basis] a_site_raw;
}

transformed parameters {
  vector[num_data] Y_hat;
  matrix[num_sites, num_basis] a_site;
  vector[num_basis] mu_a;         // global mean for spline coefficients
  

  mu_a[1] = mu_a_raw[1];
  for (i in 2:num_basis){
    mu_a[i] =  mu_a[i-1] + mu_a_raw[i] * sigma_mu_a;
  }
  
  for (j in 1:num_basis) {
    a_site[:,j] =  mu_a[j]+ sigma_a* a_site_raw[:,j];
  }
    
  for (i in 1:num_data) {
    Y_hat[i] = beta_0 + ((delta/(1+k*exp(-1*beta_A*t_ex[i])))- eta)* A[i] + dot_product(a_site[site[i]], B[:,i]);//B[:,i] will give you a vector that consists of all elements from the ith column of that matrix.
  }
}

model {
  beta_0 ~ std_normal();
  delta ~ normal(0, 5);
  k ~ normal(0, 5);
  beta_A ~ normal(0, 5);
  eta ~ normal(0, 5);
  sigma ~ student_t(3, 0, 2.5);
  to_vector(a_site_raw) ~ std_normal();
   
  mu_a_raw ~ std_normal();
  sigma_a ~ student_t(3, 0, 2.5);
  
  lambda ~ student_t(3, 0, 2.5);//large lambda encourage smooth
  sigma_mu_a ~ std_normal();
  
  
  for (i in 2: (num_basis-1)) {
    //target += -0.5 * lambda * square(a[i] - a[i-1]);//First-Order Difference
    target += -0.5 * lambda * square(mu_a[i-1] - mu_a[i] + mu_a[i+1]);
  }

  //print("Before: ", sigma);
  y ~ normal(Y_hat, sigma);
}
generated quantities {
  vector[num_t_ex] phi;
  real phi_mean;
  real phi_max;

  for (i in 1:num_t_ex) {
    phi[i] = (delta/(1+k*exp(-1 * beta_A * t_ex[i]))) - eta;
  }
   phi_mean = mean(phi);
   phi_max = ((delta/(1+k*exp(-1*beta_A*num_t_ex)))- eta);
}