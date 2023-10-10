data {
  int<lower=0> num_data;            // number of data points
  int num_basis;                    // number of basis
  vector[num_data] y;               // response variable
  matrix[num_basis, num_data] B;    // B-spline basis matrix
  int<lower=0, upper=1> A[num_data];// binary variable
  int<lower=1> num_sites;           // number of sites
  int site[num_data];               // site id for each observation

  
}

parameters {
  real<lower=0> sigma;
  real<lower=0> sigma_a;          // variance for site-level spline coefficients
  real<lower=0> tau;
  vector[num_basis] mu_a_raw;         // global mean for spline coefficients
  real beta_A;
  real  beta_0;
  real<lower=0> lambda;
  real<lower=0> sigma_mu_a; // scale for the a coefficients
  matrix[num_sites, num_basis] a_site_raw;
  vector[num_sites] beta_A_site_raw;
  
}

transformed parameters {
  vector[num_data] Y_hat;
  matrix[num_sites, num_basis] a_site;
  vector[num_basis] mu_a;         // global mean for spline coefficients
  vector[num_sites] beta_A_site = beta_A + tau * beta_A_site_raw; // Calculate the actual parameter using the raw parameter
    //print(" beta_A_site: ",  beta_A_site);
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
  
  
    
  for (i in 1:num_data) {
    Y_hat[i] = beta_0 + beta_A_site[site[i]] * A[i] + dot_product(a_site[site[i]], B[:,i]);//B[:,i] will give you a vector that consists of all elements from the ith column of that matrix.
  }
}

model {
  sigma ~ student_t(3, 0, 10);
  beta_A ~ normal(0, 10);

  to_vector(a_site_raw) ~ std_normal();
   
  tau ~ std_normal();
  mu_a_raw ~ std_normal();
  sigma_a ~ student_t(3, 0, 10);
  beta_0 ~ normal(0, 10);
  
  //print("beta_A_site: ", beta_A_site);
  //print("beta_A: ", beta_A);
 // print("tau: ", tau);
  
  
  beta_A_site_raw ~ std_normal();
  lambda ~ student_t(3, 0, 2.5);//large lambda encourage smooth
  sigma_mu_a ~ student_t(3, 0, 5);
  
 
  mu_a_raw ~ std_normal();
  
  for (i in 2: (num_basis-1)) {
    //target += -0.5 * lambda * square(a[i] - a[i-1]);//First-Order Difference
    target += -0.5 * lambda * square(mu_a[i-1] - mu_a[i] + mu_a[i+1]);
  }
  y ~ normal(Y_hat, sigma);
}

generated quantities {
  vector[num_data] y_pred_mean;

  for (i in 1:num_data) {
    y_pred_mean[i] = beta_0 + beta_A_site[site[i]] * A[i] + dot_product(a_site[site[i]], B[:,i]);
  }
}