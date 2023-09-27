data {
  int<lower=0> num_data;            // number of data points
  int num_basis;                    // number of basis
  vector[num_data] y;               // response variable
  matrix[num_basis, num_data] B;    // B-spline basis matrix
  int<lower=0, upper=1> A[num_data];// binary variable
  int<lower=1> num_sites;           // number of sites
  int site[num_data];               // site id for each observation

  // For test data
  int<lower=0> num_data_test;
  matrix[num_basis, num_data_test] B_test;
  int<lower=0, upper=1> A_test[num_data_test];
  int site_test[num_data_test];               // site id for each observation 
}

parameters {
  real<lower=0> sigma;
  real<lower=0> sigma_a;          // variance for site-level spline coefficients
  matrix[num_sites, num_basis] a_site; // site-specific spline coefficients
  vector[num_basis] mu_a;         // global mean for spline coefficients
  vector[num_sites] beta_A_site;
  real beta_A;
  real  beta_0;
  real<lower=0> lambda;
}

transformed parameters {
  vector[num_data] Y_hat;

  for (i in 1:num_data) {
    Y_hat[i] = beta_0 + beta_A_site[site[i]] * A[i] + dot_product(a_site[site[i]], B[:,i]);//B[:,i] will give you a vector that consists of all elements from the ith column of that matrix.
  }
}

model {
  sigma ~ cauchy(0, 2.5);
  sigma_a ~ cauchy(0, 2.5);
  beta_0 ~ normal(0, 1);
  beta_A_site ~ normal(beta_A, 2.5);
  beta_A ~ normal(0, 20);
  lambda ~ cauchy(0, 2.5);//large lambda encourage smooth
  for (j in 1:num_basis) {
    a_site[:,j] ~ normal(mu_a[j], sigma_a);
  }
  mu_a ~ normal(0, 2.5);
  
  for (i in 2: (num_basis-1)) {
    //target += -0.5 * lambda * square(a[i] - a[i-1]);//First-Order Difference
    target += -0.5 * lambda * square(mu_a[i-1] - mu_a[i] + mu_a[i+1]);
  }
  y ~ normal(Y_hat, sigma);
}

