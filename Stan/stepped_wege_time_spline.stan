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
}

parameters {
  real<lower=0> sigma;
  real<lower=0> sigma_a;          // variance for site-level spline coefficients
  matrix[num_sites, num_basis] a_site; // site-specific spline coefficients
  vector[num_basis] mu_a;         // global mean for spline coefficients
  real beta_A;
}

transformed parameters {
  vector[num_data] Y_hat;

  for (i in 1:num_data) {
    Y_hat[i] = beta_A * A[i] + dot_product(a_site[site[i]], B[:,i]);
  }
}

model {
  sigma ~ cauchy(0, 2.5);
  sigma_a ~ cauchy(0, 2.5);
  beta_A ~ normal(0, 10);
  for (j in 1:num_basis) {
    a_site[:,j] ~ normal(mu_a[j], sigma_a);
  }
  mu_a ~ normal(0, 1);
  
  y ~ normal(Y_hat, sigma);
}

generated quantities {
  vector[num_data_test] y_pred_test;

  for (i in 1:num_data_test) {
    y_pred_test[i] = A_test[i] * beta_A + dot_product(mu_a, B_test[:,i]);
  }
}
