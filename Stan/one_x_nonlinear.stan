data {
  int<lower=0> num_data; // number of data points
  int num_basis; //number of basis
  vector[num_data] y;   // response variable
  matrix[num_basis, num_data] B; // B-spline basis matrix
  int<lower=0, upper=1> A[num_data]; // binary variable
  //real<lower=0> lambda; // Regularization parameter
}

parameters {
 // real a0; // Intercept (if necessary)
  row_vector[num_basis] a_raw;
  real<lower=0> sigma;
  real<lower=0> tau;
  real beta_A;
}

transformed parameters {
  row_vector[num_basis] a; 
  vector[num_data] Y_hat; 
  a = a_raw*tau;  
    
  Y_hat = to_vector(A) * beta_A + to_vector(a * B);
  //Y_hat = a0 + to_vector(A) * 5 + to_vector(a * B);
}

model {
  // Priors (modify as required)
  //a0 ~ normal(0, 10);
  a_raw ~ normal(0, 1);
  tau ~ cauchy(0, 1);
  sigma ~ cauchy(0, 2);

  // Regularization term for the B-spline coefficients
  // Adding a penalty to the log posterior for large differences between consecutive B-spline coefficients. This penalty encourages adjacent coefficients to be similar to one another, thus promoting smoothness in the estimated spline function.
  //for (i in 2:num_basis) {
   // target += -0.5 * lambda * square(a[i] - a[i-1]);
 // }

  // Likelihood
  y ~ normal(Y_hat, sigma);
}
