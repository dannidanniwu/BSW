data {
  int<lower=0> num_data; // number of data points
  int num_basis; //number of basis
  vector[num_data] y;   // response variable
  matrix[num_basis, num_data] B; // B-spline basis matrix
  int<lower=0, upper=1> A[num_data]; // binary variable
  //real<lower=0> lambda; // Regularization parameter
  // For test data
  int<lower=0> num_data_test;
  matrix[num_basis, num_data_test] B_test;
  int<lower=0, upper=1> A_test[num_data_test];
}

parameters {
 // real a0; // Intercept (if necessary)
  row_vector[num_basis] a_raw;
  real<lower=0> sigma;
  real<lower=0> sqrt_tau;
  real<lower=0> lambda;
  real beta_A;
}

transformed parameters {
  real<lower=0> tau = square(sqrt_tau);
  row_vector[num_basis] a; 
  vector[num_data] Y_hat; 
  a[1] = a_raw[1];
  for (i in 2:num_basis)
    a[i] = a[i-1] + a_raw[i]*tau;

  Y_hat = to_vector(A) * beta_A + to_vector(a * B);
  //Y_hat = a0 + to_vector(A) * 5 + to_vector(a * B);
}

model {
  // Priors (modify as required)
  //a0 ~ normal(0, 10);
  a_raw ~ normal(0, 1);
  sqrt_tau ~ normal(0, 1);
  sigma ~ cauchy(0, 2.5);
  lambda ~ cauchy(0, 2.5);//large lambda encourage smooth
  // Regularization term for the B-spline coefficients
  // Adding a penalty to the log posterior for large differences between consecutive B-spline coefficients. This penalty encourages adjacent coefficients to be similar to one another, thus promoting smoothness in the estimated spline function.
  for (i in 2: (num_basis-1)) {
    //target += -0.5 * lambda * square(a[i] - a[i-1]);//First-Order Difference
    target += -0.5 * lambda * square(a[i-1] - a[i] + a[i+1]);
  }

  // Likelihood
  y ~ normal(Y_hat, sigma);
}

generated quantities {
  vector[num_data_test] y_pred_test;

  for (i in 1:num_data_test) {
    y_pred_test[i] = A_test[i] * beta_A + dot_product(to_row_vector(B_test[,i]), a);
  }
}
