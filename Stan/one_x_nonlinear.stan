data {
  int<lower=0> N;          // Number of observations
  int<lower=0> K;          // Number of basis functions for x1
  matrix[N, K] X_bspline;  // Design matrix for B-splines
  vector[N] A;             // Binary variable
  vector[N] y;             // Response variable
}

parameters {
  vector[K] beta_bspline;  // Coefficients for B-splines
  real beta_A;             // Coefficient for A
  real<lower=0> sigma;     // Residual standard deviation
}

model {
  y ~ normal(X_bspline * beta_bspline + beta_A * A, sigma);
}
