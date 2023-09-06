data {
  int<lower=0> N;          // Total number of observations
  int<lower=0> J;          // Number of clusters
  int<lower=0> K_bspline;  // Number of B-spline basis functions
  matrix[N, K_bspline] B;  // B-spline basis matrix for each observation
  vector[N] Y;             // Outcome variable
  int<lower=1, upper=J> cluster[N]; // Cluster assignment for each observation
  vector[N] A;             // Intervention indicator for each observation
}

parameters {
  vector[J] a;                        // Random intercepts for clusters
  matrix[J, K_bspline] beta_spline;   // Coefficients for B-spline basis for each cluster
  real delta;                         // Coefficient for intervention
  real<lower=0> sigma_a;              // Standard deviation for a_j
  real<lower=0> sigma_e;              // Standard deviation for e_{ijk}
  //real beta;                          // overall intercept
}

model {
  // Priors (can be adjusted as needed)
  a ~ normal(0, sigma_a);
  
  for (j in 1:J) {
    beta_spline[j] ~ normal(0, 5);
  }
  
  delta ~ normal(0, 10);         
  sigma_a ~ exponential(1/5); // Exponential with a mean of 2, suggesting most variation is expected within a range of 5.
  sigma_e ~ exponential(1/5); // Same reasoning as above.
  //beta ~ normal(0, 5);
  
  // Likelihood
  for (i in 1:N) {
    Y[i] ~ normal(a[cluster[i]] + dot_product(B[i], beta_spline[cluster[i]]) + delta * A[i], sigma_e);
  }
}
