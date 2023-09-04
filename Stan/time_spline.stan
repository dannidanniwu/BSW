// stan code for this model:Mixed effects model with fixed spline and random intercept
//Y_{ijk} = beta + a_{j} + s(k) + \delta A_{jk} + e_{ijk}
//where $s(.)$ is a smooth function of time; $Y_{ijk}$ is the (continuous) outcome of individual $i$ in cluster $j$ during time period $k$. $a_j$ is the random intercept for site $j$, and we assume that $a_j \sim N(0, \sigma^2_a)$. $A_{jk}$ is the intervention indicator for site $j$ during time period $k$. $\beta_k$ is a period-specific effect. And $e_{ijk}$ is the individual level effect, $e_{ijk} \sim N(0, \sigma^2_e)$. Using B-sline for the smooth function
data {
  int<lower=0> N;          // Total number of observations
  int<lower=0> J;          // Number of clusters
  int<lower=0> K_bspline;  // Number of B-spline basis functions
  matrix[N, K_bspline] B;  // B-spline basis matrix
  vector[N] Y;             // Outcome variable
  int<lower=1, upper=J> cluster[N]; // Cluster assignment for each observation
  vector[N] A;             // Intervention indicator for each observation
}

parameters {
  vector[J] a;                     // Random intercepts for clusters
  vector[K_bspline] beta_spline;   // Coefficients for B-spline basis
  real delta;                      // Coefficient for intervention
  real<lower=0> sigma_a;           // Standard deviation for a_j
  real<lower=0> sigma_e;           // Standard deviation for e_{ijk}
  real beta;                       //overall intercept
}

model {
  // Priors (can be adjusted as needed)
  a ~ normal(0, 5);
  beta_spline ~ normal(0, 5);   
  delta ~ normal(0, 10);         
  sigma_a ~ normal(0, 5);
  sigma_e ~ normal(0, 5);
  beta ~ normal(0, 5);
  
  // Likelihood
  Y ~ normal(beta + a[cluster] + B * beta_spline + delta * A, sigma_e);
}
