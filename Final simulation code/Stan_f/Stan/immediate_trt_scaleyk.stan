data {
  int<lower=0> num_data;            // number of data points (N)
  int num_basis;                    // number of basis functions (p)
  vector[num_data] y;               // response variable (Y_ijt)
  matrix[num_basis, num_data] B;    // B-spline basis matrix for time (B)
  int<lower=0, upper=1> A[num_data];// binary treatment variable (A_jt)
  int<lower=1> num_clusters;           // number of clusters (J)
  int cluster[num_data];               // cluster id for each observation (j)
  real range_y;
}

parameters {
  real alpha;                       // intercept term (alpha)
  real<lower=1e-3> sigma_b;         // standard deviation (std) for cluster-level spline coefficients (sigma_b)
  vector[num_basis] beta_raw;       // raw overall spline coefficients 
  real tau;                         // treatment A effect (tau)
  real<lower=1e-3> sigma_beta;      // std for the beta coefficients (sigma_beta)
  real<lower=1e-3> lambda;          // regularization parameter for smoothness (lambda)
  matrix[num_clusters, num_basis] b_cluster_raw; // raw cluster-level time effect spline coefficients 
  real<lower=1e-3> sigma;          // std of outcome
}

transformed parameters {
  vector[num_data] Y_hat;           // expected response variable (E[Y_ijt])
  matrix[num_clusters, num_basis] b_cluster; // cluster-level time effect spline coefficients (beta_b)
  vector[num_basis] beta;           // overall time effect spline coefficients (beta)

  beta[1] = beta_raw[1];            // first coefficient of beta is directly from beta_raw, which has a standard normal prior distribution
  for (m in 2:num_basis) {
    beta[m] = beta[m-1] + beta_raw[m] * sigma_beta; // constructing beta iteratively using non-central parameterization ad random walk prior
  }
  
  for (m in 1:num_basis) {
    b_cluster[:, m] = beta[m] + sigma_b * b_cluster_raw[:, m]; // constructing cluster-level coefficients using non-central parameterization and Bayesian hierarchical modeling 
  }
  

  for (i in 1:num_data) {
    Y_hat[i] = alpha + dot_product(b_cluster[cluster[i]], B[:, i]) + tau * A[i]; // B[:,i] will give you a vector that consists of all elements from the ith column of that matrix.
  }
}

model {
  alpha ~ std_normal();             // prior for intercept term (alpha)
  tau ~ normal(0, 5);               // prior for tau
  sigma_b ~ student_t(3, 0, 2.5);   // prior for cluster-level variance (sigma_b)
  sigma_beta ~ std_normal();        // prior for scale of beta coefficients (sigma_beta)
  sigma ~ student_t(3, 0, 2.5);     // prior for outcome y (sigma)
  lambda ~ student_t(3, 0, 2.5);    // prior for regularization parameter (lambda); large lambda encourage smoothness
  
  //parameterS for non-central parameterization
  to_vector(b_cluster_raw) ~ std_normal(); 
  beta_raw ~ std_normal();         
  
  for (m in 2:(num_basis-1)) {
    target += -0.5 * lambda * square(beta[m-1] - beta[m] + beta[m+1]); // smoothness penalty
  }
  
  y ~ normal(Y_hat, sigma);    // likelihood function
}

generated quantities {
  real tau_recovery;
  
   tau_recovery = range_y * tau;
   
}