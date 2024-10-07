data {
  int<lower=0> num_data;            // number of data points (N)
  int num_basis;                    // number of basis functions (p)
  int<lower=0> num_t_ex;            // number of the maximum exposure time points (T_star)
  vector[num_data] y;               // response variable (Y_ijt)
  matrix[num_basis, num_data] B;    // B-spline basis matrix for study time (B)
  matrix[num_basis, num_data] B_star; // B-spline basis matrix for exposure time (B_star)
  matrix[num_basis, num_t_ex] B_test_star; // B-spline basis matrix for unique exposure time points
  int<lower=0, upper=1> A[num_data]; // binary treatment variable (A_jt)
  int<lower=1> num_clusters;           // number of clusters (J)
  int cluster[num_data];               // cluster id for each observation (s_j)
}

parameters {
  real<lower=1e-3> sigma;           // standard deviation (std) of the observation noise (sigma)
  real<lower=1e-3> sigma_u;         // std of the random effect (sigma_u)
  real<lower=1e-3> sigma_b;         // std for cluster-level study time effect spline coefficients (sigma_b)
  vector[num_basis] beta_raw;       // raw overall study time effect spline coefficients (non-central parameterization)
  vector[num_basis] beta_star_raw;  // raw coefficients for overall time-varying treatment effect (non-central parameterization)
  real alpha;                       // intercept term (alpha)
  real<lower=1e-3> lambda;          // regularization parameter for smoothness (lambda)
  real<lower=1e-3> sigma_beta;      // standard deviation for the beta coefficients (sigma_beta)
  real<lower=1e-3> sigma_beta_star; // standard deviation for the beta star coefficients (sigma_beta_star)
  matrix[num_clusters, num_basis] beta_b_raw; // raw cluster-level study time effect spline coefficients (non-central parameterization)
  vector[num_clusters] u_raw;          // raw random effect (non-central parameterization)
}

transformed parameters {
  vector[num_data] Y_hat;           // expected response variable (E[Y_ijt])
  matrix[num_clusters, num_basis] beta_b; //  cluster-level study time effect spline coefficients (beta_b)
  vector[num_basis] beta;           // overall study time effect spline coefficients (beta)
  vector[num_basis] beta_star;      // overall time-varying treatment effect spline coefficients (beta_star)
  vector[num_clusters] u;              // random effect (u_j)

  // Constructing beta_star using non-central parameterization and random walk prior
  beta_star[1] = beta_star_raw[1];
  for (m in 2:num_basis) {
    beta_star[m] = beta_star[m-1] + beta_star_raw[m] * sigma_beta_star;
  }

  // Constructing beta using non-central parameterization and random walk prior
  beta[1] = beta_raw[1];
  for (m in 2:num_basis) {
    beta[m] = beta[m-1] + beta_raw[m] * sigma_beta;
  }
  
  // Constructing cluster-level coefficients using non-central parameterization and Bayesian hierarchical modeling
  for (m in 1:num_basis) {
    beta_b[:, m] = beta[m] + sigma_b * beta_b_raw[:, m];
  }

  // Applying random effect transformation using non-central parameterization
  u = u_raw * sigma_u;

  // Predicting response
  for (i in 1:num_data) {
    Y_hat[i] = alpha + dot_product(beta_b[cluster[i]], B[:, i]) + exp(u[cluster[i]]) * dot_product(beta_star, B_star[:, i]) * A[i];
  }
}

model {
  alpha ~ std_normal();             // prior for intercept term (alpha)
  sigma ~ student_t(3, 0, 2.5);     // prior for observation noise (sigma)
  sigma_b ~ student_t(3, 0, 2.5);   // prior for cluster-level spline coefficients (sigma_b)
  sigma_u ~ normal(0, 0.2);         // prior for random effect standard deviation (sigma_u)
  sigma_beta ~ std_normal();        // prior for the standard deviation of beta coefficients (sigma_beta)
  sigma_beta_star ~ std_normal();        // prior for the standard deviation of beta coefficients (sigma_beta_star)
  lambda ~ student_t(3, 0, 2.5);    // prior for regularization parameter (lambda)
   
  //parameters for non-central parameterization  
  beta_raw ~ std_normal();          // prior for raw global mean coefficients (beta_star)
  beta_star_raw ~ std_normal();     // prior for raw treatment effect spline coefficients (beta_star)
  to_vector(beta_b_raw) ~ std_normal(); // prior for raw cluster-level coefficients (beta_b)
  u_raw ~ std_normal();             // prior for raw random effect (u_j)

  
  // Smoothness penalty for global mean coefficients
  for (m in 2:(num_basis - 1)) {
    target += -0.5 * lambda * square(beta[m-1] - beta[m] + beta[m+1]);
  }

  // Likelihood function
  y ~ normal(Y_hat, sigma);
}

generated quantities {
  vector[num_t_ex] phi;             // effect estimates at unique exposure times (phi)
  real phi_mean;                    // mean of phi: TATE
  real phi_max;                     // phi at the maximum exposure time point: LTE
  vector[num_data] log_lik;         // log-likelihood for each observation (log_lik)
   
  // Calculating log-likelihood
  for (n in 1:num_data) { 
    log_lik[n] = normal_lpdf(y[n] | Y_hat[n], sigma); 
  } 
  
  // Calculating phi at unique exposure times
  for (i in 1:num_t_ex) {
    phi[i] = dot_product(beta_star, B_test_star[:, i]);
  }
  phi_mean = mean(phi);
  phi_max = phi[num_t_ex];
}
