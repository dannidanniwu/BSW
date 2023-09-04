#Stan code for GAM from brms package
data {
  int<lower=1> N;  // total number of observations
  vector[N] Y;  // response variable
  int<lower=1> K;  // number of population-level effects
  matrix[N, K] X;  // population-level design matrix
  // data for spline s(k, site, bs = "fs", k = 12)
  int nb_1;  // number of bases
  int knots_1[nb_1];  // This is an integer array of size nb_1, where each entry represents the number of knots for each basis. 
  // basis function matrices
  matrix[N, knots_1[1]] Zs_1_1;#This matrix represents the design matrix for the spline basis functions for the first dimension of the smooth. Each row of the matrix corresponds to an observation and each column corresponds to a basis function. The basis functions are evaluated at the data points for the first dimension of the smoother. The number of columns (basis functions) for this matrix is given by knots_1[1].
  matrix[N, knots_1[2]] Zs_1_2;#design matrix for the spline basis functions for the second dimension of the smooth.
  matrix[N, knots_1[3]] Zs_1_3;
  int prior_only;  // should the likelihood be ignored?
}
transformed data {
  int Kc = K - 1;
  matrix[N, Kc] Xc;  // centered version of X without an intercept
  vector[Kc] means_X;  // column means of X before centering
  for (i in 2:K) {
    means_X[i - 1] = mean(X[, i]);
    Xc[, i - 1] = X[, i] - means_X[i - 1];
  }
}
parameters {
  vector[Kc] b;  // population-level effects
  real Intercept;  // temporary intercept for centered predictors
  // parameters for spline s(k, site, bs = "fs", k = 12)
  // standarized spline coefficients
  vector[knots_1[1]] zs_1_1;
  // standarized spline coefficients
  vector[knots_1[2]] zs_1_2;
  // standarized spline coefficients
  vector[knots_1[3]] zs_1_3;
  real<lower=0> sds_1_1;  // standard deviations of spline coefficients
  real<lower=0> sds_1_2;  // standard deviations of spline coefficients
  real<lower=0> sds_1_3;  // standard deviations of spline coefficients
  real<lower=0> sigma;  // dispersion parameter
}
transformed parameters {
  // actual spline coefficients
  vector[knots_1[1]] s_1_1;
  // actual spline coefficients
  vector[knots_1[2]] s_1_2;
  // actual spline coefficients
  vector[knots_1[3]] s_1_3;
  // compute actual spline coefficients
  s_1_1 = sds_1_1 * zs_1_1;
  // compute actual spline coefficients
  s_1_2 = sds_1_2 * zs_1_2;
  // compute actual spline coefficients
  s_1_3 = sds_1_3 * zs_1_3;
}
model {
  // likelihood including constants
  if (!prior_only) {
    // initialize linear predictor term
    vector[N] mu = Intercept + rep_vector(0.0, N) + Zs_1_1 * s_1_1 + Zs_1_2 * s_1_2 + Zs_1_3 * s_1_3;
    target += normal_id_glm_lpdf(Y | Xc, mu, b, sigma);
  }
  // priors including constants
  target += student_t_lpdf(Intercept | 3, 2.3, 6.8);
  target += student_t_lpdf(sds_1_1 | 3, 0, 6.8)
    - 1 * student_t_lccdf(0 | 3, 0, 6.8);
  target += student_t_lpdf(sds_1_2 | 3, 0, 6.8)
    - 1 * student_t_lccdf(0 | 3, 0, 6.8);
  target += student_t_lpdf(sds_1_3 | 3, 0, 6.8)
    - 1 * student_t_lccdf(0 | 3, 0, 6.8);
  target += std_normal_lpdf(zs_1_1);
  target += std_normal_lpdf(zs_1_2);
  target += std_normal_lpdf(zs_1_3);
  target += student_t_lpdf(sigma | 3, 0, 6.8)
    - 1 * student_t_lccdf(0 | 3, 0, 6.8);
}
generated quantities {
  // actual population-level intercept
  real b_Intercept = Intercept - dot_product(means_X, b);
}