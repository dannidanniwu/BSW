library(simstudy)
library(data.table)
library(mgcv)
library(lme4)
library(splines)
library(brms)
library(cmdstanr)
library(slurmR)
library(dplyr)
# Compile the Stan model
mod <- cmdstan_model("./one_x_nonlinear.stan")

s_define <- function() {
  #x
  def <- defData(varname = "x1", formula = 0, variance = 2)
  def <- defData(def, varname = "A", formula = 0.5, dist = "binary")
  #y
  defOut <- defDataAdd(varname = "y", formula = "exp(-(x1 - 0.5)^2)  + 5 * A", variance = 0.25)
  #exp(-(u - 0.5)^2) 
  return(list(def = def, defOut = defOut)) 
}

s_generate <- function(list_of_defs) {
  
  list2env(list_of_defs, envir = environment())
  
  #--- add data generation code ---#
  #24 sites in total
  ds <- genData(200, def)
  ds <- addColumns(defOut, ds)
  #summary(ds)
  ds[] #  generated_data is a data.table
}

s_model <- function(train_data, test_data, mod) {
  #set_cmdstan_path(path = "/gpfs/share/apps/cmdstan/2.25.0")
  fitgam <- gam(y ~ A + s(x1,bs="bs",k=9+4) -1, data = train_data, method="REML")
  res_fitgam <- c(summary(fitgam)$p.coeff["A"], summary(fitgam)$se["A"])
  range <-   res_fitgam[1] + c(-1,1) * 1.96 *   res_fitgam[2]
  
  # Define the knots based on the training data
  knots_train <- quantile(train_data$x1, probs=seq(0, 1, length=11)[-c(1, 11)])
  
  B_train <- t(predict(bs(train_data$x1, degree=3, knots=quantile(train_data$x1, probs=seq(0, 1, length=11)[-c(1, 11)]))))
  
  # Incorporating the B-spline basis into the data
  #ds_with_bspline <- cbind(train_data, X_bspline)
  
  # Fitting the model using gam
  #fitgam2 <- gam(y ~ A + X_bspline - 1, data = ds_with_bspline, method="REML")
  
  # Compute the B-spline basis for the test set using the same knots from the training data
  B_test <- t(predict(bs(test_data$x1, degree=3, knots=knots_train)))
  
  
  stan_data <- list(num_data = nrow(train_data),
                    num_basis = nrow(B),
                    B = B_train,
                    A = train_data$A,
                    y = train_data$y,
                    num_data_test = nrow(test_data),
                    B_test = B_test,
                    A_test = test_data$A)
  
  # Fit the model
  fit <- mod$sample(data = stan_data,
                    refresh = 0,
                    chains = 4L,
                    parallel_chains = 4L,
                    iter_warmup = 500,
                    iter_sampling = 2500,
                    show_messages = FALSE)
  
  diagnostics_df <- as_draws_df(fit$sampler_diagnostics())
  div <- sum(diagnostics_df[, 'divergent__'])
  bayes_gam = fit$summary(variables="beta_A",
                             posterior::default_summary_measures()[1:3],
                             quantiles = ~ quantile(., probs = c(0.025, 0.975)),
                             posterior::default_convergence_measures())
  covered_bayes =   (bayes_gam$`2.5%`< 5 & 5 < bayes_gam$`97.5%`)
  
  ##############Test RMSE on the test data#################
  # Use the GAM model to predict on the test data to get the mean predictions.
  test_data$pred_gam <- predict(fitgam, newdata = test_data)
  
  
  #posterior_samples <- data.table(as_draws_df(fit$draws()))
  #a <- as.matrix(posterior_samples%>%select(paste0("a","[",1:length(knots_train)+3,"]")))
  #beta_A <- as.matrix(posterior_samples%>%select(paste0("beta_A")))
  
  posterior_predictions <-fit$summary(variables="y_pred_test",posterior::default_summary_measures()[1])
  test_data$pred_bayesian <- posterior_predictions$mean
  # Calculate prediction error for both models
  test_data$error_gam <- (test_data$y - test_data$pred_gam)^2
  test_data$error_bayesian <- (test_data$y - test_data$pred_bayesian)^2
  
  # Calculate RMSE for both models
  rmse_gam <- sqrt(mean(test_data$error_gam))
  rmse_bayesian <- sqrt(mean(test_data$error_bayesian))
  
  model_results <- data.table(lowci=range[1], upci=range[2], bayes_gam,div, covered_freq=(range[1] < 5 & 5 < range[2]),
                              covered_bayes,rmse_gam, rmse_bayesian)
  
  setnames(model_results, c("lowci.freq", "upci.freq", "covered_freq","covered_bayes","variable","est.mean.bayes",
                            "est.med.bayes","est.sd.bayes",
                            "lowci.bayes",
                            "upci.bayes","rhat","ess_bulk","ess_tail",
                            "div","rmse_gam","rmse_bayesian"))
  
  return(model_results)
}

s_single_rep <- function(list_of_defs, mod) {
  
  train_data <- s_generate(list_of_defs)
  test_data <- s_generate(list_of_defs)
  model_results <- s_model(train_data, test_data, mod)

  return(model_results)
}

s_replicate <- function(iter, mod) {
  list_of_defs = s_define()
  model_results = s_single_rep(list_of_defs, mod)
  return(data.table(iter=iter, model_results))
}


job <- lapply(1:20,function(i) s_replicate(iter=i, mod=mod))
rbindlist(job)
