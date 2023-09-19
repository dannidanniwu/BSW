library(simstudy)
library(data.table)
library(mgcv)
library(lme4)
library(splines)
library(brms)
library(cmdstanr)
library(slurmR)

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

s_model <- function(ds, mod) {
  #set_cmdstan_path(path = "/gpfs/share/apps/cmdstan/2.25.0")
  fitgam <- gam(y ~ A + s(x1,bs="bs",k=9+4) -1, data = ds, method="REML")
  res_fitgam <- c(summary(fitgam)$p.coeff["A"], summary(fitgam)$se["A"])
  range <-   res_fitgam[1] + c(-1,1) * 1.96 *   res_fitgam[2]
  
  X_bspline <- predict(bs(ds$x1, degree=3, knots=quantile(ds$x1, probs=seq(0, 1, length=11)[-c(1, 11)])))
  
  # Incorporating the B-spline basis into the data
  #ds_with_bspline <- cbind(ds, X_bspline)
  
  # Fitting the model using gam
  #fitgam <- gam(y ~ A + X_bspline - 1, data = ds_with_bspline, method="REML")
  
  
  stan_data <- list(N = nrow(ds),
                    K = ncol(X_bspline),
                    X_bspline = X_bspline,
                    A = ds$A,
                    y = ds$y)
  
  # Fit the model
  fit <- mod$sample(data = stan_data)
  
  diagnostics_df <- as_draws_df(fit$sampler_diagnostics())
  div <- sum(diagnostics_df[, 'divergent__'])
  bayes_gam = fit$summary(variables="beta_A",
                             posterior::default_summary_measures()[1:3],
                             quantiles = ~ quantile(., probs = c(0.025, 0.975)),
                             posterior::default_convergence_measures())
  covered_bayes =   (bayes_gam$`2.5%`< 5 & 5 < bayes_gam$`97.5%`)
  model_results <- data.table(lowci=range[1], upci=range[2], covered_freq=(range[1] < 5 & 5 < range[2]),
                              covered_bayes, bayes_gam,div)
  setnames(model_results, c("lowci.freq", "upci.freq", "covered_freq","covered_bayes","variable","est.mean.bayes",
                            "est.med.bayes","est.sd.bayes",
                            "lowci.bayes",
                            "upci.bayes","rhat","ess_bulk","ess_tail",
                            "div"))
  
  return(model_results)
}

s_single_rep <- function(list_of_defs, mod) {
  
  generated_data <- s_generate(list_of_defs)
  model_results <- s_model(generated_data, mod)
  
  return(model_results)
}

s_replicate <- function(iter, mod) {
  list_of_defs = s_define()
  model_results = s_single_rep(list_of_defs, mod)
  return(data.table(iter=iter, model_results))
}


job <- lapply(1:20,function(i) s_replicate(iter=i, mod=mod))
rbindlist(job)
