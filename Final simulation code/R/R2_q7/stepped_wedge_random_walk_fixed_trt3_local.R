library(simstudy)
library(data.table)
library(mgcv)
library(lme4)
library(splines)
library(cmdstanr)
library(brms)
library(parallel)
library(dplyr)
#references:https://www.rdatagen.net/post/2022-12-13-modeling-the-secular-trend-in-a-stepped-wedge-design/
#references:https://www.rdatagen.net/post/2022-11-01-modeling-secular-trend-in-crt-using-gam/
# Compile the Stan model
#mod <- cmdstan_model("./stepped_wedge_random_walk_prior_test.stan")
mod <- cmdstan_model("./immediate_trt.stan")


s_define <- function() {
  #cluster-specific intercept
  def <- defData(varname = "a", formula = 0, variance = 0.25)
  def <- defData(def, varname = "mu_b", formula = 0, dist = "nonrandom")
  def <- defData(def, varname = "s2_b", formula = 0.36, dist = "nonrandom")
  #A: trt for each cluster and time period
  #b: site-specific time period effect
  defOut <- defDataAdd(varname = "y", formula = " a + b - 0.05 * k^2 + ..coefA * A", variance = 1)
  
  return(list(def = def, defOut = defOut)) 
}

s_generate <- function(iter, coefA, ncluster, list_of_defs) {
  
  set.seed(iter)
  list2env(list_of_defs, envir = environment())
  
  #--- add data generation code ---#
  ds <- genData(ncluster, def, id = "site")#site
  ds <- addPeriods(ds, ncluster+3, "site", perName = "k") #create time periods for each site
  ds <- addCorGen(
    dtOld = ds, idvar = "site", 
    rho = 0.95, corstr = "ar1",
    dist = "normal", param1 = "mu_b", param2 = "s2_b", cnames = "b"
  )
  #assign the treatment status based on the stepped-wedge design
  #per cluster trt change per wave
  ds <- trtStepWedge(ds, "site", nWaves = ncluster, lenWaves = 1, startPer = 2, 
                     grpName = "A", perName = "k")
  ds$site <- as.factor(ds$site)
  #30 individuals per site per period and generate each individual-level outcome
  dd <- genCluster(ds, "timeID", numIndsVar = 10, level1ID = "id")
  dd <- addColumns(defOut, dd)
  dd[, normk := (k - min(k))/(max(k) - min(k))]#scale time period into range 0-1
  
  dd[] #  generated_data is a data.table
  
}

s_model <- function(train_data, coefA, mod) {
  # #######Fitting the GAM model with default penalization 
 
  # Define the knots based on the training data
  knots_train <- quantile(train_data$k, probs=seq(0, 1, length=6)[-c(1, 6)])
  
  B_train <- predict(splines::bs(train_data$k, degree=3, knots=quantile(train_data$k, probs=seq(0, 1, length=6)[-c(1, 6)])))
  colnames(B_train) <- paste0("Bspline_", 1:ncol(B_train))
  
  
  stan_data <- list(num_data = nrow(train_data),
                    num_basis = ncol(B_train),
                    B = t(B_train),
                    A = train_data$A,
                    y = train_data$y,
                    num_clusters = length(unique(train_data$site)),
                    cluster = train_data$site
                    
  )
  #mod <- cmdstan_model("./stepped_wedge_random_walk.stan")
  # Fit the Bayesian model
  fit <- mod$sample(data = stan_data,
                    refresh = 0, 
                    iter_warmup = 500,
                    iter_sampling = 2000,
                    parallel_chains  = 4,
                    show_messages = TRUE) 
  
  diagnostics_df <- as_draws_df(fit$sampler_diagnostics())
  div <- sum(diagnostics_df[, 'divergent__'])
  bayes_gam = fit$summary(variables="tau",
                          posterior::default_summary_measures()[1:3],
                          quantiles = ~ quantile(., probs = c(0.025, 0.975)),
                          posterior::default_convergence_measures())
  
  covered_bayes =   (bayes_gam$`2.5%`< coefA & coefA < bayes_gam$`97.5%`)
  
  
  
  model_results <- data.table( bayes_gam,div,
                               covered_bayes) 
  
  return(model_results)
}

s_single_rep <- function(iter, coefA, ncluster, list_of_defs, mod) {
  
  train_data <- s_generate(iter, coefA, ncluster, list_of_defs)
  
  model_results <- s_model(train_data, coefA, mod)
  
  return(model_results)
}

s_replicate <- function(iter, coefA, ncluster, mod) {
  list_of_defs = s_define()
  model_results = s_single_rep(iter, coefA, ncluster,list_of_defs, mod)
  return(data.table(iter=iter, coefA = coefA , ncluster=ncluster, model_results))
}

scenarios = expand.grid(coefA=5,ncluster=10)
# scenarios = expand.grid(coefA=5,ncluster=10)
i=1

coefA = scenarios[i,"coefA"]
ncluster= scenarios[i,"ncluster"]
# res <- replicate(1, s_replicate(iter=1,coefA = coefA,ncluster=ncluster,
#                                 mod=mod
#                                     ))
start_time <- Sys.time()

res <- lapply(
    X = 1, 
    FUN = function(x) s_replicate(iter=x,coefA = coefA,ncluster=ncluster,
                                  mod=mod))

end_time <- Sys.time()

total_time <- end_time - start_time

# Print the total time taken
print(total_time)
