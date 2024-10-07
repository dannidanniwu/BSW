#norm y and k
library(simstudy)
library(data.table)
library(mgcv)
library(lme4)
library(splines)
library(cmdstanr)
library(brms)
library(parallel)
library(dplyr)
library(slurmR)
#references:https://www.rdatagen.net/post/2022-12-13-modeling-the-secular-trend-in-a-stepped-wedge-design/
#references:https://www.rdatagen.net/post/2022-11-01-modeling-secular-trend-in-crt-using-gam/
set_cmdstan_path(path = "/gpfs/share/apps/cmdstan/2.25.0") 
# Compile the Stan model
#mod <- cmdstan_model("./stepped_wedge_random_walk_prior_test.stan")
#mod <- cmdstan_model("./immediate_trt_scaleyk.stan")
mod <- cmdstan_model("/gpfs/home/dw2625/r/Review_sim/stan/immediate_trt_scaleyk.stan")

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
  dd[, normy := (y - min(y))/(max(y) - min(y))]
  dd[] #  generated_data is a data.table
  
}

s_model <- function(train_data, coefA, mod) {
  set_cmdstan_path(path = "/gpfs/share/apps/cmdstan/2.25.0")
   
  # Define the knots based on the training data
  knots_train <- quantile(train_data$k, probs=seq(0, 1, length=6)[-c(1, 6)])
  
  B_train <- predict(splines::bs(train_data$normk, degree=3, knots=quantile(train_data$k, probs=seq(0, 1, length=6)[-c(1, 6)])))
  colnames(B_train) <- paste0("Bspline_", 1:ncol(B_train))
  
  
  stan_data <- list(num_data = nrow(train_data),
                    num_basis = ncol(B_train),
                    B = t(B_train),
                    A = train_data$A,
                    y = train_data$normy,
                    num_clusters = length(unique(train_data$site)),
                    cluster = train_data$site,
                    range_y = max(train_data$y)-min(train_data$y) 
                    
  )
  #mod <- cmdstan_model("./stepped_wedge_random_walk.stan")
  # Fit the Bayesian model
  fit <- mod$sample(data = stan_data,
                    refresh = 0, 
                    iter_warmup = 500,
                    iter_sampling = 2500,
                    adapt_delta = 0.98,#https://mc-stan.org/rstanarm/reference/adapt_delta.html
                    parallel_chains  = 4,
                    show_messages = TRUE) 
  
  diagnostics_df <- as_draws_df(fit$sampler_diagnostics())
  div <- sum(diagnostics_df[, 'divergent__'])
  num_max_tree_depth <- sum(diagnostics_df$treedepth__ >= 10)
  bayes_gam = fit$summary(variables="tau_recovery",
                          posterior::default_summary_measures()[1:3],
                          quantiles = ~ quantile(., probs = c(0.025, 0.975)),
                          posterior::default_convergence_measures())
  
  covered_bayes =   (bayes_gam$`2.5%`< coefA & coefA < bayes_gam$`97.5%`)
  
  
  
  model_results <- data.table( bayes_gam,div,   num_max_tree_depth,
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

scenarios = expand.grid(coefA=c(5),ncluster=10)
i=1
nsim=150
coefA = scenarios[i,"coefA"]
ncluster= scenarios[i,"ncluster"]




sjob <- Slurm_lapply(1:150,
               FUN=s_replicate,
               coefA = coefA,
               ncluster = ncluster,
               mod=mod,
               njobs = 90,
               tmp_path = "/gpfs/scratch/dw2625",
               job_name = "sim3",
               sbatch_opt = list(time = "15:00:00",partition = "cpu_medium", `mem-per-cpu` = "8G"),
               export = c("s_define","s_generate","s_model","s_single_rep"),
               plan = "wait",
               overwrite=TRUE)
res <- Slurm_collect(sjob) # data is a list
#res<- site_plasma_all[lapply(site_plasma_all, function(x) length(x))>1] #filter out the error message
res <- rbindlist(res) # converting list to data.table

date_stamp <- gsub("-", "", Sys.Date())
dir.create(file.path("/gpfs/home/dw2625/r/Review_sim/results/", date_stamp), showWarnings = FALSE)
save(res, file = paste0("/gpfs/home/dw2625/r/Review_sim/results/", date_stamp, "/immediate_trt_scaleyk.rda"))



