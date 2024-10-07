library(simstudy)
library(data.table)
library(mgcv)
library(lme4)
library(splines)
library(cmdstanr)
library(brms)
library(dplyr)
library(posterior)
library(slurmR)
set_cmdstan_path(path = "/gpfs/share/apps/cmdstan/2.25.0") 
mod <- cmdstan_model("/gpfs/home/dw2625/r/Review_sim/stan/time_varying_trt.stan")

# logistic <- function(a, x, x0, k,b) {
#   b + a / (1 + exp(k * (x-x0)))
# }

s_generate <- function(iter,  ncluster) {
  
  set.seed(iter)
  def <- defData(varname = "a", formula = 0, variance = 0.25)
  def <- defData(def, varname = "mu_b", formula = 0, dist = "nonrandom")
  def <- defData(def, varname = "s2_b", formula = 0.36, dist = "nonrandom")
  #A: trt for each cluster and time period
  #b: site-specific time period effect
  #t:exposure time
  deft <- defDataAdd(varname = "t", formula = "ifelse((k-startTrt)>0,(k-startTrt),0)", dist= "nonrandom")
 
  # Modify the treatment effect to increase, plateau, then decrease
  
  defteffect <- defCondition(
    condition = "t == 0",  
    formula = "0", 
    dist = "nonrandom")
  
  defteffect <- defCondition(defteffect,
                             condition = "t < 5 & t>0",  
                             formula = "5 / (1 + exp(-5 * (t-1)))", 
                             dist = "nonrandom")
  
  defteffect <- defCondition(defteffect,
                             condition = " t>=5",  
                             formula = "2.5 + 5 / (1 + exp(2 * (t-5)))", 
                             dist = "nonrandom")
  
  
  # Modify the outcome variable formula
  defOut <- defDataAdd(varname = "y", 
                       formula = "a + b - 0.05 * k^2 + t_effect * A", 
                       variance = 1)
  
  #--- add data generation code ---#
  ds <- genData(ncluster, def, id = "site")#site
  ds <- addPeriods(ds, ncluster+3+10, "site", perName = "k") #create time periods for each site
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
  dd <- addColumns(deft, dd)
  dd <- addCondition(defteffect, dd, newvar = "t_effect")
  dd <- addColumns(defOut, dd)
  dd[, normk := (k - min(k))/(max(k) - min(k))]#scale time period into range 0-1
  
  dd[] #  generated_data is a data.table
  
}

s_model <- function(train_data, mod) {
  set_cmdstan_path(path = "/gpfs/share/apps/cmdstan/2.25.0")

  t_ex_mx=max(train_data$t)
  #true value of TATE
  t_values <- 1:t_ex_mx
  true_eff = unique(train_data[t%in%t_values,t_effect])
  #true value of LTE
  ######Bayesian estimation of  TATE & LTE#############
  # Define the knots based on the training data
  knots_train <- quantile(train_data$k, probs=seq(0, 1, length=6)[-c(1, 6)])
  
  B_train <- predict(splines::bs(train_data$k, degree=3, knots=quantile(train_data$k, probs=seq(0, 1, length=6)[-c(1, 6)])))
  
  knots_t <- quantile(train_data$t, probs=seq(0, 1, length=6)[-c(1, 6)])
  
  B_t <- predict(splines::bs(train_data$t, degree=3, knots=quantile(train_data$t, probs=seq(0, 1, length=6)[-c(1, 6)])))
  #unique(B_t)
  # Compute the B-spline basis for the test set using the same knots from the training data
  B_test_t <- predict(splines::bs(c(0: t_ex_mx), degree=3, knots=knots_t))[-1,]#-1:remove t=0
  
  
  stan_data <- list(num_data = nrow(train_data),
                    num_basis = ncol(B_train),
                    num_t_ex = nrow(B_test_t),
                    B = t(B_train),
                    B_star = t(B_t),
                    B_test_star=t(B_test_t),
                    A = train_data$A,
                    y = train_data$y,
                    num_clusters = length(unique(train_data$site)),
                    cluster = train_data$site
                    
  )
  # Fit the Bayesian model
  fit <- mod$sample(data = stan_data,
                    refresh = 0, 
                    iter_warmup = 500,
                    iter_sampling = 2000,
                    parallel_chains = 4,
                    show_messages = TRUE) 
  
  diagnostics_df <- as_draws_df(fit$sampler_diagnostics())
  div <- sum(diagnostics_df[, 'divergent__'])
  num_max_tree_depth <- sum(diagnostics_df$treedepth__ >= 10)
  
  bayes_phi = fit$summary(variables="phi",
                          posterior::default_summary_measures()[1:3],
                          quantiles = ~ quantile(., probs = c(0.025, 0.975)))
  
  
  model_results <- data.table(true_eff,bayes_phi,
                              div,num_max_tree_depth)
  return(model_results)
  
  
}

s_single_rep <- function(iter, ncluster, mod) {
  
  train_data <- s_generate(iter, ncluster)
  
  model_results <- s_model(train_data, mod)
  
  return(model_results)
}

s_replicate <- function(iter, ncluster, mod) {
  model_results = s_single_rep(iter, ncluster, mod)
  return(data.table(iter=iter, ncluster=ncluster, model_results))
}


scenarios = expand.grid(ncluster=10)
i=1
nsim=150
ncluster= scenarios[i,"ncluster"]

# args <- commandArgs(trailingOnly = TRUE)
# iter <- as.numeric(args[1])  # Get the iteration number from the job array index
# 
# result <- s_replicate(iter=iter,  ncluster=ncluster, mod=mod)


sjob <- Slurm_lapply(1:150,
                     FUN=s_replicate,
                     ncluster = ncluster,
                     mod=mod,
                     njobs = 90,
                     tmp_path = "/gpfs/scratch/dw2625",
                     job_name = "sim15",
                     sbatch_opt = list(time = "48:00:00",partition = "cpu_medium", `mem-per-cpu` = "8G"),
                     export = c("s_generate","s_model","s_single_rep"),
                     plan = "wait",
                     overwrite=TRUE)
res <- Slurm_collect(sjob) # data is a list
#res<- site_plasma_all[lapply(site_plasma_all, function(x) length(x))>1] #filter out the error message
res <- rbindlist(res) # converting list to data.table

date_stamp <- gsub("-", "", Sys.Date())
dir.create(file.path("/gpfs/home/dw2625/r/Review_sim/results/", date_stamp), showWarnings = FALSE)
save(res, file = paste0("/gpfs/home/dw2625/r/Review_sim/results/", date_stamp, "/inc_dec_trt_m2_pointwise_ci.rda"))



