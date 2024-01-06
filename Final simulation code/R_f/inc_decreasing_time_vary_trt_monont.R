library(simstudy)
library(data.table)
library(mgcv)
library(lme4)
library(splines)
library(cmdstanr)
library(brms)
library(dplyr)
library(loo)

#references:https://www.rdatagen.net/post/2022-12-13-modeling-the-secular-trend-in-a-stepped-wedge-design/
#references:https://www.rdatagen.net/post/2022-11-01-modeling-secular-trend-in-crt-using-gam/
set_cmdstan_path(path = "/gpfs/share/apps/cmdstan/2.25.0") 
# Compile the Stan model
#mod <- cmdstan_model("./stepped_wedge_random_walk_prior_test.stan")
#mod <- cmdstan_model("./stepped_wedge_time_vary_trt.stan")
mod <- cmdstan_model("/gpfs/data/troxellab/danniw/r/BS/vary_t_monotone_model.stan")


s_generate <- function(iter, ncluster) {
  
  set.seed(iter)
  
  def <- defData(varname = "a", formula = 0, variance = 0.25)
  def <- defData(def, varname = "mu_b", formula = 0, dist = "nonrandom")
  def <- defData(def, varname = "s2_b", formula = 0.36, dist = "nonrandom")
  #A: trt for each cluster and time period
  #b: site-specific time period effect
  #t:exposure time
  deft <- defDataAdd(varname = "t", formula = "max(k-startTrt,0) ", dist= "nonrandom")
  
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
  true_tate = mean(unique(train_data[t%in%t_values,t_effect]))
  #true value of LTE
  true_lte = unique(train_data[t==t_ex_mx ,"t_effect"])
  ######Bayesian estimation of  TATE & LTE#############
  # Define the knots based on the training data
  knots_train <- quantile(train_data$k, probs=seq(0, 1, length=6)[-c(1, 6)])
  
  B_train <- predict(splines::bs(train_data$k, degree=3, knots=quantile(train_data$k, probs=seq(0, 1, length=6)[-c(1, 6)])))
  
  
  monot_data <- list(num_data = nrow(train_data),
                     num_basis = ncol(B_train),
                     num_t_ex= t_ex_mx,
                     t_ex = train_data$t,
                     B = t(B_train),
                     c=rep(1,t_ex_mx),
                     A = train_data$A,
                     y = train_data$y,
                     num_sites = length(unique(train_data$site)),
                     site = train_data$site
                     
  )
  # Fit the Bayesian model
  fit <- mod$sample(data = monot_data,
                    refresh = 0, 
                    iter_warmup = 500,
                    iter_sampling = 2000,
                    parallel_chains = 4,
                    show_messages = TRUE) 
  
  diagnostics_df <- as_draws_df(fit$sampler_diagnostics())
  div <- sum(diagnostics_df[, 'divergent__'])
  num_max_tree_depth <- sum(diagnostics_df$treedepth__ >= 10)
  bayes_tate = fit$summary(variables="phi_mean",
                          posterior::default_summary_measures()[1:3],
                          quantiles = ~ quantile(., probs = c(0.025, 0.975)),
                          posterior::default_convergence_measures())
  covered_bayes_tate =   (bayes_tate$`2.5%`<   true_tate  &   true_tate  < bayes_tate$`97.5%`)
  
  #LTE from Bayesian model
  bayes_lte = fit$summary(variables="delta",
                          posterior::default_summary_measures()[1:3],
                          quantiles = ~ quantile(., probs = c(0.025, 0.975)),
                          posterior::default_convergence_measures())
  covered_bayes_lte =   (bayes_lte$`2.5%`<   true_lte  &   true_lte  < bayes_lte$`97.5%`)
  
  loo_result <- fit$loo()
  
  loo_data <- data.table(
    elpd_loo = loo_result$estimate['elpd_loo', 'Estimate'],
    elpd_se = loo_result$estimate['elpd_loo', 'SE'],
    p_loo = loo_result$estimate['p_loo', 'Estimate'],
    p_loo_se = loo_result$estimate['p_loo', 'SE'],
    looic = loo_result$estimate['looic', 'Estimate'],
    looic_se = loo_result$estimate['looic', 'SE'],
    mcse_elpd_loo = loo_result$mcse['elpd_loo', 'mcse']
  )
  
  
  model_results <- data.table(true_tate,true_lte,bayes_tate,bayes_lte,covered_bayes_tate,covered_bayes_lte,
                              div,num_max_tree_depth, loo_data)
  return(model_results)
  
  
}

s_single_rep <- function(iter, ncluster, mod) {
  
  train_data <- s_generate(iter, ncluster)
  
  model_results <- s_model(train_data,  mod)
  
  return(model_results)
}

s_replicate <- function(iter,  ncluster, mod) {
  model_results = s_single_rep(iter,  ncluster, mod)
  return(data.table(iter=iter, ncluster=ncluster, model_results))
}


scenarios = expand.grid(ncluster=10)
i=1

ncluster= scenarios[i,"ncluster"]

args <- commandArgs(trailingOnly = TRUE)
iter <- as.numeric(args[1])  # Get the iteration number from the job array index

result <- s_replicate(iter=iter, ncluster=ncluster, mod=mod)


# sjob <- Slurm_lapply(1:nsim,
#                FUN=s_replicate,
#                coefA = coefA,
#                ncluster = ncluster,
#                mod=mod,
#                njobs = 30,
#                tmp_path = "/gpfs/scratch/dw2625",
#                job_name = "BS_140",
#                sbatch_opt = list(time = "24:00:00",partition = "cpu_short", `mem-per-cpu` = "10G"),
#                export = c("s_define","s_generate","s_model","s_single_rep"),
#                plan = "wait",
#                overwrite=TRUE)
# res <- Slurm_collect(sjob) # data is a list
# #res<- site_plasma_all[lapply(site_plasma_all, function(x) length(x))>1] #filter out the error message
# res <- rbindlist(res) # converting list to data.table

#date_stamp <- gsub("-", "", Sys.Date())
dir.create(file.path("/gpfs/data/troxellab/danniw/scratch/monotone_inc_decreasing_time_vary_trt"), showWarnings = FALSE)

# Write the result to a file
output_filename <- paste0("output_iter_", iter, ".RData")
#save(result, file=output_filename)

save(result, file = paste0("/gpfs/data/troxellab/danniw/scratch/monotone_inc_decreasing_time_vary_trt/", output_filename))

#result <- cbind(res)

#save(res, file = paste0("./scenarios_fixed_coefA",coefA,"ncluster",ncluster,".rda"))


