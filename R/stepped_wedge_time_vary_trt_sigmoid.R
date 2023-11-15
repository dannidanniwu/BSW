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
#sigmoid_mod <- cmdstan_model("./sigmoid_effect_curve_model.stan")
sigmoid_mod <- cmdstan_model("/gpfs/home/dw2625/r/BS/sigmoid_effect_curve_model.stan")
s_define <- function() {
  #cluster-specific intercept
  def <- defData(varname = "a", formula = 0, variance = 0.25)
  def <- defData(def, varname = "mu_b", formula = 0, dist = "nonrandom")
  def <- defData(def, varname = "s2_b", formula = 0.36, dist = "nonrandom")
  #A: trt for each cluster and time period
  #b: site-specific time period effect
  #t:exposure time
  deft <- defDataAdd(varname = "t", formula = "max(k-startTrt,0) ", dist= "nonrandom")
  
  defteffect <- defDataAdd(varname = "t_effect", formula = " (0 * (t == 0) +(t!=0)*(..coefA * 1 / (1 + 2*exp(-t))))", dist = "nonrandom")
  
  defOut <- defDataAdd(varname = "y", formula = " a + b - 0.05 * k^2 + (0 * (t == 0) +(t != 0)*(..coefA * 1 / (1 + 2*exp(-t))))*A", variance = 1)
  
  return(list(def = def, deft=deft, defteffect=defteffect, defOut = defOut)) 
}

s_generate <- function(iter, coefA, ncluster, list_of_defs) {
  
  set.seed(iter)
  list2env(list_of_defs, envir = environment())
  
  #--- add data generation code ---#
  ds <- genData(ncluster, def, id = "site")#site
  ds <- addPeriods(ds, ncluster+3+5, "site", perName = "k") #create time periods for each site
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
  dd <- addColumns(defteffect, dd)
  dd <- addColumns(defOut, dd)
  dd[, normk := (k - min(k))/(max(k) - min(k))]#scale time period into range 0-1
  
  dd[] #  generated_data is a data.table
  
}

s_model <- function(train_data, coefA, sigmoid_mod) {
  set_cmdstan_path(path = "/gpfs/share/apps/cmdstan/2.25.0")

  t_ex_mx=max(train_data$t)
  #true value of TATE
  t_values <- 1:t_ex_mx
  true_tate = sum(coefA * 1 / (1 + 2*exp(-t_values)))/t_ex_mx
  #true value of LTE
  true_lte = coefA * 1 / (1 + 2*exp(-t_ex_mx))
  ######Bayesian estimation of  TATE & LTE#############
  knots_train <- quantile(train_data$k, probs=seq(0, 1, length=6)[-c(1, 6)])
  
  B_train <- predict(splines::bs(train_data$k, degree=3, knots=quantile(train_data$k, probs=seq(0, 1, length=6)[-c(1, 6)])))
  
  #############Bayesian monontone_effect_curve_model#################
  sigmoid_data <- list(num_data = nrow(train_data),
                    num_basis = ncol(B_train),
                    num_t_ex= t_ex_mx,
                    t_ex = train_data$t,
                    B = t(B_train),
                    A = train_data$A,
                    y = train_data$y,
                    num_sites = length(unique(train_data$site)),
                    site = train_data$site
                    
  )
  fit_sigmoid <- sigmoid_mod$sample(data = sigmoid_data,
                    refresh = 0, 
                    iter_warmup = 500,
                    iter_sampling = 2000,
                    parallel_chains = 4,
                    show_messages = TRUE) 
  
  sigmoid_diagnostics_df <- as_draws_df(fit_sigmoid$sampler_diagnostics())
  div_sigmoid <- sum(sigmoid_diagnostics_df[, 'divergent__'])
  num_max_tree_depth_sigmoid <- sum(sigmoid_diagnostics_df$treedepth__ >= 10)
  
  bayes_sigmoid_tate = fit_sigmoid$summary(variables="phi_mean",
                           posterior::default_summary_measures()[1:3],
                           quantiles = ~ quantile(., probs = c(0.025, 0.975)),
                           posterior::default_convergence_measures())
  covered_bayes_sigmoid_tate =   (bayes_sigmoid_tate$`2.5%`<   true_tate  &   true_tate  < bayes_sigmoid_tate$`97.5%`)
  
  
  bayes_sigmoid_lte = fit_sigmoid$summary(variables="phi_max",
                           posterior::default_summary_measures()[1:3],
                           quantiles = ~ quantile(., probs = c(0.025, 0.975)),
                           posterior::default_convergence_measures())
  
  
  
  model_results <- data.table(true_tate,true_lte,bayes_sigmoid_tate,bayes_sigmoid_lte,
                              covered_bayes_sigmoid_tate,covered_bayes_sigmoid_lte,div_sigmoid,
                              num_max_tree_depth_sigmoid)
  return(model_results)
  
  
}

s_single_rep <- function(iter, coefA, ncluster, list_of_defs, sigmoid_mod) {
  
  train_data <- s_generate(iter, coefA, ncluster, list_of_defs)
  
  model_results <- s_model(train_data, coefA, sigmoid_mod )
  
  return(model_results)
}

s_replicate <- function(iter, coefA, ncluster, sigmoid_mod) {
  list_of_defs = s_define()
  model_results = s_single_rep(iter, coefA, ncluster,list_of_defs, sigmoid_mod)
  return(data.table(iter=iter, coefA = coefA , ncluster=ncluster, model_results))
}


scenarios = expand.grid(coefA= 5,ncluster=10)
i=1
nsim=150
coefA = scenarios[i,"coefA"]
ncluster= scenarios[i,"ncluster"]
# res <- replicate(1, s_replicate(iter=1,coefA = coefA,ncluster=ncluster,
#                                 sigmoid_mod =sigmoid_mod 
#                                     ))
# 
res <- parallel::mclapply(
    X = 1 : nsim,
    FUN = function(x) s_replicate(iter=x,coefA = coefA,ncluster=ncluster,
                                  sigmoid_mod =sigmoid_mod ),
    mc.cores = 4)

# 
# 
# 
# sjob <- Slurm_lapply(1:nsim,
#                FUN=s_replicate,
#                coefA = coefA,
#                ncluster = ncluster,
#                sigmoid_mod =sigmoid_mod ,
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

date_stamp <- gsub("-", "", Sys.Date())
dir.create(file.path("/gpfs/home/dw2625/r/BS/", date_stamp), showWarnings = FALSE)
save(res, file = paste0("/gpfs/home/dw2625/r/BS/", date_stamp, "/sigmoid_time_vary_trt",coefA,"ncluster",ncluster,".rda"))

#result <- cbind(res)

#save(res, file = paste0("./scenarios_fixed_coefA",coefA,"ncluster",ncluster,".rda"))


