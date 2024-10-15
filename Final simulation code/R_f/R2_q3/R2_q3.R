#\item \comments{Figure 5. I suggest including results for the Bayesian cluster-specific time-varying effect model in figure 5 so we can see what the ``penalty'' is for including random treatment effects when none exist. Clearly this is a model choice the analyst must make without knowing the true answer.}

#\todo{use Figure 3 scenarios 2(time varying effect without cluster specific), fit the Bayesian cluster specific time-varying effect model.}
library(simstudy)
library(data.table)
library(mgcv)
library(lme4)
library(splines)
library(cmdstanr)
library(brms)
library(dplyr)
library(loo)
library(slurmR)
library(fda)
#references:https://www.rdatagen.net/post/2022-12-13-modeling-the-secular-trend-in-a-stepped-wedge-design/
#references:https://www.rdatagen.net/post/2022-11-01-modeling-secular-trend-in-crt-using-gam/
set_cmdstan_path(path = "/gpfs/share/apps/cmdstan/2.25.0") 
# Compile the Stan model
#mod <- cmdstan_model("./stepped_wedge_random_walk_prior_test.stan")
#mod <- cmdstan_model("./pen_spline_prior.stan")
mod <- cmdstan_model("/gpfs/home/dw2625/r/Review_sim/stan/pen_spline_prior.stan")

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
  true_tate = mean(unique(train_data[t%in%t_values,t_effect]))
  #true value of LTE
  true_lte = mean(train_data[t == t_ex_mx, t_effect])
  ######Bayesian estimation of  TATE & LTE#############
  # Define the knots based on the training data
  knots_train <- quantile(train_data$k, probs=seq(0, 1, length=6)[-c(1, 6)])
  #basis matrix for study time
  B_train <- predict(splines::bs(train_data$k, degree=3, knots=quantile(train_data$k, probs=seq(0, 1, length=6)[-c(1, 6)])))
  
  knots_t <- quantile(train_data$t, probs=seq(0, 1, length=6)[-c(1, 6)])
  #basis matrix for exposure time
  B_t <- predict(splines::bs(train_data$t, degree=3, knots=quantile(train_data$t, probs=seq(0, 1, length=6)[-c(1, 6)])))
  #unique(B_t)
  # Compute the B-spline basis for the test set using the same knots from the training data
  B_test_t <- predict(splines::bs(c(0: t_ex_mx), degree=3, knots=knots_t))[-1,]#-1:remove t=0
  
  
  # Create the penalty matrix using the same number of basis functions
  B_basis <- create.bspline.basis(rangeval=range(train_data$k), nbasis=7, norder=4)
  
  # Compute the penalty matrix for the second derivative
  P <- bsplinepen(B_basis, Lfdobj=2)
  
  B_basis_t <- create.bspline.basis(rangeval=range(train_data$t), nbasis=7, norder=4)
  
  # Compute the penalty matrix for the second derivative
  P_t <- bsplinepen(B_basis_t, Lfdobj=2)
  
  P <- P + diag(1e-6, nrow(P))
  P_t <- P_t + diag(1e-6, nrow(P_t))
  
  
  stan_data <- list(num_data = nrow(train_data),
                    num_basis = ncol(B_train),
                    num_t_ex = nrow(B_test_t),
                    B = t(B_train),
                    B_star = t(B_t),
                    B_test_star=t(B_test_t),
                    A = train_data$A,
                    y = train_data$y,
                    num_clusters = length(unique(train_data$site)),
                    cluster = train_data$site,
                    P = P,
                    P_t = P_t
                    
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
  bayes_tate = fit$summary(variables="phi_mean",
                           posterior::default_summary_measures()[1:3],
                           quantiles = ~ quantile(., probs = c(0.025, 0.975)),
                           posterior::default_convergence_measures())
  covered_bayes_tate =   (bayes_tate$`2.5%`<   true_tate  &   true_tate  < bayes_tate$`97.5%`)
  
  #LTE from Bayesian model
  bayes_lte = fit$summary(variables="phi_max",
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
  
  ############GAM estimation
  #The model includes the smooth for t only when A=1
  # fitgam <- mgcv::bam(y ~ s(t, by=A) +  s(k) + s(k, site, bs = "fs"), data = train_data, method="fREML")
  # newdata <- data.frame(t = t_values, A = 1,k=1,site=1)
  # predictions <- predict(fitgam, newdata, type = "terms")
  # smooth_t_values <- predictions[,"s(t):A"]
  # gam_tate <- mean(smooth_t_values)
  # gam_lte <- smooth_t_values[t_ex_mx]
  
  model_results <- data.table(true_tate,true_lte,bayes_tate,bayes_lte,covered_bayes_tate,covered_bayes_lte,
                              div,num_max_tree_depth,loo_data)
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
# res <- replicate(1, s_replicate(iter=1,ncluster=ncluster,
#                                 mod=mod
#                                     ))
# 
# res <- parallel::mclapply(
#   X = 1 : nsim, 
#   FUN = function(x) s_replicate(iter=x,coefA = coefA,ncluster=ncluster,
#                                 mod=mod),
#   mc.cores = 4)


job <- Slurm_lapply(
  X = 1:150,
  FUN =s_replicate,
  ncluster=ncluster,
  mod=mod,
  njobs = 90,
  mc.cores = 4L,
  job_name = "sim_1",
  tmp_path = "/gpfs/scratch/dw2625",
  plan = "wait",
  sbatch_opt = list(time = "48:00:00", partition = "cpu_medium", `mem-per-cpu` = "8G"),
  export = c("s_generate", "s_model","s_single_rep"),
  overwrite = TRUE
)

res <- Slurm_collect(job)
res <- rbindlist(res)

date_stamp <- gsub("-", "", Sys.Date())
dir.create(file.path("/gpfs/home/dw2625/r/Review_sim/results/", date_stamp), showWarnings = FALSE)
save(res, file = paste0("/gpfs/home/dw2625/r/Review_sim/results/", date_stamp, "/time_varying_gen_m3fit_update.rda"))




