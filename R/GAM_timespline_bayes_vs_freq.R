#Data generation
#In this data generation process, the time effect will not be explicitly smooth, but the underlying covariance structure used to generate the period effects will induce some level of smoothness. 
#references:https://www.rdatagen.net/post/2022-12-13-modeling-the-secular-trend-in-a-stepped-wedge-design/
#references:https://www.rdatagen.net/post/2022-11-01-modeling-secular-trend-in-crt-using-gam/
#brms:https://github.com/paul-buerkner/brms
library(simstudy)
library(data.table)
library(mgcv)
library(lme4)
library(splines)
library(brms)
library(cmdstanr)
library(slurmR)


set_cmdstan_path(path = "/gpfs/share/apps/cmdstan/2.25.0")
mod <- cmdstan_model("/gpfs/data/troxellab/danniw/r/time_spline.stan");
#mod <- cmdstan_model("/gpfs/data/troxellab/danniw/r/cluster_specific_time_spline.stan");
#mod <- cmdstan_model("./time_spline.stan");
#mod <- cmdstan_model("./cluster_specific_time_spline.stan");
s_define <- function() {
  #cluster-specific intercept
  def <- defData(varname = "a", formula = 0, variance = 9)
  #mean of cluster-time effect
  def <- defData(def, varname = "mu_b", formula = 0, dist = "nonrandom")
  #variance of cluster-time effect
  def <- defData(def, varname = "s2_b", formula = 9, dist = "nonrandom")
  #b:cluster-time effect
  #A: trt for each cluster and time period
  defOut <- defDataAdd(varname = "y", formula = "a + b + 5 * A", variance = 4)
  
  return(list(def = def, defOut = defOut)) 
}

s_generate <- function(list_of_defs) {
  
  list2env(list_of_defs, envir = environment())
  
  #--- add data generation code ---#
  #24 sites in total
  ds <- genData(24, def, id = "site")#24 site
  ds <- addPeriods(ds, 25, "site", perName = "k") #create 25 periods for each site
  #cluster-time effect, a vector of site-specific time period effect,
  #assume each cluster ~MVN(0,\Sigma_b)
  #b:period-specific effects for each site
  ds <- addCorGen(dtOld = ds, idvar = "site", 
                  rho = 0.8, corstr = "ar1",#covariance matrix 
                  dist = "normal", param1 = "mu_b", param2 = "s2_b", cnames = "b")
  #assign the treatment status based on the stepped-wedge design
  #per cluster trt change per wave
  ds <- trtStepWedge(ds, "site", nWaves = 24, lenWaves = 1, startPer = 1, 
                     grpName = "A", perName = "k")
  ds$site <- as.factor(ds$site)
  #30 individuals per site per period and generate each individual-level outcome
  dd <- genCluster(ds, "timeID", numIndsVar = 30, level1ID = "id")
  dd <- addColumns(defOut, dd)
  dd[, normk := (k - min(k))/(max(k) - min(k))]#scale time period into range 0-1
  
  dd[] #  generated_data is a data.table
}

# dp <- dd[, .(avg = mean(y)), keyby = .(A, site, k)]
# 
# ggplot(data = dp, aes(x = k, y = avg)) +
#   geom_line(aes(group = site, color = factor(A))) +
#   scale_color_manual(values = c("#a53e3e", "#3ea5a5"),
#                      guide = guide_legend(reverse = TRUE),
#                      labels = c("no treatment", "treatment")) +
#   ylab("average Y") +
#   xlab("period (k)") +
#   theme(panel.grid = element_blank(),
#         legend.title = element_blank()) 

# ggplot(data = dd, aes(x = k, y = y)) +
#   geom_point(aes(color = factor(A, labels = c("Control", "Intervention"))), 
#              size = 0.1) +
#   scale_color_manual(values = c("#d07b7c", "#7ba7d0")) +
#   facet_wrap(~site, ncol = 8) +
#   theme(panel.grid = element_blank(),
#         legend.title = element_blank(),
#         axis.text.y = element_text(size = 7),
#         strip.text = element_text(size = 8),
#         axis.text.x = element_blank(),
#         axis.ticks.x = element_blank()) +
#   xlim(c(0, 24)) +
#   guides(color = guide_legend(override.aes = list(size = 2)))
s_model <- function(dd, mod) {
  set_cmdstan_path(path = "/gpfs/share/apps/cmdstan/2.25.0")
  fitgam <- gam(y ~ A + s(k, site, bs = "fs", k = 5), data = dd, method="REML")
  res_fitgam <- c(summary(fitgam)$p.coeff["A"], summary(fitgam)$se["A"])
  #fit a Bayesian version
  #stan code: time_spline.stan
  # df: Compute B-spline basis matrix
  # Decide on the number of basis functions
  N <- nrow(dd)
  J <- length(unique(dd$site))
  B <- splines::bs(dd$k, df = 5, intercept = TRUE)
  # matplot(dd$timeID, B, type = 'l', lty = 1, 
  #         xlab = "Time", ylab = "B-spline Basis", 
  #         main = "B-spline Basis Functions")
  
  J = length(unique(dd$site))
  
  K_bspline = ncol(B)
  
  studydata <- list(
    N = N, J=J, K_bspline=K_bspline,
    B=B, Y=dd$y, cluster = as.integer(dd$site), A=dd$A)
  
  fit <- mod$sample(
    data = studydata,
    refresh = 0,
    chains = 4L,
    parallel_chains = 4L,
    iter_warmup = 500,
    iter_sampling = 2500,
    show_messages = FALSE)
  
  diagnostics_df <- as_draws_df(fit$sampler_diagnostics())
  div <- sum(diagnostics_df[, 'divergent__'])
  tree_hit <- sum(diagnostics_df$treedepth__ == 10)
  res_bayesgam = fit$summary(variables="delta",
                          posterior::default_summary_measures()[1:3],
                          quantiles = ~ quantile(., probs = c(0.025, 0.975)),
                          posterior::default_convergence_measures())
  
  model_results <- data.table(t(res_fitgam), res_bayesgam,div)
  setnames(model_results, c("est.freq", "se.freq", "variable","est.mean.bayes",
                            "est.med.bayes","est.sd.bayes",
                            "lowci.bayes",
                            "upci.bayes","rhat","ess_bulk","ess_tail",
                            "div"))
  return(model_results) # model_results is a data.table
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


job <- Slurm_lapply(X=1:200,
                    FUN = s_replicate,
                    mod=mod,
                    njobs = 90,
                    mc.cores = 4L,
                    job_name = "Stw_4",
                    tmp_path = "/gpfs/data/troxellab/danniw/scratch",
                    plan = "wait",
                    sbatch_opt = list(time = "12:00:00", partition = "cpu_short", `mem-per-cpu` = "8G"),
                    export = c("s_define","s_generate","s_model","s_single_rep"),
                    overwrite = TRUE)

results <- Slurm_collect(job)
results_agg <- rbindlist(results)

date_stamp <- gsub("-", "", Sys.Date())
dir.create(file.path("/gpfs/data/troxellab/danniw/data/Stw/", date_stamp), showWarnings = FALSE)
save(results_agg, file = paste0("/gpfs/data/troxellab/danniw/data/Stw/", date_stamp, "/GAM_timespline_bayes_freq_v3.rda"))
