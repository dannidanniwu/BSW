#Data generation
#In this data generation process, the time effect will not be explicitly smooth, but the underlying covariance structure used to generate the period effects will induce some level of smoothness. 
#references:https://www.rdatagen.net/post/2022-12-13-modeling-the-secular-trend-in-a-stepped-wedge-design/
#references:https://www.rdatagen.net/post/2022-11-01-modeling-secular-trend-in-crt-using-gam/
library(simstudy)
library(data.table)
library(mgcv)
library(lme4)
library(splines)
library(brms)
library(rstan)
library(slurmR)

s_define <- function() {
  #cluster-specific intercept
  def <- defData(varname = "a", formula = 0, variance = 9)
  #mean of cluster-time effect
  def <- defData(def, varname = "mu_b", formula = 0, dist = "nonrandom")
  #variance of cluster-time effect
  def <- defData(def, varname = "s2_b", formula = 9, dist = "nonrandom")
  #b:cluster-time effect
  #A: trt for each cluster and time period
  defOut <- defDataAdd(varname = "y", formula = "a + b + 5 * A", variance = 25)
  
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
s_model <- function(dd) {
  
  fitgam <- gam(y ~ A + s(k, site, bs = "fs", k = 5), data = dd, method="REML")
  res_fitgam <- c(summary(fitgam)$p.coeff["A"], summary(fitgam)$se["A"])
  
  #fit a Bayesian version
  #https://fromthebottomoftheheap.net/2018/04/21/fitting-gams-with-brms/
  #res_bayesgam$model: can be used to check the stan code
  #every parameter is summarized using the mean (Estimate) and the standard deviation
  #(Est.Error) of the posterior distribution as well as two-sided 95% credible intervals
  #(l-95% CI and u-95% CI) based on quantiles
  bayesgam = brm(y ~ A +  s(k, site, bs = "fs", k = 5), data = dd, family = gaussian(), #cores = 4,
            iter = 2500, warmup = 500,  refresh = 0,
            control = list(adapt_delta = 0.9))
  sparams = get_sampler_params(bayesgam$fit, inc_warmup=FALSE)
  div = sum(sapply(sparams, function(x) sum(x[, 'divergent__'])))
  res_bayesgam <- c(summary(bayesgam)$fixed$Estimate[2], summary(bayesgam)$fixed$`l-95% CI`[2],
                    summary(bayesgam)$fixed$`u-95% CI`[2],median(extract(bayesgam$fit, pars = 'b_A')[[1]]), 
                    summary(bayesgam)$fixed$Rhat[2],div)
  
  model_results <- data.table(t(res_fitgam), t(res_bayesgam))
  setnames(model_results, c("est.freq", "se.freq", "est.mean.bayes", "lowci.bayes",
                            "upci.bayes","est.med.bayes.gam","rhat","div"))
  return(model_results) # model_results is a data.table
}

s_single_rep <- function(list_of_defs) {
  
  generated_data <- s_generate(list_of_defs)
  model_results <- s_model(generated_data)
  
  return(model_results)
}

s_replicate <- function(iter) {
  list_of_defs = s_define()
  model_results = s_single_rep(list_of_defs)
  return(data.table(iter=iter, model_results))
}


job <- Slurm_lapply(X=1:200,
                    FUN = s_replicate,
                    njobs = 90,
                    mc.cores = 4L,
                    job_name = "Stw_2",
                    tmp_path = "/gpfs/data/troxellab/danniw/scratch",
                    plan = "wait",
                    sbatch_opt = list(time = "24:00:00", partition = "cpu_short", `mem-per-cpu` = "8G"),
                    export = c("s_define","s_generate","s_model","s_single_rep"),
                    overwrite = TRUE)

results <- Slurm_collect(job)
results_agg <- rbindlist(results)

date_stamp <- gsub("-", "", Sys.Date())
dir.create(file.path("/gpfs/data/troxellab/danniw/data/Stw/", date_stamp), showWarnings = FALSE)
save(results_agg, file = paste0("/gpfs/data/troxellab/danniw/data/Stw/", date_stamp, "/GAM_bayes_freq.rda"))
