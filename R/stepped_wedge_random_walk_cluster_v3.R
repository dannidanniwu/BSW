library(simstudy)
library(data.table)
library(mgcv)
library(lme4)
library(splines)
library(cmdstanr)
library(brms)
library(slurmR)
library(dplyr)
#references:https://www.rdatagen.net/post/2022-12-13-modeling-the-secular-trend-in-a-stepped-wedge-design/
#references:https://www.rdatagen.net/post/2022-11-01-modeling-secular-trend-in-crt-using-gam/
set_cmdstan_path(path = "/gpfs/share/apps/cmdstan/2.25.0") 
# Compile the Stan model
#mod <- cmdstan_model("./stepped_wedge_random_walk.stan")
mod <- cmdstan_model("/gpfs/data/troxellab/danniw/r/BS/stepped_wedge_random_walk.stan")

s_define <- function() {
  #cluster-specific intercept
  def <- defData(varname = "a", formula = 0, variance = 0.25)
  def <- defData(def, varname = "mu_b", formula = 0, dist = "nonrandom")
  def <- defData(def, varname = "s2_b", formula = 0.4, dist = "nonrandom")
  
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
  set_cmdstan_path(path = "/gpfs/share/apps/cmdstan/2.25.0")
  # #######Fitting the GAM model with default penalization 
  fitsat <- lmer(y ~ A + as.factor(k) + (as.factor(k)|site), data = train_data)
  coef_A <- fixef(fitsat)["A"]
  se_A <- summary(fitsat)$coefficients["A", "Std. Error"]
  res_fitsat <- c(coef_A, se_A)
  range_sat <-   res_fitsat[1] + c(-1,1) * 1.96 *   res_fitsat[2]
  
  fitnotime <- lm(y ~ A, data = train_data)
  res_fitnotime <- c(coef(fitnotime)["A"], summary(fitnotime)$coefficients["A", "Std. Error"])
  range_notime <-   res_fitnotime[1] + c(-1,1) * 1.96 *   res_fitnotime[2]
  
  fitlntime <- lmer(y ~ A + k +  (k|site), data = train_data)
  coef_A <- fixef(fitlntime)["A"]
  se_A <- summary(fitlntime)$coefficients["A", "Std. Error"]
  res_lntime <- c(coef_A, se_A)
  range_lntime <-   res_lntime[1] + c(-1,1) * 1.96 *   res_lntime[2]
  
  
  fitgam <- mgcv::bam(y ~ A +  s(k) + s(k, site, bs = "fs"), data = train_data, method="fREML")
  res_fitgam <- c(summary(fitgam)$p.coeff["A"], summary(fitgam)$se["A"])
  range <-   res_fitgam[1] + c(-1,1) * 1.96 *   res_fitgam[2]
  
  # Fitting a more complex  model using gam
  fitgam2 <- mgcv::bam(y ~ A +  s(k) + s(site, A, bs = "re") + s(k, site,  bs = "fs"), data = train_data, method="fREML")
  res_fitgam2 <- c(summary(fitgam2)$p.coeff["A"], summary(fitgam2)$se["A"])
  range2 <-   res_fitgam2[1] + c(-1,1) * 1.96 *   res_fitgam2[2]
  
  # Define the knots based on the training data
  knots_train <- quantile(train_data$k, probs=seq(0, 1, length=6)[-c(1, 6)])
  
  B_train <- predict(splines::bs(train_data$k, degree=3, knots=quantile(train_data$k, probs=seq(0, 1, length=6)[-c(1, 6)])))
  colnames(B_train) <- paste0("Bspline_", 1:ncol(B_train))
  
  
  stan_data <- list(num_data = nrow(train_data),
                    num_basis = ncol(B_train),
                    B = t(B_train),
                    A = train_data$A,
                    y = train_data$y,
                    num_sites = length(unique(train_data$site)),
                    sigma_beta_A = 5,
                    site = train_data$site
                    
  )
  #mod <- cmdstan_model("./stepped_wedge_random_walk.stan")
  # Fit the Bayesian model
  fit <- mod$sample(data = stan_data,
                    refresh = 0, 
                    iter_warmup = 500,
                    iter_sampling = 2000,
                    show_messages = TRUE) 
  
  diagnostics_df <- as_draws_df(fit$sampler_diagnostics())
  div <- sum(diagnostics_df[, 'divergent__'])
  bayes_gam = fit$summary(variables="beta_A",
                          posterior::default_summary_measures()[1:3],
                          quantiles = ~ quantile(., probs = c(0.025, 0.975)),
                          posterior::default_convergence_measures())
  
  covered_bayes =   (bayes_gam$`2.5%`< coefA & coefA < bayes_gam$`97.5%`)
  
  stan_data_sd10 <- list(num_data = nrow(train_data),
                    num_basis = ncol(B_train),
                    B = t(B_train),
                    A = train_data$A,
                    y = train_data$y,
                    num_sites = length(unique(train_data$site)),
                    sigma_beta_A = 10,
                    site = train_data$site
                    
  )
  
  fit_sd10 <- mod$sample(data = stan_data,
                    refresh = 0, 
                    iter_warmup = 500,
                    iter_sampling = 2000,
                    show_messages = TRUE) 
  
  diagnostics_df_sd10 <- as_draws_df(fit_sd10$sampler_diagnostics())
  div_sd10 <- sum(diagnostics_df_sd10[, 'divergent__'])
  bayes_gam_sd10 = fit_sd10$summary(variables="beta_A",
                          posterior::default_summary_measures()[1:3],
                          quantiles = ~ quantile(., probs = c(0.025, 0.975)),
                          posterior::default_convergence_measures())
  
  bayes_gam_sd10 =bayes_gam_sd10[,-1]
  covered_bayes_sd10 =   (bayes_gam_sd10$`2.5%`< coefA & coefA < bayes_gam_sd10$`97.5%`)
  
  #Fit a frequentist linear model with the same basis as the Bayesian model, but no penalization
  # Incorporating the B-spline basis into the data
  ds_with_bspline <- cbind(train_data,  B_train)
  # Fitting the model using gam
  bspline_terms <- paste(colnames(B_train), collapse = " + ")
  formula_str <- paste("y ~ A + s(site, A, bs = 're') +", bspline_terms)
  fitgam3 <- mgcv::bam(as.formula(formula_str), data = ds_with_bspline, method="fREML")
  res_fitgam3 <- c(summary(fitgam3)$p.coeff["A"], summary(fitgam3)$se["A"])
  range3 <-   res_fitgam3[1] + c(-1,1) * 1.96 *   res_fitgam3[2]
  
  
  model_results <- data.table(est_gam_freq= res_fitgam[1], se_gam_freq=res_fitgam[2], 
                              gam_lowci=range[1], gam_upci=range[2], bayes_gam,div, bayes_gam_sd10,
                              div_sd10,
                              est_gam_freq_rdn = res_fitgam2[1], se_gam_freq_rdn=res_fitgam2[2], 
                              gam_rdn_lowci=range2[1], gam_rdn_upci=range2[2],
                              est_gam_freq_np = res_fitgam3[1], se_gam_freq_np=res_fitgam3[2], 
                              gam_np_lowci=range3[1], gam_np_upci=range3[2],
                              covered_gam_freq=(range[1] < coefA & coefA < range[2]),
                              covered_bayes,covered_bayes_sd10,
                              covered_gam_rdn_freq=(range2[1] < coefA & coefA < range2[2]),
                              covered_gam_np_freq=(range3[1] < coefA & coefA < range3[2]),
                              est_sat = res_fitsat[1],
                              se_sat = res_fitsat[2],
                              sat_lowci = range_sat[1],
                              sat_upci = range_sat[2],
                              est_notime = res_fitnotime[1],
                              se_notime = res_fitnotime[2],
                              notime_lowci = range_notime[1],
                              notime_upci = range_notime[2],
                              est_lntime = res_lntime[1],
                              se_lntime = res_lntime[2],
                              lntime_lowci = range_lntime[1],
                              lntime_upci = range_lntime[2],
                              covered_sat = res_fitsat[1] < coefA & coefA < res_fitsat[2],
                              covered_notime = res_fitnotime[1] < coefA & coefA < res_fitnotime[2],
                              covered_lntime = res_lntime[1] < coefA & coefA < res_lntime[2]
                              
  ) 
  
  setnames(model_results, c("est_gam_freq", "se_gam_freq", "lowci_freq", "upci_freq", "variable", "est_mean_bayes",
                            "est_med_bayes", "est_sd_bayes", "lowci_bayes", "upci_bayes", "rhat", "ess_bulk", "ess_tail", "div", 
                            "est_mean_bayes10",
                            "est_med_bayes10","est_sd_bayes10",
                            "lowci_bayes10",
                            "upci_bayes10","rhat_bayes10","ess_bulk_bayes10","ess_tail_bayes10",
                            "div_bayes10",
                            "est_gam_rdn_freq", "se_gam_rdn_freq", "lowci_freq_rdn", "upci_freq_rdn", "est_gam_np_freq", "se_gam_np_freq", 
                            "lowci_freq_np", "upci_freq_np", "covered_freq", "covered_bayes","covered_bayes10", "covered_gam_rdn_freq",
                            "covered_gam_np_freq",
                            "est_sat", "se_sat", "sat_lowci", "sat_upci", "est_notime", "se_notime", "notime_lowci", "notime_upci", 
                            "est_lntime", "se_lntime", "lntime_lowci", "lntime_upci","covered_sat",
                            "covered_notime","covered_lntime"))
  
  model_results <- model_results%>%
    mutate(across(-c(variable,covered_freq, covered_bayes,covered_bayes10,covered_gam_rdn_freq, covered_gam_np_freq), round, 3))
  
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

scenarios = expand.grid(coefA=seq(0, 1.25, length.out = 6),ncluster=10)


i=6
coefA = scenarios[i,"coefA"]
ncluster= scenarios[i,"ncluster"]
# res <- replicate(1, s_replicate(iter=1,coefA = coefA,ncluster=ncluster,
#                                 mod=mod
#                                     ))

sjob <- Slurm_lapply(1:160,
               FUN=s_replicate,
               coefA = coefA,
               ncluster = ncluster,
               mod=mod,
               njobs = 90,
               tmp_path = "/gpfs/data/troxellab/danniw/scratch",
               job_name = "BS_117",
               sbatch_opt = list(time = "12:00:00",partition = "cpu_short", `mem-per-cpu` = "8G"),
               export = c("s_define","s_generate","s_model","s_single_rep"),
               plan = "wait",
               overwrite=TRUE)
res <- Slurm_collect(sjob) # data is a list
#res<- site_plasma_all[lapply(site_plasma_all, function(x) length(x))>1] #filter out the error message
res <- rbindlist(res) # converting list to data.table

date_stamp <- gsub("-", "", Sys.Date())
dir.create(file.path("/gpfs/data/troxellab/danniw/r/BS/", date_stamp), showWarnings = FALSE)
save(res, file = paste0("/gpfs/data/troxellab/danniw/r/BS/", date_stamp, "/scenarios",i,".rda"))




