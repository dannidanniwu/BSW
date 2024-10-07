library(simstudy)
library(data.table)
library(mgcv)
library(lme4)
library(splines)
library(cmdstanr)
library(brms)
library(dplyr)
library(loo)



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

s_model <- function(train_data,  mod) {
  
  t_ex_mx=max(train_data$t)
  #true value of TATE
  t_values <- 1:t_ex_mx
  true_eff = unique(train_data[t%in%t_values,t_effect])
  true_tate = mean(unique(train_data[t%in%t_values,t_effect]))
  #true value of LTE
  true_lte = unique(train_data[t==t_ex_mx ,"t_effect"])
  ######Bayesian estimation of  TATE & LTE#############
  # Define the knots based on the training data
  knots_train <- quantile(train_data$k, probs=seq(0, 1, length=6)[-c(1, 6)])
  
  boundary_knots <- range(train_data$k)
  
  B_train_ns <- splines::ns(train_data$k, knots=knots_train, Boundary.knots=boundary_knots)
  
  
  knots_t <- quantile(train_data$t, probs=seq(0, 1, length=6)[-c(1, 6)])
  
  t_boundary_knots <- range(train_data$t)
  
  B_t <- splines::ns(train_data$t, knots=knots_t, Boundary.knots=t_boundary_knots)
  
  
  fit <- lmer(y ~ B_train_ns + B_t:A + (1|site), data= train_data)
  
  
  coef_fixed <- fixef(fit)
  # Coefficients for the time effect (B_t:A)
  coef_time_effect <- coef_fixed[grep(":A", names(coef_fixed))]
  
  
  B_test_t <- splines::ns(c(1: t_ex_mx), knots=knots_t, Boundary.knots=t_boundary_knots)
  
  predicted_tau <- B_test_t[,1:length(coef_time_effect)] %*%  coef_time_effect
  
  ncs_tate = mean(predicted_tau)
  ncs_lte =  predicted_tau[t_ex_mx,] 
  # Step 1: Extract the variance-covariance matrix of the fixed effects
  vcov_matrix <- vcov(fit)
  
  se <- sqrt(diag(B_test_t[,1:length(coef_time_effect)] %*% vcov_matrix[grep(":A", names(coef_fixed)), grep(":A", names(coef_fixed))] %*% t(B_test_t[,1:length(coef_time_effect)])))
  
  lower_bound <- predicted_tau - 1.96 * se
  upper_bound <- predicted_tau + 1.96 * se
  
  results <- data.frame(
    time = 1:t_ex_mx,
    true_eff = true_eff,
    predicted = predicted_tau,
    lower_ci = lower_bound,
    upper_ci = upper_bound
  )
  
  model_results <- data.table(true_tate,true_lte,ncs_tate,ncs_lte,results)
  return(model_results)
  
  
}

s_single_rep <- function(iter,  ncluster) {
  
  train_data <- s_generate(iter,  ncluster)
  
  model_results <- s_model(train_data)
  
  return(model_results)
}

s_replicate <- function(iter, ncluster) {
  model_results = s_single_rep(iter, ncluster)
  return(data.table(iter=iter, ncluster=ncluster, model_results))
}


scenarios = expand.grid(ncluster=10)
i=1
nsim=150
ncluster= scenarios[i,"ncluster"]

start_time <- Sys.time()

res <- lapply(
  X = 1:nsim, 
  FUN = function(x) s_replicate(iter=x, ncluster=ncluster))
res <- rbindlist(res)
save(res, file = "output/R1_q4/inc_dec_time_vary_trt_ncs.rda")

