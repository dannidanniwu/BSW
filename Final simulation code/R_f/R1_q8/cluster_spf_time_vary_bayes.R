library(simstudy)
library(data.table)
library(mgcv)
library(lme4)
library(splines)
library(rstan)
library(brms)
library(dplyr)
library(loo)
library(priorsense)

mod <- stan_model("./cluster_specific_time_varying_trt_priorsense.stan")


s_generate <- function(iter,  ncluster) {
  
  set.seed(iter)
  def <- defData(varname = "a", formula = 0, variance = 0.25)
  def <- defData(def, varname = "rdn", formula = 0, variance = 0.04)
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
  
  defteffect_rdn <- defDataAdd(varname = "t_effect", 
                               formula = "exp(rdn)*t_effect_overall", 
                               dist= "nonrandom")
  
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
  dd <- addCondition(defteffect, dd, newvar = "t_effect_overall")
  dd <- addColumns(defteffect_rdn, dd)
  dd <- addColumns(defOut, dd)
  dd[, normk := (k - min(k))/(max(k) - min(k))]#scale time period into range 0-1
  
  dd[] #  generated_data is a data.table
  
}


iter=1#nsim=150
ncluster= 10
train_data <- s_generate(iter, ncluster)

t_ex_mx=max(train_data$t)
#true value of TATE
t_values <- 1:t_ex_mx
true_tate = mean(unique(train_data[t%in%t_values,t_effect_overall]))
#true value of LTE
true_lte = unique(train_data[t==t_ex_mx ,"t_effect_overall"])
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
                  cluster = as.numeric(train_data$site)
                  
)
# Fit the Bayesian model
fit <- sampling(mod, 
                data = stan_data,        # The data list for Stan
                chains = 4,              # Number of chains
                warmup = 500,            # Warmup iterations (burn-in)
                iter = 2500,             # Total iterations (warmup + sampling)
                cores = 4,               # Number of cores for parallel chains
                refresh = 0,             # Refresh frequency to control message output
                verbose = TRUE)          # Show messages 
setwd("C:/My PC (DESKTOP-52UHO48)/Desktop/Bayesian Spline/Bayesian stepped wedge/Bayesian stepped wedge/Bayesian-stepped-wedge/Final simulation code/Code for paper revision/Stan")
res <- powerscale_sensitivity(fit)
save(res,file="./cluster_spf_time_vary_bayes_priorsense_res.rda")