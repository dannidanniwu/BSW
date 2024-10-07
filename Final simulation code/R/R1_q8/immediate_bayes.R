library(simstudy)
library(data.table)
library(mgcv)
library(lme4)
library(splines)
library(rstan)
library(brms)
library(parallel)
library(dplyr)
library(priorsense)

# Compile the Stan model
#mod <- cmdstan_model("./stepped_wedge_random_walk_prior_test.stan")
#mod <- cmdstan_model("./stepped_wedge_random_walk_fixed_trt.stan")
mod <- stan_model("./immediate_trt_priorsense.stan")


s_generate <- function(iter, coefA, ncluster) {
  
  set.seed(iter)
  def <- defData(varname = "a", formula = 0, variance = 0.25)
  def <- defData(def, varname = "mu_b", formula = 0, dist = "nonrandom")
  def <- defData(def, varname = "s2_b", formula = 0.36, dist = "nonrandom")
  #A: trt for each cluster and time period
  #b: site-specific time period effect
  defOut <- defDataAdd(varname = "y", formula = " a + b - 0.05 * k^2 + ..coefA * A", variance = 1)
  
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


i=2
nsim=150
coefA = 5
ncluster= 10
iter=1
train_data <- s_generate(iter, coefA, ncluster)
  
  # Define the knots based on the training data
  knots_train <- quantile(train_data$k, probs=seq(0, 1, length=6)[-c(1, 6)])
  
  B_train <- predict(splines::bs(train_data$k, degree=3, knots=quantile(train_data$k, probs=seq(0, 1, length=6)[-c(1, 6)])))
  colnames(B_train) <- paste0("Bspline_", 1:ncol(B_train))
  
  
  stan_data <- list(num_data = nrow(train_data),
                    num_basis = ncol(B_train),
                    B = t(B_train),
                    A = train_data$A,
                    y = train_data$y,
                    num_clusters = length(unique(train_data$site)),
                    cluster = as.numeric(train_data$site)
                    
  )
  #mod <- cmdstan_model("./stepped_wedge_random_walk.stan")
  # Fit the Bayesian model
  fit <- sampling(mod, 
                  data = stan_data,        # The data list for Stan
                  chains = 4,              # Number of chains
                  warmup = 500,            # Warmup iterations (burn-in)
                  iter = 2500,             # Total iterations (warmup + sampling)
                  cores = 4,               # Number of cores for parallel chains
                  refresh = 0,             # Refresh frequency to control message output
                  verbose = TRUE)          # Show messages 
  
  res <- powerscale_sensitivity(fit)


  