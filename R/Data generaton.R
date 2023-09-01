#Data generation
#In this data generation process, the time effect will not be explicitly smooth, but the underlying covariance structure used to generate the period effects will induce some level of smoothness. 
#references:https://www.rdatagen.net/post/2022-12-13-modeling-the-secular-trend-in-a-stepped-wedge-design/
#references:https://www.rdatagen.net/post/2022-11-01-modeling-secular-trend-in-crt-using-gam/
library(simstudy)
library(ggplot2)
library(data.table)
library(mgcv)
library(lme4)
library(splines)

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
  
  return(dd) #  generated_data is a data.table
}

dp <- dd[, .(avg = mean(y)), keyby = .(A, site, k)]

ggplot(data = dp, aes(x = k, y = avg)) +
  geom_line(aes(group = site, color = factor(A))) +
  scale_color_manual(values = c("#a53e3e", "#3ea5a5"),
                     guide = guide_legend(reverse = TRUE),
                     labels = c("no treatment", "treatment")) +
  ylab("average Y") +
  xlab("period (k)") +
  theme(panel.grid = element_blank(),
        legend.title = element_blank()) 