
load("C:/My PC (DESKTOP-52UHO48)/Desktop/Bayesian Spline/Bayesian stepped wedge/Bayesian stepped wedge/Bayesian-stepped-wedge/Final simulation code/Code for paper revision/output/R2_q5/immediate_trt_scaleyk.rda")
head(res)

performance <- function(true_value, est) {
  bias = (mean(est) - true_value)#/true_value *100
  rmse = sqrt(mean((true_value - est)^2))
  return(list(bias=bias, rmse=rmse))
}

performance = performance(unique(res$coefA), res$median)
# $bias
# [1] -0.0210033
# 
# $rmse
# [1] 0.1499974
# > mean(res$covered_bayes)
# [1] 0.94
####div <1%
# > res <- res%>%filter(div<=100)
# > dim(res)
# [1] 148  15
# > performance = performance(unique(res$coefA), res$median)
# > performance
# $bias
# [1] -0.01849814
# 
# $rmse
# [1] 0.1485978
# 
# > mean(res$covered_bayes)
# [1] 0.9459459

load("C:/My PC (DESKTOP-52UHO48)/Desktop/Bayesian Spline/Bayesian stepped wedge/Bayesian stepped wedge/Bayesian-stepped-wedge/Final simulation code/Code for paper revision/output/R2_q5/cluster_specific_trt_scaley.rda")

colnames(res) <- c( 
  "iter","ncluster","true_tate","true_lte", "variable_1", "mean_tate","median_tate",            
  "sd_tate","2.5%_tate","97.5%_tate","rhat_tate","ess_bulk_tate","ess_tail_tate","variable", "mean_lte",              
  "median_lte","sd_lte","2.5%_lte","97.5%_lte","rhat_lte","ess_bulk_lte","ess_tail_lte","covered_bayes_tate",
  "covered_bayes_lte","div","num_max_tree_depth")

dim(res)
bayes_tate_pf = performance(unique(res$true_tate), res$median_tate)
bayes_lte_pf = performance(unique(res$true_lte), res$median_lte)
bayes_tate_pf 
bayes_lte_pf

library(dplyr)

load("C:/My PC (DESKTOP-52UHO48)/Desktop/Bayesian Spline/Bayesian stepped wedge/Bayesian stepped wedge/Bayesian-stepped-wedge/Final simulation code/Code for paper revision/output/R2_q5/cluster_specific_trt_scaleykt.rda")
colnames(res) <- c( 
  "iter","ncluster","true_tate","true_lte", "variable_1", "mean_tate","median_tate",            
  "sd_tate","2.5%_tate","97.5%_tate","rhat_tate","ess_bulk_tate","ess_tail_tate","variable", "mean_lte",              
  "median_lte","sd_lte","2.5%_lte","97.5%_lte","rhat_lte","ess_bulk_lte","ess_tail_lte","covered_bayes_tate",
  "covered_bayes_lte","div","num_max_tree_depth")
bayes_tate_pf = performance(unique(res$true_tate), res$median_tate)
bayes_lte_pf = performance(unique(res$true_lte), res$median_lte)
bayes_tate_pf 
bayes_lte_pf

