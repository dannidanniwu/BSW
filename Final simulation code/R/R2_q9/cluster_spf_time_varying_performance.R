library(dplyr)
load("C:/Users/Danni/OneDrive - NYU Langone Health/Bayesian stepped wedge/Bayesian-stepped-wedge/Data/res_dec_increase_cluster_specific_t_vary_m3.rda")
colnames(res) <- c("iter","ncluster","true_tate","true_lte", "variable_1", "mean_tate","median_tate",            
                   "sd_tate","2.5%_tate","97.5%_tate","rhat_tate","ess_bulk_tate","ess_tail_tate","variable", "mean_lte",              
                   "median_lte","sd_lte","2.5%_lte","97.5%_lte","rhat_lte","ess_bulk_lte","ess_tail_lte","covered_bayes_tate",
                   "covered_bayes_lte","div","num_max_tree_depth","elpd_loo",
                   "elpd_se", "p_loo","p_loo_se","looic", "looic_se")

dim(res)


performance <- function(true_value, est) {
  bias = (mean(est) - true_value)#/true_value *100
  rmse = sqrt(mean((true_value - est)^2))
  return(list(bias=bias, rmse=rmse))
}

gam_tate_pf = performance(unique(res$true_tate), res$gam_tate)
gam_lte_pf = performance(unique(res$true_lte), res$gam_lte)
gam_tate_pf 
gam_lte_pf

bayes_tate_pf = performance(unique(res$true_tate), res$median_tate)
bayes_lte_pf = performance(unique(res$true_lte), res$median_lte)
bayes_tate_pf 
bayes_lte_pf

covered_bayes_tate <- mean(res$covered_bayes_tate)
covered_bayes_lte <- mean(res$covered_bayes_lte)
covered_bayes_tate
covered_bayes_lte

summary(res$rhat_lte)
summary(res$looic)


load("C:/Users/Danni/OneDrive - NYU Langone Health/Bayesian stepped wedge/Bayesian-stepped-wedge/Data/res_cluster_spf_t_vary_our_mono.rda")

colnames(res) <- c("iter","ncluster","true_tate","true_lte", "variable_1", "mean_tate","median_tate",            
                    "sd_tate","2.5%_tate","97.5%_tate","rhat_tate","ess_bulk_tate","ess_tail_tate","variable", "mean_lte",              
                    "median_lte","sd_lte","2.5%_lte","97.5%_lte","rhat_lte","ess_bulk_lte","ess_tail_lte","covered_bayes_tate",
                    "covered_bayes_lte","div","num_max_tree_depth","elpd_loo",
                   "elpd_se", "p_loo","p_loo_se","looic", "looic_se")
#res <- res %>%filter(div_monot <=80 & num_max_tree_depth_monot <=80)
dim(res)

monot_tate_pf = performance(unique(res$true_tate), res$median_tate)
monot_lte_pf = performance(unique(res$true_lte), res$median_lte)
monot_tate_pf 
monot_lte_pf


mean(res$covered_bayes_monot_tate)
mean(res$covered_bayes_monot_lte)
summary(res$rhat_lte)
summary(res$looic)

load("/Users/danni/Library/CloudStorage/OneDrive-NYULangoneHealth/Bayesian stepped wedge/Bayesian-stepped-wedge/Final simulation code/Data/res_dec_increase_cluster_specific_t_vary_m2.rda")
colnames(res) <- c("iter","ncluster","true_tate","true_lte", "variable_1", "mean_tate","median_tate",            
                   "sd_tate","2.5%_tate","97.5%_tate","rhat_tate","ess_bulk_tate","ess_tail_tate","variable", "mean_lte",              
                   "median_lte","sd_lte","2.5%_lte","97.5%_lte","rhat_lte","ess_bulk_lte","ess_tail_lte","covered_bayes_tate",
                   "covered_bayes_lte","gam_tate","gam_lte","div","num_max_tree_depth","elpd_loo",
                   "elpd_se", "p_loo","p_loo_se","looic", "looic_se")

dim(res)


performance <- function(true_value, est) {
  bias = (mean(est) - true_value)#/true_value *100
  rmse = sqrt(mean((true_value - est)^2))
  return(list(bias=bias, rmse=rmse))
}


bayes_tate_pf = performance(unique(res$true_tate), res$median_tate)
bayes_lte_pf = performance(unique(res$true_lte), res$median_lte)
bayes_tate_pf 
bayes_lte_pf

covered_bayes_tate <- mean(res$covered_bayes_tate)
covered_bayes_lte <- mean(res$covered_bayes_lte)
covered_bayes_tate
covered_bayes_lte

summary(res$rhat_lte)
summary(res$looic)