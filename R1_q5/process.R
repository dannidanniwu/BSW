library(dplyr)
#model 3
############update date genration m3
load("C:/My PC (DESKTOP-52UHO48)/Desktop/Bayesian Spline/Bayesian stepped wedge/Bayesian stepped wedge/Bayesian-stepped-wedge/Final simulation code/Code for paper revision/output/R1_q5/time_varying_gen_m3fit_update.rda")

colnames(res) <- c( "iter","ncluster","true_tate","true_lte", "variable_1", "mean_tate","median_tate",            
                    "sd_tate","2.5%_tate","97.5%_tate","rhat_tate","ess_bulk_tate","ess_tail_tate","variable", "mean_lte",              
                    "median_lte","sd_lte","2.5%_lte","97.5%_lte","rhat_lte","ess_bulk_lte","ess_tail_lte","covered_bayes_tate",
                    "covered_bayes_lte","div","num_max_tree_depth","elpd_loo",
                    "elpd_se", "p_loo","p_loo_se","looic", "looic_se")




performance <- function(true_value, est) {
  bias = (mean(est) - true_value)#/true_value #*100
  rmse = sqrt(mean((true_value - est)^2))
  return(list(bias=bias, rmse=rmse))
}

true_tate = res$true_tate
true_lte = res$true_lte
bayes_tate_pf = performance(unique(res$true_tate), res$median_tate)
bayes_lte_pf = performance(unique(res$true_lte), res$median_lte)
bayes_tate_pf 
bayes_lte_pf


summary(res$rhat_lte)
summary(res$looic)
#model 2

load("C:/My PC (DESKTOP-52UHO48)/Desktop/Bayesian Spline/Bayesian stepped wedge/Bayesian stepped wedge/Bayesian-stepped-wedge/Final simulation code/Data/res_inc_decrease_t_vary_m2.rda")
colnames(res) <- c( "iter","ncluster","true_tate","true_lte", "variable_1", "mean_tate","median_tate",            
                    "sd_tate","2.5%_tate","97.5%_tate","rhat_tate","ess_bulk_tate","ess_tail_tate","variable", "mean_lte",              
                    "median_lte","sd_lte","2.5%_lte","97.5%_lte","rhat_lte","ess_bulk_lte","ess_tail_lte","covered_bayes_tate",
                    "covered_bayes_lte","gam_tate","gam_lte","div","num_max_tree_depth","elpd_loo",
                    "elpd_se", "p_loo","p_loo_se","looic", "looic_se")

bayes_tate_pf = performance(unique(true_tate), res$median_tate)
bayes_lte_pf = performance(unique(true_lte), res$median_lte)
bayes_tate_pf 
bayes_lte_pf


summary(res$rhat_lte)
summary(res$looic)


load("C:/My PC (DESKTOP-52UHO48)/Desktop/Bayesian Spline/Bayesian stepped wedge/Bayesian stepped wedge/Bayesian-stepped-wedge/Final simulation code/Data/res_dec_increase_t_vary_mono.rda")
colnames(res) <- c( "iter","ncluster","true_tate","true_lte", "variable_1", "mean_tate","median_tate",            
                    "sd_tate","2.5%_tate","97.5%_tate","rhat_tate","ess_bulk_tate","ess_tail_tate","variable", "mean_lte",              
                    "median_lte","sd_lte","2.5%_lte","97.5%_lte","rhat_lte","ess_bulk_lte","ess_tail_lte","covered_bayes_monot_tate",
                    "covered_bayes_monot_lte","div_monot","num_max_tree_depth_monot","elpd_loo",
                    "elpd_se", "p_loo","p_loo_se","looic", "looic_se")


monot_tate_pf = performance(unique(true_tate), res$median_tate)
monot_lte_pf = performance(unique(true_lte), res$median_lte)
monot_tate_pf 
monot_lte_pf


summary(res$rhat_lte)
summary(res$looic)
