####time vary

##increase
#monotone increase
load("C:/My PC (DESKTOP-52UHO48)/Desktop/Bayesian Spline/Bayesian stepped wedge/Bayesian stepped wedge/Bayesian-stepped-wedge/Final simulation code/Data/res_increase_t_vary_mono.rda")
colnames(res) <- c( "iter","coefA","ncluster","true_tate","true_lte", "variable_1", "mean_tate","median_tate",            
                    "sd_tate","2.5%_tate","97.5%_tate","rhat_tate","ess_bulk_tate","ess_tail_tate","variable", "mean_lte",              
                    "median_lte","sd_lte","2.5%_lte","97.5%_lte","rhat_lte","ess_bulk_lte","ess_tail_lte","covered_bayes_monot_tate",
                    "covered_bayes_monot_lte","div_monot","num_max_tree_depth_monot","elpd_loo",
                    "elpd_se", "p_loo","p_loo_se","looic", "looic_se")
##m2 
load("C:/My PC (DESKTOP-52UHO48)/Desktop/Bayesian Spline/Bayesian stepped wedge/Bayesian stepped wedge/Bayesian-stepped-wedge/Final simulation code/Data/res_increase_t_vary_our_m2.rda")
colnames(res) <- c( "iter","coefA","ncluster","true_tate","true_lte", "variable_1", "mean_tate","median_tate",            
                    "sd_tate","2.5%_tate","97.5%_tate","rhat_tate","ess_bulk_tate","ess_tail_tate","variable", "mean_lte",              
                    "median_lte","sd_lte","2.5%_lte","97.5%_lte","rhat_lte","ess_bulk_lte","ess_tail_lte","covered_bayes_tate",
                    "covered_bayes_lte","gam_tate","gam_lte","div","num_max_tree_depth","elpd_loo",
                    "elpd_se", "p_loo","p_loo_se","looic", "looic_se")

#dec inc
load("C:/My PC (DESKTOP-52UHO48)/Desktop/Bayesian Spline/Bayesian stepped wedge/Bayesian stepped wedge/Bayesian-stepped-wedge/Final simulation code/Data/res_dec_increase_t_vary_mono.rda")
colnames(res) <- c( "iter","coefA","ncluster","true_tate","true_lte", "variable_1", "mean_tate","median_tate",            
                    "sd_tate","2.5%_tate","97.5%_tate","rhat_tate","ess_bulk_tate","ess_tail_tate","variable", "mean_lte",              
                    "median_lte","sd_lte","2.5%_lte","97.5%_lte","rhat_lte","ess_bulk_lte","ess_tail_lte","covered_bayes_monot_tate",
                    "covered_bayes_monot_lte","div_monot","num_max_tree_depth_monot","elpd_loo",
                    "elpd_se", "p_loo","p_loo_se","looic", "looic_se")


##m2 
load("C:/My PC (DESKTOP-52UHO48)/Desktop/Bayesian Spline/Bayesian stepped wedge/Bayesian stepped wedge/Bayesian-stepped-wedge/Final simulation code/Data/res_inc_decrease_t_vary_m2.rda")
colnames(res) <- c( "iter","ncluster","true_tate","true_lte", "variable_1", "mean_tate","median_tate",            
                    "sd_tate","2.5%_tate","97.5%_tate","rhat_tate","ess_bulk_tate","ess_tail_tate","variable", "mean_lte",              
                    "median_lte","sd_lte","2.5%_lte","97.5%_lte","rhat_lte","ess_bulk_lte","ess_tail_lte","covered_bayes_tate",
                    "covered_bayes_lte","gam_tate","gam_lte","div","num_max_tree_depth","elpd_loo",
                    "elpd_se", "p_loo","p_loo_se","looic", "looic_se")