library(dplyr)

load("C:/My PC (DESKTOP-52UHO48)/Desktop/Bayesian Spline/Bayesian stepped wedge/Bayesian stepped wedge/Bayesian-stepped-wedge/Final simulation code/Code for paper revision/R/R2_q3/inc_dec_time_varying_pen_spline.rda")
head(res)
colnames(res) <- c( "iter","ncluster","true_tate","true_lte", "variable_1", "mean_tate","median_tate",            
                    "sd_tate","2.5%_tate","97.5%_tate","rhat_tate","ess_bulk_tate","ess_tail_tate","variable", "mean_lte",              
                    "median_lte","sd_lte","2.5%_lte","97.5%_lte","rhat_lte","ess_bulk_lte","ess_tail_lte","covered_bayes_tate",
                    "covered_bayes_lte","div","num_max_tree_depth","elpd_loo",
                    "elpd_se", "p_loo","p_loo_se","looic", "looic_se")




performance <- function(true_value, est) {
  bias = abs((mean(est) - true_value))#/true_value #*100
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

tate_data_2 <- data.frame(
  Measure ="TATE",
  Scenario ="2",
  Model = c("Bayesian model (6)", "Monotone effect curve model"),
  Bias = abs(c(-0.3083283, -7.55285)),
  RMSE = c(0.9335343, 10.63354),
  LOO = c(6763, 6842)
)

# Define the data for LTE scenario
lte_data_2 <- data.frame(
  Measure ="LTE",
  Scenario ="2",
  Model = c("Bayesian model (6)", "Monotone effect curve model"),
  Bias = abs(c(-0.2787317, -9.749612)),
  RMSE = c(1.698069, 14.52556),
  LOO = c(6763, 6842)
)