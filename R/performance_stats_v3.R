library(dplyr)
coverage <- function(x, a, b) {
  range <- a + c(-1,1) * 1.96 * b
  (range[1] < x & x < range[2])
}

performance <- function(true_value, est, se) {
  
  bias = mean(est) - true_value
  rmse = sqrt(mean((true_value - est)^2))
  true_value.se = sd(est)
  est.se = mean(se)
  coverage = mean(mapply("coverage", true_value, est, se)) * 100
  
  return(list(bias = bias, rmse = rmse, true_value.se = true_value.se, 
              est.se = est.se, coverage = coverage, type1error = 100-coverage))
}

bayes_coverage <- function(true_value, lowci, upci){
  (lowci < true_value & true_value < upci)
}

bayes_performance <- function(true_value, est, lowci, upci, se){
  bias = mean(est) - true_value
  rmse = sqrt(mean((true_value - est)^2))
  est.se = mean(se)
  bay_coverage = mean(mapply(bayes_coverage, true_value, lowci, upci)) * 100
  
  return(list(bias = bias, rmse = rmse, est.se=est.se,coverage = bay_coverage, type1error = 100-bay_coverage))
}

power_func <- function(lowci){
  mean(lowci>0)*100
}

res <- res %>%filter(div<=80)

dim(res)
# Example usage (assuming `results_agg` is properly defined before these function calls):
truth=0
truth=0.25
truth=0.75
truth=0.5
truth=1
truth=5
truth=4
freq_sat = performance(truth, res$est_sat, res$se_sat)
round(unlist(freq_sat),3)

freq_p = performance(truth, res$est_gam_freq, res$se_gam_freq)
round(unlist(freq_p),3)
bayes_p_med = bayes_performance(truth, res$est_med_bayes, res$lowci_bayes, res$upci_bayes,res$est_sd_bayes)
round(unlist(bayes_p_med),3)
bayes_p_mean = bayes_performance(truth, res$est_mean_bayes, res$lowci_bayes, res$upci_bayes,res$est_sd_bayes)
round(unlist(bayes_p_mean),3)



freq_notime  = performance(truth, res$est_notime , res$se_notime )
round(unlist(freq_notime ),truth)

freq_lntime = performance(truth, res$est_lntime, res$se_lntime)
round(unlist(freq_lntime),3)


#bias
bias_compare <- c(
  freq_sat$bias,
  freq_notime$bias,
  freq_lntime$bias,
  freq_p$bias,
  bayes_p_med$bias

)
round(bias_compare,3)

# Extract Type I error rates
type1_errors <- c(
  freq_sat$type1error,
  freq_notime$type1error,
  freq_lntime$type1error,
  freq_p$type1error,
  bayes_p_med$type1error

)
round(type1_errors,3)

coverage_true <- c(
  freq_sat$coverage,
  freq_notime$coverage,
  freq_lntime$coverage,
  freq_p$coverage,
  bayes_p_med$coverage
)


round(coverage_true,3)

rmse_value <- c(
  freq_sat$rmse,
  freq_notime$rmse,
  freq_lntime$rmse,
  freq_p$rmse,
  bayes_p_med$rmse
)


round(rmse_value,3)
#power

freq_sat = power_func(res$sat_lowci)
round(unlist(freq_sat),3)

freq_notime  = power_func(res$notime_lowci)
round(unlist(freq_notime ),3)

freq_lntime = power_func(res$lntime_lowci)
round(unlist(freq_lntime),3)

freq_p = power_func(res$lowci_freq)
round(unlist(freq_p),3)


bayes_p_med = power_func(res$lowci_bayes)
round(unlist(bayes_p_med),3)

bayes_p_sd2_5 = power_func(res$lowci_bayes2_5)
round(unlist(bayes_p_sd2_5),3)
bayes_p_sd1 = power_func(res$lowci_bayes1)
round(unlist(bayes_p_sd1),3)

power_res <- c(round(unlist(freq_sat),3),
           round(unlist(freq_notime ),3),
           round(unlist(freq_lntime),3),
           round(unlist(freq_p),3),
           
           round(unlist(bayes_p_med),3)
           )
power_res
#power

bayes_cri_sd2_5 = res$upci_bayes2_5-res$lowci_bayes2_5
summary(bayes_cri_sd2_5)
bayes_cri_sd1 = res$upci_bayes1 -res$lowci_bayes1
summary(bayes_cri_sd1)
bayes_cri_sd5 = res$upci_bayes -res$lowci_bayes
summary(bayes_cri_sd5)
freq_ci =res$upci_freq -res$lowci_freq
summary(freq_ci)


