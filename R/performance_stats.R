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

res <- res %>%filter(div <= 60)
dim(res)
# Example usage (assuming `results_agg` is properly defined before these function calls):
truth=0
truth=0.25
truth=0.75
truth=0.5
truth=1
freq_sat = performance(truth, res$est_sat, res$se_sat)
round(unlist(freq_sat),3)

freq_p = performance(truth, res$est_gam_freq, res$se_gam_freq)
round(unlist(freq_p),3)
bayes_p_med = bayes_performance(truth, res$est_med_bayes, res$lowci_bayes, res$upci_bayes,res$est_sd_bayes)
round(unlist(bayes_p_med),3)
bayes_p_mean = bayes_performance(truth, res$est_mean_bayes, res$lowci_bayes, res$upci_bayes,res$est_sd_bayes)
round(unlist(bayes_p_mean),3)

bayes_p_sd10 = bayes_performance(truth, res$est_med_bayes10, res$lowci_bayes10, res$upci_bayes10,res$est_sd_bayes10)
round(unlist(bayes_p_sd10),3)

freq_rdn = performance(truth, res$est_gam_rdn_freq, res$se_gam_rdn_freq)
round(unlist(freq_rdn),3)

freq_sat = performance(truth, res$est_sat, res$se_sat)
round(unlist(freq_sat),3)

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
  freq_rdn$bias,
  bayes_p_med$bias,
  bayes_p_sd10$bias
)
round(bias_compare,3)

# Extract Type I error rates
type1_errors <- c(
  freq_sat$type1error,
  freq_notime$type1error,
  freq_lntime$type1error,
  freq_p$type1error,
  freq_rdn$type1error,
  bayes_p_med$type1error,
  bayes_p_sd10$type1error
)
round(type1_errors,3)

coverage_true <- c(
  freq_sat$coverage,
  freq_notime$coverage,
  freq_lntime$coverage,
  freq_p$coverage,
  freq_rdn$coverage,
  bayes_p_med$coverage,
  bayes_p_sd10$coverage
)
round(coverage_true,3)

rmse_value <- c(
  freq_sat$rmse,
  freq_notime$rmse,
  freq_lntime$rmse,
  freq_p$rmse,
  freq_rdn$rmse,
  bayes_p_med$rmse,
  bayes_p_sd10$rmse
)
round(rmse_value,3)
#power


res <- res %>%filter(div <= 60)
dim(res)
freq_sat = power_func(res$sat_lowci)
round(unlist(freq_sat),3)

freq_notime  = power_func(res$notime_lowci)
round(unlist(freq_notime ),3)

freq_lntime = power_func(res$lntime_lowci)
round(unlist(freq_lntime),3)

freq_p = power_func(res$lowci_freq)
round(unlist(freq_p),3)

freq_rdn = power_func(res$lowci_freq_rdn)
round(unlist(freq_rdn),3)

bayes_p_med = power_func(res$lowci_bayes)
round(unlist(bayes_p_med),3)

bayes_p_sd10 = power_func(res$lowci_bayes10)
round(unlist(bayes_p_sd10),3)

power <- c(round(unlist(freq_sat),3),
           round(unlist(freq_notime ),3),
           round(unlist(freq_lntime),3),
           round(unlist(freq_p),3),
           round(unlist(freq_rdn),3),
           round(unlist(bayes_p_med),3),
           round(unlist(bayes_p_sd10),3)
           )
power
#power


