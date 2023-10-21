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

power <- function(lowci){
  mean(lowci>0)*100
}

res <- res %>%filter(div <= 100)
dim(res)
# Example usage (assuming `results_agg` is properly defined before these function calls):

freq_p = performance(2, res$est_gam_freq, res$se_gam_freq)
round(unlist(freq_p),3)
bayes_p_med = bayes_performance(2, res$est_med_bayes, res$lowci_bayes, res$upci_bayes,res$est_sd_bayes)
round(unlist(bayes_p_med),3)
bayes_p_mean = bayes_performance(2, res$est_mean_bayes, res$lowci_bayes, res$upci_bayes,res$est_sd_bayes)
round(unlist(bayes_p_mean),3)

freq_rdn = performance(2, res$est_gam_rdn_freq, res$se_gam_rdn_freq)
round(unlist(freq_rdn),3)

freq_sat = performance(2, res$est_sat, res$se_sat)
round(unlist(freq_sat),3)

freq_notime  = performance(2, res$est_notime , res$se_notime )
round(unlist(freq_notime ),3)

freq_lntime = performance(2, res$est_lntime, res$se_lntime)
round(unlist(freq_lntime),3)


#type 1 error rate



freq_sat = power(res$sat_lowci)
round(unlist(freq_sat),3)

freq_notime  = power(res$notime_lowci)
round(unlist(freq_notime ),3)

freq_lntime = power(res$lntime_lowci)
round(unlist(freq_lntime),3)

freq_p = power(res$lowci_freq)
round(unlist(freq_p),3)

freq_rdn = power(res$lowci_freq_rdn)
round(unlist(freq_rdn),3)

bayes_p_med = power(res$lowci_bayes)
round(unlist(bayes_p_med),3)

#power


