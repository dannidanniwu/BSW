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
              est.se = est.se, coverage = coverage))
}

bayes_coverage <- function(true_value, lowci, upci){
  (lowci < true_value & true_value < upci)
}

bayes_performance <- function(true_value, est, lowci, upci, se){
  bias = mean(est) - true_value
  rmse = sqrt(mean((true_value - est)^2))
  est.se = mean(se)
  bay_coverage = mean(mapply(bayes_coverage, true_value, lowci, upci)) * 100
  
  return(list(bias = bias, rmse = rmse, est.se=est.se,coverage = bay_coverage))
}

# Example usage (assuming `results_agg` is properly defined before these function calls):
freq_p = performance(5, results_agg$est.freq, results_agg$se.freq)
round(unlist(freq_p),3)
bayes_p_med = bayes_performance(5, results_agg$est.med.bayes, results_agg$lowci.bayes, results_agg$upci.bayes,results_agg$est.sd.bayes)
round(unlist(bayes_p_med),3)
bayes_p_mean = bayes_performance(5, results_agg$est.mean.bayes, results_agg$lowci.bayes, results_agg$upci.bayes,results_agg$est.sd.bayes)
round(unlist(bayes_p_mean),3)

