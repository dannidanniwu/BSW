load("C:/My PC (DESKTOP-52UHO48)/Desktop/Bayesian Spline/Bayesian stepped wedge/Bayesian stepped wedge/Bayesian-stepped-wedge/Final simulation code/Code for paper revision/output/R1_q4/inc_time_vary_trt_ncs.rda")

summary(res)

library(data.table)


res$covered <- with(res, lower_ci < true_eff & true_eff < upper_ci)


# Calculate the mean of 'covered' by 'variable'
coverage_rate <- res[, .(mean_coverage = round(mean(covered)*100,2)), by = time]

# View the coverage rate for each variable
print(coverage_rate)

library(dplyr)

# Assuming your data frame is called df
df <- res%>%select(iter, true_tate, true_lte,ncs_tate,ncs_lte)
df_unique <- df %>% distinct(iter, .keep_all = TRUE)

performance <- function(true_value, est) {
  bias = abs((mean(est) - true_value))#/true_value *100
  rmse = sqrt(mean((true_value - est)^2))
  return(list(bias=bias, rmse=rmse))
}

performance_tate = performance(unique(df_unique$true_tate), df_unique$ncs_tate)
performance_lte = performance(unique(df_unique$true_lte), df_unique$ncs_lte)

performance_tate 

performance_lte


load("C:/My PC (DESKTOP-52UHO48)/Desktop/Bayesian Spline/Bayesian stepped wedge/Bayesian stepped wedge/Bayesian-stepped-wedge/Final simulation code/Code for paper revision/output/R1_q4/inc_dec_time_vary_trt_ncs.rda")
res$covered <- with(res, lower_ci < true_eff & true_eff < upper_ci)


# Calculate the mean of 'covered' by 'variable'
coverage_rate <- res[, .(mean_coverage = round(mean(covered)*100,2)), by = time]

# View the coverage rate for each variable
print(coverage_rate)

# Assuming your data frame is called df
df <- res%>%select(iter, true_tate, t_effect,ncs_tate,ncs_lte)
df_unique <- df %>% distinct(iter, .keep_all = TRUE)

performance_tate = performance(unique(df_unique$true_tate), df_unique$ncs_tate)
performance_lte = performance(unique(df_unique$t_effect), df_unique$ncs_lte)

performance_tate 

performance_lte

load("C:/My PC (DESKTOP-52UHO48)/Desktop/Bayesian Spline/Bayesian stepped wedge/Bayesian stepped wedge/Bayesian-stepped-wedge/Final simulation code/Code for paper revision/output/R1_q4/inc_dec_time_cluster_spc_vary_trt_ncs.rda")
res$covered <- with(res, lower_ci < true_eff & true_eff < upper_ci)


# Calculate the mean of 'covered' by 'variable'
coverage_rate <- res[, .(mean_coverage = round(mean(covered)*100,2)), by = time]

# View the coverage rate for each variable
print(coverage_rate)

df <- res%>%select(iter, true_tate, t_effect_overall,ncs_tate,ncs_lte)
df_unique <- df %>% distinct(iter, .keep_all = TRUE)

performance_tate = performance(unique(df_unique$true_tate), df_unique$ncs_tate)
performance_lte = performance(unique(df_unique$t_effect), df_unique$ncs_lte)

performance_tate 

performance_lte
