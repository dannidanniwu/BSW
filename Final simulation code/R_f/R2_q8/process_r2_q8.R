load("C:/My PC (DESKTOP-52UHO48)/Desktop/Bayesian Spline/Bayesian stepped wedge/Bayesian stepped wedge/Bayesian-stepped-wedge/Final simulation code/Code for paper revision/output/R2_q8/inc_dec_trt_m2_pointwise_ci.rda")
head(res)
summary(res)

library(data.table)

res$covered <- with(res, `2.5%` < true_eff & true_eff < `97.5%`)

# Calculate the mean of 'covered' by 'variable'
coverage_rate <- res[, .(mean_coverage = round(mean(covered)*100,2)), by = variable]

# View the coverage rate for each variable
print(coverage_rate)
#############ncs model
load("C:/My PC (DESKTOP-52UHO48)/Desktop/Bayesian Spline/Bayesian stepped wedge/Bayesian stepped wedge/Bayesian-stepped-wedge/Final simulation code/Code for paper revision/output/R2_q8/inc_dec_cluster_specific_trt_m3_pointwise_ci.rda")
head(res)
summary(res)

res$covered <- with(res, `2.5%` < true_eff & true_eff < `97.5%`)

# Calculate the mean of 'covered' by 'variable'
coverage_rate <- res[, .(mean_coverage = round(mean(covered)*100,2)), by = variable]
print(coverage_rate)

load("C:/My PC (DESKTOP-52UHO48)/Desktop/Bayesian Spline/Bayesian stepped wedge/Bayesian stepped wedge/Bayesian-stepped-wedge/Final simulation code/Code for paper revision/output/R1_q4/inc_dec_time_cluster_spc_vary_trt_ncs.rda")
res$covered <- with(res, lower_ci < true_eff & true_eff < upper_ci)


# Calculate the mean of 'covered' by 'variable'
coverage_rate_ncs <- res[, .(mean_coverage = round(mean(covered)*100,2)), by = time]

# View the coverage rate for each variable
print(coverage_rate_ncs)

all <- cbind(coverage_rate,coverage_rate_ncs)

write.csv(all, file ='./m3_ncs_compare.csv')
