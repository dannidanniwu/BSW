library(simstudy)
library(data.table)
library(mgcv)
library(lme4)
library(splines)
library(brms)
library(cmdstanr)
library(slurmR)
library(dplyr)
#references:https://www.rdatagen.net/post/2022-12-13-modeling-the-secular-trend-in-a-stepped-wedge-design/
#references:https://www.rdatagen.net/post/2022-11-01-modeling-secular-trend-in-crt-using-gam/

# Compile the Stan model
#mod <- cmdstan_model("./stepped_wedge_time_spline.stan")
#mod <- cmdstan_model("./stepped_wedge_time_spline_penalized.stan")
mod <- cmdstan_model("./stepped_wedge_test_hier.stan")

s_define <- function() {
  #cluster-specific intercept
  def <- defData(varname = "a", formula = 0, variance = 1)
  def2 <- defDataAdd(varname = "b", formula = "(k - 0.5)^2", variance =0.4)
  #A: trt for each cluster and time period
  defOut <- defDataAdd(varname = "y", formula = "a + b + 2 * A", variance = 1)
  
  return(list(def = def, def2 =def2, defOut = defOut)) 
}

s_generate <- function(iter, list_of_defs) {
  
  set.seed(iter)
  list2env(list_of_defs, envir = environment())
  
  #--- add data generation code ---#
  ds <- genData(10, def, id = "site")#site
  ds <- addPeriods(ds, 11, "site", perName = "k") #create time periods for each site
  ds <- addColumns(def2, ds)
  #assign the treatment status based on the stepped-wedge design
  #per cluster trt change per wave
  ds <- trtStepWedge(ds, "site", nWaves = 10, lenWaves = 1, startPer = 1, 
                     grpName = "A", perName = "k")
  ds$site <- as.factor(ds$site)
  #30 individuals per site per period and generate each individual-level outcome
  dd <- genCluster(ds, "timeID", numIndsVar = 25, level1ID = "id")
  dd <- addColumns(defOut, dd)
  dd[, normk := (k - min(k))/(max(k) - min(k))]#scale time period into range 0-1
  
  dd[] #  generated_data is a data.table
  
}

iter=1
list_of_defs = s_define()
train_data <- s_generate(iter,list_of_defs)

fitgam <- bam(y ~ A +  s(k) + s(k, site, bs = "fs"), data = train_data, method="fREML")
res_fitgam <- c(summary(fitgam)$p.coeff["A"], summary(fitgam)$se["A"])
range <-   res_fitgam[1] + c(-1,1) * 1.96 *   res_fitgam[2]

# Fitting a more complex  model using gam
fitgam2 <- mgcv::bam(y ~ A +  s(k) + s(site, A, bs = "re") + s(k, site,  bs = "fs"), data = train_data, method="fREML")
res_fitgam2 <- c(summary(fitgam2)$p.coeff["A"], summary(fitgam2)$se["A"])
range2 <-   res_fitgam2[1] + c(-1,1) * 1.96 *   res_fitgam2[2]

# Define the knots based on the training data
knots_train <- quantile(train_data$k, probs=seq(0, 1, length=6)[-c(1, 6)])

B_train <- predict(splines::bs(train_data$k, degree=3, knots=quantile(train_data$k, probs=seq(0, 1, length=6)[-c(1, 6)])))
colnames(B_train) <- paste0("Bspline_", 1:ncol(B_train))


stan_data <- list(num_data = nrow(train_data),
                num_basis = ncol(B_train),
                B = t(B_train),
                A = train_data$A,
                y = train_data$y,
                num_sites = length(unique(train_data$site)),
                site = train_data$site
                
)

# Fit the Bayesian model
fit <- mod$sample(data = stan_data,
                refresh = 0,
                show_messages = FALSE)

diagnostics_df <- as_draws_df(fit$sampler_diagnostics())
div <- sum(diagnostics_df[, 'divergent__'])
bayes_gam = fit$summary(variables="beta_A",
                      posterior::default_summary_measures()[1:3],
                      quantiles = ~ quantile(., probs = c(0.025, 0.975)),
                      posterior::default_convergence_measures())
covered_bayes =   (bayes_gam$`2.5%`< 2 & 2 < bayes_gam$`97.5%`)

# Use the GAM model to predict on the test data to get the mean predictions.
train_data$pred_gam <- predict(fitgam, newdata = train_data)
train_data$pred_gam2 <- predict(fitgam2, newdata = train_data)

posterior_predictions <-fit$summary(variables="y_pred_mean",posterior::default_summary_measures()[1])
train_data$pred_bayesian <- posterior_predictions$mean


model_results <- data.table(est_gam_freq= res_fitgam[1], se_gam_freq=res_fitgam[2], 
                          gam_lowci=range[1], gam_upci=range[2], bayes_gam,div, 
                          est_gam_freq_rdn = res_fitgam2[1], se_gam_freq_rdn=res_fitgam2[2], 
                          gam_rdn_lowci=range2[1], gam_rdn_upci=range2[2],
                          covered_gam_freq=(range[1] < 2 & 2 < range[2]),
                          covered_bayes,covered_gam_rdn_freq=(range2[1] < 2 & 2 < range2[2])
) %>%
mutate(across(-c(variable,covered_gam_freq, covered_bayes,covered_gam_rdn_freq), round, 3))

setnames(model_results, c("est_gam_freq","se_gam_freq","lowci_freq", "upci_freq","variable","est_mean_bayes",
                        "est_med_bayes","est_sd_bayes",
                        "lowci_bayes",
                        "upci_bayes","rhat","ess_bulk","ess_tail",
                        "div","est_gam_rdn_freq","se_gam_rdn_freq","lowci_freq_rdn", "upci_freq_rdn",
                        "covered_freq","covered_bayes","covered_gam_rdn_freq"))


####frequentist plot
dd = train_data
dp <- dd[, .(avg = mean(y)), keyby = .(A, site, k)]

ggplot(data = dp, aes(x = k, y = avg)) +
  geom_line(aes(group = site, color = factor(A))) +
  scale_color_manual(values = c("#a53e3e", "#3ea5a5"),
                     guide = guide_legend(reverse = TRUE),
                     labels = c("no treatment", "treatment")) +
  ylab("average Y") +
  xlab("period (k)") +
  theme(panel.grid = element_blank(),
        legend.title = element_blank())

ggplot(data = dd, aes(x = k, y = y)) +
  geom_point(aes(color = factor(A, labels = c("Control", "Intervention"))),
             size = 0.1) +
  scale_color_manual(values = c("#d07b7c", "#7ba7d0")) +
  facet_wrap(~site, ncol = 8) +
  theme(panel.grid = element_blank(),
        legend.title = element_blank(),
        axis.text.y = element_text(size = 7),
        strip.text = element_text(size = 8),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  xlim(c(0, 24)) +
  guides(color = guide_legend(override.aes = list(size = 2)))



dp <- dd[, .(avg =  pred_bayesian), keyby = .(A, site, k)]

ggplot(data = dp, aes(x = k, y = avg)) +
  geom_line(aes(group = site, color = factor(A))) +
  scale_color_manual(values = c("#a53e3e", "#3ea5a5"),
                     guide = guide_legend(reverse = TRUE),
                     labels = c("no treatment", "treatment")) +
  ylab("average Y") +
  xlab("period (k)") +
  theme(panel.grid = element_blank(),
        legend.title = element_blank())


