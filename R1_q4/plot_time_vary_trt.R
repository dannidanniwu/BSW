library(ggplot2)
library(gridExtra)
library(ggpubr)

my_theme <- theme_classic() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 10),
    legend.position = "bottom",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10)
  )
# Assuming data is entered as follows (replace with actual values):
tate_data <- data.frame(
  Measure ="TATE",
  Scenario ="1",
  Model = c("Bayesian model (2)", "Monotone effect curve model","NCS"),
  Bias = abs(c(-0.2561728, -0.1420227,0.01368202)),
  RMSE = c(0.7913839, 1.297454, 0.5351386),
  LOO = c(5246, 5261,NA)
)

lte_data <- data.frame(
  Measure ="LTE",
  Scenario ="1",
  Model = c("Bayesian model (2)", "Monotone effect curve model","NCS"),
  Bias = abs(c(-0.2765687, 0.7749753, 0.5420786)),
  RMSE = c(1.454923, 1.95226, 1.122187),
  LOO = c(5246, 5261, NA)
)

# Define the data for TATE scenario
tate_data_2 <- data.frame(
  Measure ="TATE",
  Scenario ="2",
  Model = c("Bayesian model (2)", "Monotone effect curve model","NCS"),
  Bias = abs(c(-0.3083283, -7.55285,0.0665564)),
  RMSE = c(0.9335343, 10.63354,0.5700118),
  LOO = c(6763, 6842, NA)
)

# Define the data for LTE scenario
lte_data_2 <- data.frame(
  Measure ="LTE",
  Scenario ="2",
  Model = c("Bayesian model (2)", "Monotone effect curve model","NCS"),
  Bias = abs(c(-0.2787317, -9.749612, 0.4019669)),
  RMSE = c(1.698069, 14.52556, 1.179419),
  LOO = c(6763, 6842, NA)
)

# Combine the data into one data frame
combined_data <- rbind(tate_data, lte_data, tate_data_2, lte_data_2)

combined_data$Model <- factor(combined_data$Model, levels=c("Monotone effect curve model","NCS","Bayesian model (2)"))


data_subset1 <- subset(combined_data, Scenario == "1")
data_subset2 <- subset(combined_data, Scenario == "2")

# Define the colors for each model
colors <- c("Bayesian model (2)" = "blue", "Monotone effect curve model" = "orange","NCS" = "red")

# Create a plot for each metric
p_bias_plot1 <- ggplot(data_subset1, aes(x = Measure, y = Bias, fill = Model)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_manual(values = colors) +
  labs( y = "Bias")+ ylim(0,10)+
  my_theme

p_bias_plot2 <- ggplot(data_subset2, aes(x = Measure, y = Bias, fill = Model)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_manual(values = colors) +
  labs( y = "Bias")+ ylim(0,10)+
  my_theme


p_rmse_plot1 <- ggplot(data_subset1, aes(x = Measure, y = RMSE, fill = Model)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_manual(values = colors) +
  labs(y = "RMSE")+ ylim(0,15)+
  my_theme

p_rmse_plot2 <- ggplot(data_subset2, aes(x = Measure, y = RMSE, fill = Model)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_manual(values = colors) +
  labs(y = "RMSE")+ ylim(0,15)+
  my_theme

# Create separate plots
# p_loo_plot1 <- ggplot(data_subset1, aes(x = Measure, y = LOO, fill = Model)) +
#   geom_bar(stat = "identity", position = position_dodge()) +
#   scale_fill_manual(values = colors) +
#   labs(y = "LOOCV") +
#   coord_cartesian(ylim = c(5200, 5280)) +
#   my_theme
# 
# p_loo_plot2 <- ggplot(data_subset2, aes(x = Measure, y = LOO, fill = Model)) +
#   geom_bar(stat = "identity", position = position_dodge()) +
#   scale_fill_manual(values = colors) +
#   labs(y = "LOOCV") +
#   coord_cartesian(ylim = c(6700, 6850)) +
#   my_theme

load("C:/My PC (DESKTOP-52UHO48)/Desktop/Bayesian Spline/Bayesian stepped wedge/Bayesian stepped wedge/Bayesian-stepped-wedge/Final simulation code/Data/res_increase_t_vary_mono.rda")
colnames(res) <- c( "iter","coefA","ncluster","true_tate","true_lte", "variable_1", "mean_tate","median_tate",            
                    "sd_tate","2.5%_tate","97.5%_tate","rhat_tate","ess_bulk_tate","ess_tail_tate","variable", "mean_lte",              
                    "median_lte","sd_lte","2.5%_lte","97.5%_lte","rhat_lte","ess_bulk_lte","ess_tail_lte","covered_bayes_monot_tate",
                    "covered_bayes_monot_lte","div_monot","num_max_tree_depth_monot","elpd_loo",
                    "elpd_se", "p_loo","p_loo_se","looic", "looic_se")

mono_loo <- res$looic

mono_tate <- res$median_tate
mono_lte <- res$median_lte

##m2 
load("C:/My PC (DESKTOP-52UHO48)/Desktop/Bayesian Spline/Bayesian stepped wedge/Bayesian stepped wedge/Bayesian-stepped-wedge/Final simulation code/Data/res_increase_t_vary_our_m2.rda")
colnames(res) <- c( "iter","coefA","ncluster","true_tate","true_lte", "variable_1", "mean_tate","median_tate",            
                    "sd_tate","2.5%_tate","97.5%_tate","rhat_tate","ess_bulk_tate","ess_tail_tate","variable", "mean_lte",              
                    "median_lte","sd_lte","2.5%_lte","97.5%_lte","rhat_lte","ess_bulk_lte","ess_tail_lte","covered_bayes_tate",
                    "covered_bayes_lte","gam_tate","gam_lte","div","num_max_tree_depth","elpd_loo",
                    "elpd_se", "p_loo","p_loo_se","looic", "looic_se")
# true_tate = unique(res$true_tate)
# true_lte = unique(res$true_lte)
true_tate = 4.738699
true_lte = 4.999997
m2_loo <- res$looic

m2_tate <- res$median_tate
m2_lte <- res$median_lte

##ncs
load("C:/My PC (DESKTOP-52UHO48)/Desktop/Bayesian Spline/Bayesian stepped wedge/Bayesian stepped wedge/Bayesian-stepped-wedge/Final simulation code/Code for paper revision/output/R1_q4/inc_time_vary_trt_ncs.rda")

summary(res)
ncs_tate <- res$ncs_tate
ncs_lte <- res$ncs_lte
ncs_loo <- NA

#estimation
model <- c(rep("ncs", length(ncs_tate)), rep("m2", length(m2_tate)), rep("mono", length(mono_tate)))
estimation <- c(rep("TATE", length(ncs_tate) + length(m2_tate) + length(mono_tate)),
                rep("LTE", length(ncs_lte) + length(m2_lte) + length(mono_lte)))
value <- c(ncs_tate, m2_tate, mono_tate, ncs_lte, m2_lte, mono_lte)

data <- data.frame(Model = factor(model, levels = c( "mono", "ncs","m2")),
                   Estimation = factor(estimation, levels = c("LTE","TATE")),
                   Value = value)

p_est_plot1 <- ggplot(data, aes(x = Estimation, y = Value, fill = Model)) +
  geom_boxplot(outlier.shape = 4) +
  geom_hline(aes(yintercept = 4.738699), linetype = "dashed", color = "red", size = 1, show.legend = TRUE) +
  geom_hline(aes(yintercept = 4.999997), linetype = "dashed", color = "blue", size = 1, show.legend = TRUE) +
  labs(
    x = "Estimation",
    y = "Posterior Median",
    fill = "Model")+my_theme
rm(model)
rm(value)
rm(estimation)
rm(data)
#looic
model <- c(rep("ncs", length(ncs_loo)), rep("m2", length(m2_loo)), rep("mono", length(mono_loo)))

value <- c(ncs_loo,m2_loo,mono_loo)

df_long <- data.frame(Model = factor(model, levels = c( "mono","ncs","m2")),
                      value = value)

p_loo_plot1 <- ggplot(df_long, aes(x = Model, y = value, fill=Model)) +
  geom_boxplot(outlier.shape = 4) +
  labs(
    x = "Model",
    y = "LOOIC",
    fill = "Model")+my_theme 

rm(model)
rm(value)
rm(res)
rm(df_long)
#dec inc data generation
#mono
load("C:/My PC (DESKTOP-52UHO48)/Desktop/Bayesian Spline/Bayesian stepped wedge/Bayesian stepped wedge/Bayesian-stepped-wedge/Final simulation code/Data/res_dec_increase_t_vary_mono.rda")
colnames(res) <- c( "iter","ncluster","true_tate","true_lte", "variable_1", "mean_tate","median_tate",            
                    "sd_tate","2.5%_tate","97.5%_tate","rhat_tate","ess_bulk_tate","ess_tail_tate","variable", "mean_lte",              
                    "median_lte","sd_lte","2.5%_lte","97.5%_lte","rhat_lte","ess_bulk_lte","ess_tail_lte","covered_bayes_monot_tate",
                    "covered_bayes_monot_lte","div_monot","num_max_tree_depth_monot","elpd_loo",
                    "elpd_se", "p_loo","p_loo_se","looic", "looic_se")

mono_loo <- res$looic

mono_tate <- res$median_tate
mono_lte <- res$median_lte
rm(res)
#m2
load("C:/My PC (DESKTOP-52UHO48)/Desktop/Bayesian Spline/Bayesian stepped wedge/Bayesian stepped wedge/Bayesian-stepped-wedge/Final simulation code/Data/res_inc_decrease_t_vary_m2.rda")
colnames(res) <- c("iter","ncluster","true_tate","true_lte", "variable_1", "mean_tate","median_tate",            
                   "sd_tate","2.5%_tate","97.5%_tate","rhat_tate","ess_bulk_tate","ess_tail_tate","variable", "mean_lte",              
                   "median_lte","sd_lte","2.5%_lte","97.5%_lte","rhat_lte","ess_bulk_lte","ess_tail_lte","covered_bayes_tate",
                   "covered_bayes_lte","gam_tate","gam_lte","div","num_max_tree_depth","elpd_loo",
                   "elpd_se", "p_loo","p_loo_se","looic", "looic_se")
# true_tate = unique(res$true_tate)
# true_lte = unique(res$true_lte)
true_tate = 3.033328
true_lte = 2.5
m2_loo <- res$looic

m2_tate <- res$median_tate
m2_lte <- res$median_lte
rm(res)
#ncs
load("C:/My PC (DESKTOP-52UHO48)/Desktop/Bayesian Spline/Bayesian stepped wedge/Bayesian stepped wedge/Bayesian-stepped-wedge/Final simulation code/Code for paper revision/output/R1_q4/inc_dec_time_vary_trt_ncs.rda")
summary(res)
ncs_tate <- res$ncs_tate
ncs_lte <- res$ncs_lte
ncs_loo <- NA

#estimation
model <- c(rep("ncs", length(ncs_tate)), rep("m2", length(m2_tate)), rep("mono", length(mono_tate)))
estimation <- c(rep("TATE", length(ncs_tate) + length(m2_tate) + length(mono_tate)),
                rep("LTE", length(ncs_lte) + length(m2_lte) + length(mono_lte)))
value <- c(ncs_tate, m2_tate, mono_tate, ncs_lte, m2_lte, mono_lte)

data <- data.frame(Model = factor(model, levels = c( "mono", "ncs","m2")),
                   Estimation = factor(estimation, levels = c("LTE","TATE")),
                   Value = value)

p_est_plot2 <- ggplot(data, aes(x = Estimation, y = Value, fill = Model)) +
  geom_boxplot(outlier.shape = 4) +
  geom_hline(aes(yintercept = 3.033328), linetype = "dashed", color = "red", size = 1, show.legend = TRUE) +
  geom_hline(aes(yintercept = 2.5), linetype = "dashed", color = "blue", size = 1, show.legend = TRUE) +
  labs(
    x = "Estimation",
    y = "Posterior Median",
    fill = "Model")+my_theme
rm(model)
rm(value)

#looic
model <- c(rep("ncs", length(ncs_loo)), rep("m2", length(m2_loo)), rep("mono", length(mono_loo)))

value <- c(ncs_loo,m2_loo,mono_loo)

df_long <- data.frame(Model = factor(model, levels = c( "mono","ncs","m2")),
                      value = value)

p_loo_plot2 <- ggplot(df_long, aes(x = Model, y = value, fill=Model)) +
  geom_boxplot(outlier.shape = 4) +
  labs(
    x = "Model",
    y = "LOOIC",
    fill = "Model")+my_theme 

rm(model)
rm(value)

# Arrange the plots in a grid
combined_plot <- ggarrange(p_rmse_plot1,p_rmse_plot2,p_bias_plot1,p_bias_plot2,p_est_plot1,p_est_plot2,p_loo_plot1, p_loo_plot2, nrow=4, ncol=2, common.legend = TRUE, legend="bottom")
combined_plot

ggsave(
  filename = "./t_vary_eff_perform.eps",
  plot = combined_plot,
  width = 11.969,       # width in inches
  height = 11.969,      # height in inches
  units = "in",         # units for width and height
  dpi = 3000            # high resolution
)