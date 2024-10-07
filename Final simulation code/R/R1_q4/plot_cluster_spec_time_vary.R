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


tate_data <- data.frame(
  Measure = rep("TATE", 4),
  Scenario = rep("1", 4),
  Model = c("Bayesian model (6)", "Monotone effect curve model", "Bayesian model (7)","NCS"),
  LOOIC = c(6781, 6920, 6770,NA),
  Bias = abs(c(-0.4427766, -8.718332, -0.02817164,0.03514932)),
  RMSE = c(1.108692, 11.15854,0.7057979,0.8318129)
)


lte_data <- data.frame(
  Measure = rep("LTE", 4),
  Scenario = rep("1", 4),
  Model = c("Bayesian model (6)", "Monotone effect curve model", "Bayesian model (7)","NCS"),
  LOOIC = c(6781, 6920, 6770,NA),
  Bias = abs(c(-0.4526795, -11.50264, 0.2896472,0.2334521)),
  RMSE = c(2.049236, 15.42879, 1.387395,1.598591)
)


combined_data <- rbind(tate_data, lte_data)
combined_data$Model <- factor(combined_data$Model, levels=c("Monotone effect curve model", "Bayesian model (6)","NCS","Bayesian model (7)"))
# Define the colors for each model
colors <- c("Bayesian model (6)" = "blue", "Monotone effect curve model" = "red", "Bayesian model (7)" = "black", "NCS"="purple")

# Create a plot for each metric
p_bias <- ggplot(combined_data, aes(x = Measure, y = Bias, fill = Model)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_manual(values = colors) +
  labs( y = "Bias")+ 
  my_theme

p_rmse <- ggplot(combined_data, aes(x = Measure, y = RMSE, fill = Model)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_manual(values = colors) +
  labs(y = "RMSE")+ 
  my_theme

# p_loo <- ggplot(combined_data, aes(x = Measure, y = LOOIC, fill = Model)) +
#   geom_bar(stat = "identity", position = position_dodge()) +
#   scale_fill_manual(values = colors) +
#   labs(y = "LOOCV") +
#   coord_cartesian(ylim = c(6750, 6950)) +
#   my_theme

###bayes spf time vary model
load("C:/My PC (DESKTOP-52UHO48)/Desktop/Bayesian Spline/Bayesian stepped wedge/Bayesian stepped wedge/Bayesian-stepped-wedge/Final simulation code/Data/res_dec_increase_cluster_specific_t_vary_m3.rda")
colnames(res) <- c( 
  "iter","ncluster","true_tate","true_lte", "variable_1", "mean_tate","median_tate",            
  "sd_tate","2.5%_tate","97.5%_tate","rhat_tate","ess_bulk_tate","ess_tail_tate","variable", "mean_lte",              
  "median_lte","sd_lte","2.5%_lte","97.5%_lte","rhat_lte","ess_bulk_lte","ess_tail_lte","covered_bayes_tate",
  "covered_bayes_lte","gam_tate","gam_lte","div","num_max_tree_depth","elpd_loo",
  "elpd_se", "p_loo","p_loo_se","looic", "looic_se")
m3_loo <- res$looic

m3_tate <- res$median_tate
m3_lte <- res$median_lte
rm(res)

load("C:/My PC (DESKTOP-52UHO48)/Desktop/Bayesian Spline/Bayesian stepped wedge/Bayesian stepped wedge/Bayesian-stepped-wedge/Final simulation code/Data/res_dec_increase_cluster_specific_t_vary_m2.rda")
colnames(res) <- c( 
  "iter","ncluster","true_tate","true_lte", "variable_1", "mean_tate","median_tate",            
  "sd_tate","2.5%_tate","97.5%_tate","rhat_tate","ess_bulk_tate","ess_tail_tate","variable", "mean_lte",              
  "median_lte","sd_lte","2.5%_lte","97.5%_lte","rhat_lte","ess_bulk_lte","ess_tail_lte","covered_bayes_tate",
  "covered_bayes_lte","div","num_max_tree_depth","elpd_loo",
  "elpd_se", "p_loo","p_loo_se","looic", "looic_se")
true_tate = 3.033328
true_lte = 2.5
m2_loo <- res$looic

m2_tate <- res$median_tate
m2_lte <- res$median_lte
rm(res)

######mono model
load("C:/My PC (DESKTOP-52UHO48)/Desktop/Bayesian Spline/Bayesian stepped wedge/Bayesian stepped wedge/Bayesian-stepped-wedge/Final simulation code/Data/res_cluster_spf_t_vary_our_mono.rda")
colnames(res) <- c("iter","ncluster","true_tate","true_lte", "variable_1", "mean_tate","median_tate",            
                   "sd_tate","2.5%_tate","97.5%_tate","rhat_tate","ess_bulk_tate","ess_tail_tate","variable", "mean_lte",              
                   "median_lte","sd_lte","2.5%_lte","97.5%_lte","rhat_lte","ess_bulk_lte","ess_tail_lte","covered_bayes_monot_tate",
                   "covered_bayes_monot_lte","div_monot","num_max_tree_depth_monot","elpd_loo",
                   "elpd_se", "p_loo","p_loo_se","looic", "looic_se")

mono_loo <- res$looic

mono_tate <- res$median_tate
mono_lte <- res$median_lte
rm(res)

load("C:/My PC (DESKTOP-52UHO48)/Desktop/Bayesian Spline/Bayesian stepped wedge/Bayesian stepped wedge/Bayesian-stepped-wedge/Final simulation code/Code for paper revision/output/R1_q4/inc_dec_time_cluster_spc_vary_trt_ncs.rda")

ncs_tate <- res$ncs_tate
ncs_lte <- res$ncs_lte
ncs_loo <- NA

model <- c(rep("ncs", length(ncs_tate)), rep("m2", length(m2_tate)), rep("mono", length(mono_tate)),rep("m3", length(m3_tate)))
estimation <- c(rep("TATE", length(ncs_tate) + length(m2_tate) + length(mono_tate)+length(m3_tate)),
                rep("LTE", length(ncs_lte) + length(m2_lte) + length(mono_lte)+length(m3_lte)))
value <- c(ncs_tate, m2_tate, mono_tate,m3_tate, ncs_lte, m2_lte, mono_lte,m3_lte)

data <- data.frame(Model = factor(model, levels = c( "mono","m2","ncs","m3")),
                   Estimation = factor(estimation, levels = c("LTE","TATE")),
                   Value = value)

p_est_plot <- ggplot(data, aes(x = Estimation, y = Value, fill = Model)) +
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
model <- c( rep("m2", length(m2_loo)), rep("mono", length(mono_loo)),rep("m3", length(m3_loo)))

value <- c(m2_loo,mono_loo,m3_loo)

df_long <- data.frame(Model = factor(model, levels = c( "mono","m2","m3")),
                      value = value)

p_loo_plot <- ggplot(df_long, aes(x = Model, y = value, fill=Model)) +
  geom_boxplot(outlier.shape = 4) +
  labs(
    x = "Model",
    y = "LOOIC",
    fill = "Model")+my_theme+ylim(6500,7100)

rm(model)
rm(value)


# Arrange the plots in a grid
combined_plot <- ggarrange( p_rmse,p_bias,p_est_plot,p_loo_plot,  nrow=2, ncol=2, common.legend = TRUE, legend="bottom")
combined_plot

ggsave(
  filename = "./cluster_spf_t_vary_eff_perform.eps",
  plot = combined_plot,
  width = 11.969,       # width in inches
  height = 11.969,      # height in inches
  units = "in",         # units for width and height
  dpi = 3000            # high resolution
)