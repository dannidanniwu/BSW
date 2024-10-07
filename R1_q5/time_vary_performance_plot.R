library(ggplot2)
library(gridExtra)
library(ggpubr)
library(tidyr)

my_theme <- theme_classic() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 10),
    legend.position = "bottom",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10)
  )


# Define the data for TATE scenario
tate_data_2 <- data.frame(
  Measure ="TATE",
  Scenario ="2",
  Model = c("Bayesian model (3)", "Bayesian model (2)", "Monotone effect curve model"),
  Bias = abs(c(-0.07506146,-0.3083283, -7.55285)),
  RMSE = c(0.9053765, 0.9335343, 10.63354),
  LOO = c(6770, 6763, 6842)
)

# Define the data for LTE scenario
lte_data_2 <- data.frame(
  Measure ="LTE",
  Scenario ="2",
  Model = c("Bayesian model (3)", "Bayesian model (2)", "Monotone effect curve model"),
  Bias = abs(c(0.1772518, -0.2787317, -9.749612)),
  RMSE = c(1.744707, 1.698069, 14.52556),
  LOO = c(6770, 6763, 6842)
)

# Combine the data into one data frame
combined_data <- rbind(tate_data_2, lte_data_2)

combined_data$Model <- factor(combined_data$Model, levels=c("Monotone effect curve model", "Bayesian model (2)","Bayesian model (3)"))



data_subset2 <- subset(combined_data, Scenario == "2")

# Define the colors for each model
colors <- c("Bayesian model (2)" = "blue", "Monotone effect curve model" = "red", "Bayesian model (3)" = "purple")



p_bias_plot2 <- ggplot(data_subset2, aes(x = Measure, y = Bias, fill = Model)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_manual(values = colors) +
  labs( y = "Bias")+ ylim(0,10)+
  my_theme

p_rmse_plot2 <- ggplot(data_subset2, aes(x = Measure, y = RMSE, fill = Model)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_manual(values = colors) +
  labs(y = "RMSE")+ ylim(0,15)+
  my_theme


p_loo_plot2 <- ggplot(data_subset2, aes(x = Measure, y = LOO, fill = Model)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_manual(values = colors) +
  labs(y = "LOOCV") +
  coord_cartesian(ylim = c(6700, 6850)) +
  my_theme

load("C:/My PC (DESKTOP-52UHO48)/Desktop/Bayesian Spline/Bayesian stepped wedge/Bayesian stepped wedge/Bayesian-stepped-wedge/Final simulation code/Code for paper revision/output/R1_q5/time_varying_gen_m3fit_update.rda")
colnames(res) <- c("iter","ncluster","true_tate","true_lte", "variable_1", "mean_tate","median_tate",            
                   "sd_tate","2.5%_tate","97.5%_tate","rhat_tate","ess_bulk_tate","ess_tail_tate","variable", "mean_lte",              
                   "median_lte","sd_lte","2.5%_lte","97.5%_lte","rhat_lte","ess_bulk_lte","ess_tail_lte","covered_bayes_tate",
                   "covered_bayes_lte","div","num_max_tree_depth","elpd_loo",
                   "elpd_se", "p_loo","p_loo_se","looic", "looic_se")

model3_div <- res$div
m3_loo <- res$looic

m3_tate <- res$median_tate
m3_lte <- res$median_lte

true_tate = unique(res$true_tate)
true_lte = unique(res$true_lte)

rm(res)
load("C:/My PC (DESKTOP-52UHO48)/Desktop/Bayesian Spline/Bayesian stepped wedge/Bayesian stepped wedge/Bayesian-stepped-wedge/Final simulation code/Data/res_inc_decrease_t_vary_m2.rda")
colnames(res) <- c("iter","ncluster","true_tate","true_lte", "variable_1", "mean_tate","median_tate",            
                                    "sd_tate","2.5%_tate","97.5%_tate","rhat_tate","ess_bulk_tate","ess_tail_tate","variable", "mean_lte",              
                                    "median_lte","sd_lte","2.5%_lte","97.5%_lte","rhat_lte","ess_bulk_lte","ess_tail_lte","covered_bayes_tate",
                                    "covered_bayes_lte","gam_tate","gam_lte","div","num_max_tree_depth","elpd_loo",
                                    "elpd_se", "p_loo","p_loo_se","looic", "looic_se")
model2_div <- res$div
m2_loo <- res$looic

m2_tate <- res$median_tate
m2_lte <- res$median_lte

rm(res)
load("C:/My PC (DESKTOP-52UHO48)/Desktop/Bayesian Spline/Bayesian stepped wedge/Bayesian stepped wedge/Bayesian-stepped-wedge/Final simulation code/Data/res_dec_increase_t_vary_mono.rda")
colnames(res) <- c("iter","ncluster","true_tate","true_lte", "variable_1", "mean_tate","median_tate",            
                                    "sd_tate","2.5%_tate","97.5%_tate","rhat_tate","ess_bulk_tate","ess_tail_tate","variable", "mean_lte",              
                                    "median_lte","sd_lte","2.5%_lte","97.5%_lte","rhat_lte","ess_bulk_lte","ess_tail_lte","covered_bayes_tate",
                                    "covered_bayes_lte","div","num_max_tree_depth","elpd_loo",
                                    "elpd_se", "p_loo","p_loo_se","looic", "looic_se")
monomod_div <- res$div
mono_loo <- res$looic

mono_tate <- res$median_tate
mono_lte <- res$median_lte



#estimation
model <- c(rep("m3", length(m3_tate)), rep("m2", length(m2_tate)), rep("mono", length(mono_tate)))
estimation <- c(rep("TATE", length(m3_tate) + length(m2_tate) + length(mono_tate)),
                rep("LTE", length(m3_lte) + length(m2_lte) + length(mono_lte)))
value <- c(m3_tate, m2_tate, mono_tate, m3_lte, m2_lte, mono_lte)

data <- data.frame(Model = factor(model, levels = c( "mono","m2", "m3")),
                   Estimation = factor(estimation, levels = c("LTE","TATE")),
                   Value = value)

p_est <- ggplot(data, aes(x = Estimation, y = Value, fill = Model)) +
  geom_boxplot(outlier.shape = 4) +
  geom_hline(aes(yintercept = true_tate), linetype = "dashed", color = "red", size = 1, show.legend = TRUE, labels = 'TATE') +
  geom_hline(aes(yintercept = true_lte), linetype = "dashed", color = "blue", size = 1, show.legend = TRUE, labels ='LTE') +
  labs(
       x = "Estimation",
       y = "Posterior Median",
       fill = "Model")+my_theme
rm(model)
rm(value)
######################

model <- c(rep("m3", length(model3_div)), rep("m2", length(model2_div)), rep("mono", length(monomod_div)))
value <- c(model3_div,model2_div,monomod_div)*100/8000

df_div <- data.frame(Model = factor(model, levels = c( "mono","m2", "m3")),
                   value = value)
# Create the box plot
p_div <- ggplot(df_div, aes(x = Model, y = value, fill=Model)) +
  geom_boxplot(outlier.shape = 4) +
  geom_hline(yintercept = 1, color = "red", linetype = "dashed") +
  labs(
       x = "Model",
       y = "Divergent Transitions(%)",
       fill = "Model")+my_theme
rm(model)
rm(value)
#########################################

model <- c(rep("m3", length(m3_loo)), rep("m2", length(m2_loo)), rep("mono", length(mono_loo)))
value <- c(m3_loo,m2_loo,mono_loo)

df_long <- data.frame(Model = factor(model, levels = c( "mono","m2", "m3")),
                      value = value)

# Create the box plot
p_loo <- ggplot(df_long, aes(x = Model, y = value, fill=Model)) +
  geom_boxplot(outlier.shape = 4) +
  labs(
    x = "Model",
    y = "LOOIC",
    fill = "Model")+my_theme 

# Arrange the plots in a grid
combined_plot <- ggarrange(p_rmse_plot2,p_bias_plot2, p_est, p_loo, nrow=2, ncol=2, common.legend = TRUE, legend="bottom")
combined_plot

ggsave(
  filename = "./m3_genscen2.eps",
  plot = combined_plot,
  width = 11.969,       # width in inches
  height = 11.969,      # height in inches
  units = "in",         # units for width and height
  dpi = 3000            # high resolution
)