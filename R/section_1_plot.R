# Load the necessary libraries
library(ggplot2)
library(reshape2)

# Input the data for the first table (Coverage of true value)
coverage_data <- data.frame(
  TreatmentEffect = c(0, 0.25, 0.5, 0.75, 1),
  `Model (2)` = c(78.808, 76.812, 76.923, 75.510, 76.433),
  `Model (3)` = c(87.417, 86.232, 87.821, 86.395, 87.898),
  `Model (4)` = c(0, 0, 0, 0, 0),
  `Model (5)` = c(90.728, 90.580, 89.744, 89.796, 90.446),
  `Model (6)` = c(94.040, 94.203, 93.590, 93.197, 93.631),
  `Bayesian Model (1a)` = c(99.338, 98.551, 99.359, 98.639, 99.363),
  `Bayesian Model (1b)` = c(98.675, 98.551, 99.359, 98.639, 98.726)
)

# Convert to long format for plotting
coverage_long <- melt(coverage_data, id.vars = 'TreatmentEffect')

# Plot the data
ggplot(coverage_long, aes(x = TreatmentEffect, y = value, colour = variable)) +
  geom_line(size = 1.5) +  # Increase line size
  geom_point(size = 3) +   # Increase point size
  scale_color_manual(values = c("red", "blue", "green", "orange", "purple", "brown", "cyan"),
                     labels = c("Model (2)","Model (3)", "Model (4)", "Model (5)", "Model (6)", "Bayesian Model (1a)", "Bayesian Model (1b)")) +
  labs(title = "Coverage of True Value", x = "Treatment Effect", y = "Coverage (%)") +
  theme_minimal()+
  theme(
    axis.title.x = element_text(size = 14),   # Increase x-axis label size
    axis.title.y = element_text(size = 14),   # Increase y-axis label size
    legend.title = element_text(size = 12),   # Increase legend title size
    legend.text = element_text(size = 12)     # Increase legend text size
  )

# Repeat similar steps for the other two tables (RMSE and the first table if needed)

# Load the necessary libraries
library(ggplot2)
library(reshape2)

# Input the data for the RMSE table
# Input the data for the RMSE table
rmse_data <- data.frame(
  TreatmentEffect = c(0, 0.25, 0.5, 0.75, 1),
  `Model (2)` = c(0.177, 0.182, 0.182, 0.184, 0.180),
  `Model (3)` = c(0.139, 0.147, 0.144, 0.146, 0.143),
  `Model (4)` = c(3.139, 3.125, 3.119, 3.137, 3.119),
  `Model (5)` = c(0.127, 0.136, 0.135, 0.135, 0.134),
  `Model (6)` = c(0.126, 0.135, 0.133, 0.133, 0.132),
  `Bayesian Model (1a)` = c(0.161, 0.168, 0.166, 0.162, 0.166),
  `Bayesian Model (1b)` = c(0.161, 0.168, 0.166, 0.162, 0.166)
)

# Convert to long format for plotting
rmse_long <- melt(rmse_data, id.vars = 'TreatmentEffect')

# Plot the data using ggplot2
ggplot(rmse_long, aes(x = TreatmentEffect, y = value, colour = variable)) +
  geom_line(size = 1.5) +  # Increase line size
  geom_point(size = 3) +   # Increase point size
  scale_color_manual(values = c("red", "blue", "green", "orange", "purple", "brown", "cyan"),
                     labels = c("Model (2)","Model (3)", "Model (4)", "Model (5)", "Model (6)", "Bayesian Model (1a)", "Bayesian Model (1b)")) +
  labs(title = "RMSE of Estimations", x = "Treatment Effect", y = "RMSE") +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 14),   # Increase x-axis label size
    axis.title.y = element_text(size = 14),   # Increase y-axis label size
    legend.title = element_text(size = 12),   # Increase legend title size
    legend.text = element_text(size = 12)     # Increase legend text size
  )

# Run the plotting code to generate the plot


# Note: Replace '... enter RMSE data here' with the actual data from your RMSE table.

# Enter the data for the first table (power)
power_data <- data.frame(
  Treatment_effect = c(0, 0.25, 0.5, 0.75, 1),
  `Model (2)` = c(21.192, 54.348, 94.231, 100, 100),
  `Model (3)` = c(12.583, 57.246, 97.436, 100, 100),
  `Model (4)` = c(100, 0, 0, 0, 0),
  `Model (5)` = c(9.272, 51.44, 98.71, 100, 100),
  `Model (6)` = c(5.96, 42.754, 93.59, 95.238, 100),
  `Bayesian Model (1a)` = c(0.662, 4.348, 46.795, 95.238, 100),
  `Bayesian Model (1b)` = c(1.325, 6.522, 46.154, 93.197, 100)
)

# Convert the data to long format for plotting with ggplot2
library(reshape2)
power_long <- melt(power_data, id.vars = 'Treatment_effect', variable.name = 'Model', value.name = 'Power')

# Plotting power
library(ggplot2)
# Create the ggplot object
ggplot(power_long, aes(x = Treatment_effect, y = Power, color = Model)) +
  geom_line(size = 1.5) +  # Increase line size
  geom_point(size = 3) +   # Increase point size
  ggtitle('Power and Type 1 error rates') +          # Increase title size
  scale_color_manual(values = c("red", "blue", "green", "orange", "purple", "brown", "cyan"),
                     labels = c("Model (2)","Model (3)", "Model (4)", "Model (5)", "Model (6)", "Bayesian Model (1a)", "Bayesian Model (1b)")) +
  xlab('Treatment Effect') +
  ylab('% of simulations with 95% CrI/CIs not including zero') +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 14),   # Increase x-axis label size
    axis.title.y = element_text(size = 14),   # Increase y-axis label size
    legend.title = element_text(size = 12),   # Increase legend title size
    legend.text = element_text(size = 12)     # Increase legend text size
  )

# Repeat similar steps for the other two tables (coverage and RMSE)

