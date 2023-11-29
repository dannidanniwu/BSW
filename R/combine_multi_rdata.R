# Set the working directory to the folder containing the RData files
setwd("~/Library/CloudStorage/OneDrive-NYULangoneHealth/Bayesian stepped wedge/Bayesian-stepped-wedge/Data/cluster_spf_t_vary_m1")

# Initialize an empty list
data_list <- list()

# Loop through the files
for (i in 1:150) {
  file_name <- paste0("./output_iter_", i, ".RData")
  load(file_name)
  # Assuming each file contains one main variable, you append it to the list
  # Replace 'your_data_variable' with the actual name of the variable in your .rda files
  data_list[[i]] <- result
}

# Now data_list contains the data from all RData files

res <- rbindlist(data_list)
summary(res)
save(res,file="res_cluster_spf_t_vary_our_m.rda")
