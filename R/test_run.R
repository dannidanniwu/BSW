library(simstudy)
library(data.table)
library(mgcv)
library(lme4)
library(splines)
library(brms)
library(cmdstanr)
library(slurmR)
library(dplyr)
s_replicate <- function(iter) {
 
  return(data.table(iter=iter))
}

scenarios = seq(0, 2, length.out = 4)
results.aggregated2 <- vector("list", length= length(scenarios))

for (r in 1: length(scenarios)){
  iter = scenarios[r]
  
  # res <- replicate(1, s_replicate(iter=1,coefA = coefA,
  #                                 mod=mod
  #                                     ))
  
  sjob <- Slurm_lapply(1:5, 
                       FUN=s_replicate, 
                       njobs = 5, 
                       tmp_path = "/gpfs/data/troxellab/danniw/scratch", 
                       job_name = "BS_110", 
                       sbatch_opt = list(time = "4:00:00",partition = "cpu_short", `mem-per-cpu` = "8G"), 
                       plan = "wait", 
                       overwrite=TRUE) 
  res <- Slurm_collect(sjob) # data is a list
  #res<- site_plasma_all[lapply(site_plasma_all, function(x) length(x))>1] #filter out the error message
  res <- rbindlist(res) # converting list to data.table
  
  date_stamp <- gsub("-", "", Sys.Date())
  dir.create(file.path("/gpfs/data/troxellab/danniw/r/BS/", date_stamp), showWarnings = FALSE)
  save(res, file = paste0("/gpfs/data/troxellab/danniw/r/BS/", date_stamp, "/test",r,".rda"))
  results.aggregated2[[r]] <- res
}

date_stamp <- gsub("-", "", Sys.Date())
dir.create(file.path("/gpfs/data/troxellab/danniw/r/BS/", date_stamp), showWarnings = FALSE)
save(results.aggregated2, file = paste0("/gpfs/data/troxellab/danniw/r/BS/", date_stamp, "/all_test",r,".rda"))

