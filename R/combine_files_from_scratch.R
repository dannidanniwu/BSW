#got data after running interim.R
library("collapse")
library("dplyr")
library("data.table")
library("tidyverse")

##sim16
site_plasma_all <- list.files(pattern = "^03-answer") %>%
  map(readRDS)

plasma <- lapply(site_plasma_all, function(x) Filter(Negate(is.null), x))
p2 <- lapply(plasma,function(x) x[lapply(x,length)>1])
p3 <- p2[lapply(p2,function(x) length(x))>1]
res<- unlist2d(p3, idcols = "replicate",DT = TRUE)#the first colums: the first level list;the second column:the second level list
save(res,file="C:/Users/Danni/OneDrive - NYU Langone Health/Bayesian TBI/04092022/n1600/3y_5cov_1600_betas.rda")

#save(res,file="./6y_10cov_n200_iter100_comsym.rda")
#save(res, file = "./mixed_hier_ord_bi_3c_2b_cov_rdn_gen_v3.rda")
save(res,file="./stepped_wedge_cluster_specific_treatment_effect_24cluster.rda")
