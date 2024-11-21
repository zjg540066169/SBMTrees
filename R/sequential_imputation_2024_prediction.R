library(Rcpp)
library(RcppParallel)
setThreadOptions(stackSize = 3e+8)
library(lme4)
#library(optimx)
library(dplyr)
library(readr)
sourceCpp("./src/BMTrees_MCMC.cpp")


source("R/simulate.R")
#source("R/evaluation.R")
#sourceCpp("./src/BMTrees_MCMC.cpp")

# model
# 1 nDP residual
# 2 cDP random effects
# 3 double DP



folder_name = function(nonlinear, nonrandeff, nonresidual){
   fold = ""
   if(nonlinear){
     fold = paste0(fold, "nonlinear_")
   }else{
     fold = paste0(fold, "linear_")
   }
   if(nonrandeff){
     fold = paste0(fold, "nonrandeff_")
   }else{
     fold = paste0(fold, "randeff_")
   }
   if(nonresidual){
     fold = paste0(fold, "nonresidual")
   }else{
     fold = paste0(fold, "residual")
   }
   fold
}



args = commandArgs(trailingOnly=TRUE)

ags = as.integer(args[1])
ags = 1
seed = as.integer(ags[1])
print(paste0("seed:", seed))
path_ = "./simulation_prediction/"
i = as.integer(ags)[1]
dir.create(paste0(path_, i))
scenario = paste0(path_, i, "/")

path = paste0(scenario, folder_name(T, T, T))
dir.create(path)
data = simulation(1000, n_obs_per_sub = 5, nonlinear = T, nonrandeff = T, nonresidual = T, alligned = F, seed = i)
X = data$X_O
subject_id = data$subject_id
Y = data$Y_O
Z = data$Z

#write_csv(data$X_O, paste0(path, "/X_O.csv"))
#write_csv(data$X_mis, paste0(path, "/X_mis.csv"))
write_csv(as.data.frame(data$Y), paste0(path, "/Y.csv"))
write_csv(as.data.frame(data$Y_O), paste0(path, "/Y_O.csv"))

#write_csv(as.data.frame(data$Z), paste0(path, "/Z.csv"))
#write_csv(as.data.frame(data$subject_id), paste0(path, "/subject_id.csv"))
#write_csv(as.data.frame(data$trajectory), paste0(path, "/trajectory.csv"))



nburn = 3000
npost = 4000
skip = 200


pred = BMTrees_MCMC(!is.na(data$Y), nburn, npost, Y, as.matrix(X),Z, subject_id, FALSE, F, NULL, 1e-40, CDP_residual = F, CDP_re = F, 0, 0.99)












dir.create(paste0(path, "/0/"))
write_csv(as.data.frame(pred$post_y_sample_test), paste0(path, "/0/test.csv"))
write_csv(as.data.frame(pred$post_y_sample), paste0(path, "/0/train.csv"))


pred = BMTrees_MCMC(!is.na(data$Y), nburn, npost, Y, as.matrix(X),Z, subject_id, FALSE, T, NULL, 1e-40, CDP_residual = T, CDP_re = F, 0, 0.99)

dir.create(paste0(path, "/1/"))
write_csv(as.data.frame(pred$post_y_sample_test), paste0(path, "/1/test.csv"))
write_csv(as.data.frame(pred$post_y_sample), paste0(path, "/1/train.csv"))


pred = BMTrees_MCMC(!is.na(data$Y), nburn, npost, Y, as.matrix(X),Z, subject_id, FALSE, T, NULL, 1e-40, CDP_residual = F, CDP_re = T, 0, 0.99)

dir.create(paste0(path, "/2/"))
write_csv(as.data.frame(pred$post_y_sample_test), paste0(path, "/2/test.csv"))
write_csv(as.data.frame(pred$post_y_sample), paste0(path, "/2/train.csv"))


pred = BMTrees_MCMC(!is.na(data$Y), nburn, npost, Y, as.matrix(X),Z, subject_id, FALSE, T, NULL, 1e-40, CDP_residual = T, CDP_re = T, 0, 0.99)

dir.create(paste0(path, "/3/"))
write_csv(as.data.frame(pred$post_y_sample_test), paste0(path, "/3/test.csv"))
write_csv(as.data.frame(pred$post_y_sample), paste0(path, "/3/train.csv"))



#imp = sequential_imputation_R(X, Y, type, Z, subject_id, model = 3, nburn = nburn, npost = npost, n_model_burn = 0L, skip = skip, verbose = TRUE, seed = seed, tol = 1e-20, ncores = 0)
#write_imputed_data(imp$imputed_data, paste0(path, "/3/"))

#imp = sequential_imputation_R(X, Y, type, Z, subject_id, model = 1, nburn = nburn, npost = npost, n_model_burn = 0L, skip = skip, verbose = TRUE, seed = seed, tol = 1e-20, ncores = 0)
#write_imputed_data(imp$imputed_data, paste0(path, "/1/"))

#imp = sequential_imputation_R(X, Y, type, Z, subject_id, model = 3, nburn = nburn, npost = npost, n_model_burn = 0L, skip = skip, verbose = TRUE, seed = seed, tol = 1e-20, ncores = 0)
#write_imputed_data(imp$imputed_data, paste0(path, "/3/"))
