library(Rcpp)
library(RcppParallel)
library(lme4)
library(optimx)
library(tidyverse)

source("R/sequential_imputation_R.R")
source("R/simulate.R")
source("R/evaluation.R")
#sourceCpp("./src/BMTrees_MCMC.cpp")

# model
# 1 nDP residual
# 2 cDP random effects
# 3 double DP



write_imputed_data = function(imputed_X, dic){
  if(!dir.exists(dic)){
    dir.create(dic)
  }
  for (i in 1:dim(imputed_X)[1]) {
    document_path = paste(dic, "/", i, ".csv", sep = "")
    write_csv(as.data.frame(imputed_X[i,,]), document_path)
  }
}

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
seed = as.integer(ags[1])
print(paste0("seed:", seed))
path_ = "./simulation_new/"
i = as.integer(ags)[1]
dir.create(paste0(path_, i))
scenario = paste0(path_, i, "/")

path = paste0(scenario, folder_name(F, T, T))
dir.create(path)
data = simulation(1000, n_obs_per_sub = 5, nonlinear = F, nonrandeff = T, nonresidual = T, alligned = F, seed = 1)
X = data$X_mis
subject_id = data$subject_id
Y = data$Y
Z = data$Z

write_csv(data$X_O, paste0(path, "/X_O.csv"))
write_csv(data$X_mis, paste0(path, "/X_mis.csv"))
write_csv(as.data.frame(data$Y), paste0(path, "/Y.csv"))
write_csv(as.data.frame(data$Y_O), paste0(path, "/Y_O.csv"))

write_csv(as.data.frame(data$Z), paste0(path, "/Z.csv"))
write_csv(as.data.frame(data$subject_id), paste0(path, "/subject_id.csv"))
write_csv(as.data.frame(data$trajectory), paste0(path, "/trajectory.csv"))


type = c(1, 1, 0, 0, 0, 0, 0)
nburn = 3000
npost = 4000
skip = 200
imp = sequential_imputation(X, Y, type, Z, subject_id, model = 2, nburn = nburn, npost = npost, skip = skip, verbose = TRUE, seed = seed, tol = 1e-20, ncores = 0)
write_imputed_data(imp$imputed_data, paste0(path, "/2/"))

imp = sequential_imputation_R(X, Y, type, Z, subject_id, model = 0, nburn = nburn, npost = npost, n_model_burn = 0L, skip = skip, verbose = TRUE, seed = seed, tol = 1e-20, ncores = 0)
write_imputed_data(imp$imputed_data, paste0(path, "/0/"))

imp = sequential_imputation_R(X, Y, type, Z, subject_id, model = 1, nburn = nburn, npost = npost, n_model_burn = 0L, skip = skip, verbose = TRUE, seed = seed, tol = 1e-20, ncores = 0)
write_imputed_data(imp$imputed_data, paste0(path, "/1/"))
