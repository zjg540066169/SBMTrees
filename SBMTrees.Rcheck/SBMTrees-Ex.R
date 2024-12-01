pkgname <- "SBMTrees"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
base::assign(".ExTimings", "SBMTrees-Ex.timings", pos = 'CheckExEnv')
base::cat("name\tuser\tsystem\telapsed\n", file=base::get(".ExTimings", pos = 'CheckExEnv'))
base::assign(".format_ptime",
function(x) {
  if(!is.na(x[4L])) x[1L] <- x[1L] + x[4L]
  if(!is.na(x[5L])) x[2L] <- x[2L] + x[5L]
  options(OutDec = '.')
  format(x[1L:3L], digits = 7L)
},
pos = 'CheckExEnv')

### * </HEADER>
library('SBMTrees')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("BMTrees_prediction")
### * BMTrees_prediction

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: BMTrees_prediction
### Title: Bayesian Trees Mixed-Effects Models for Predicting Longitudinal
###   Outcomes
### Aliases: BMTrees_prediction

### ** Examples

## Not run: 
##D data = simulation_prediction(n_subject = 800, seed = 123, nonlinear = TRUE, 
##D nonrandeff = TRUE, nonresidual = TRUE) 
##D model = BMTrees_prediction(data$X_train, data$Y_train, data$Z_train, 
##D data$subject_id_train, data$X_test, data$Z_test, data$subject_id_test, model = "BMTrees", 
##D binary = FALSE, nburn = 3000L, npost = 4000L, skip = 1L, verbose = TRUE, seed = 123)
##D model$post_predictive_y_test
##D model$post_sigma
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("BMTrees_prediction", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("SBMTrees-package")
### * SBMTrees-package

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: SBMTrees-package
### Title: Sequential Imputation with Bayesian Trees Mixed-Effects Models
### Aliases: SBMTrees-package SBMTrees
### Keywords: SBMTrees, longitudinal missing data, Bayesian non-parametric
###   methods, sequential imputation, mixed-effects models, multiple
###   imputation

### ** Examples

  ## Not run: 
##D      ## Example of predicting longitudinal outcomes
##D      data <- simulation_prediction(n_subject = 800, seed = 123, nonlinear = TRUE,
##D              nonrandeff = TRUE, nonresidual = TRUE)
##D      X_train <- data$X_train
##D      Y_train <- data$Y_train
##D      Z_train <- data$Z_train
##D      subject_id_train <- data$subject_id_train
##D 
##D      X_test <- data$X_test
##D      Z_test <- data$Z_test
##D      subject_id_test <- data$subject_id_test
##D 
##D      model <- BMTrees_prediction(X_train, Y_train, Z_train, subject_id_train, 
##D               X_test, Z_test, subject_id_test, model = "BMTrees")
##D      model$post_predictive_y_test
##D      
##D      
##D      
##D      data2 = simulation_imputation(n_subject = 800, seed = 123, nonrandeff = TRUE, 
##D              nonresidual = TRUE, alligned = F) 
##D      X_mis = data2$X_mis # get missing covariates
##D      Y_mis = data2$Y_mis # get missing outcomes
##D      Z = data2$Z # get random predictors
##D      subject_id = data2$subject_id  # get subject id
##D 
##D      model2 = sequential_imputation(X_mis, Y_mis, Z, subject_id, rep(0, 9), F, 
##D              model = "BMTrees", nburn = 3000L, npost = 4000L, skip = 200L, 
##D              verbose = TRUE, seed = 123)
##D      model2$imputed_data
##D      model2$imputed_data[,,10]
##D   
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("SBMTrees-package", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("apply_locf_nocb")
### * apply_locf_nocb

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: apply_locf_nocb
### Title: Impute Missing Values Using LOCF and NOCB
### Aliases: apply_locf_nocb

### ** Examples

X <- matrix(c(NA, 2, NA, 4, 5, NA, 7, 8, NA, NA), nrow = 5, byrow = TRUE)
subject_id <- c(1, 1, 1, 2, 2)
apply_locf_nocb(X, subject_id)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("apply_locf_nocb", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("sequential_imputation")
### * sequential_imputation

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: sequential_imputation
### Title: Sequential Imputation for Missing Data
### Aliases: sequential_imputation

### ** Examples

## Not run: 
##D data <- simulation_imputation(n_subject = 800, seed = 123, nonrandeff = TRUE, 
##D         nonresidual = TRUE, alligned = FALSE) 
##D model <- sequential_imputation(data$X_mis, data$Y_mis, data$Z, data$subject_id, 
##D         rep(0, 9), binary_outcome = FALSE, model = "BMTrees", nburn = 3000L, 
##D         npost = 4000L, skip = 200L, verbose = TRUE, seed = 123)
##D model$imputed_data
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("sequential_imputation", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("simulation_imputation")
### * simulation_imputation

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: simulation_imputation
### Title: Simulate Longitudinal Data with Missingness
### Aliases: simulation_imputation

### ** Examples

simulated_data <- simulation_imputation(
  n_subject = 800,
  seed = 123,
  nonrandeff = TRUE,
  nonresidual = TRUE,
  alligned = FALSE
)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("simulation_imputation", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("simulation_prediction")
### * simulation_prediction

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: simulation_prediction
### Title: Simulate Longitudinal Data for Prediction
### Aliases: simulation_prediction

### ** Examples

  # Generate data with nonlinear associations and non-normal random effects and residuals
  data <- simulation_prediction(
    n_subject = 800,
    seed = 123,
    nonlinear = TRUE,
    nonrandeff = TRUE,
    nonresidual = TRUE
  )
  # Access training and testing data
  X_train <- data$X_train
  Y_train <- data$Y_train
  Z_train <- data$Z_train
  subject_id_train <- data$subject_id_train

  X_test <- data$X_test
  Y_test <- data$Y_test
  Z_test <- data$Z_test
  subject_id_test <- data$subject_id_test
  
  Y_test_true = data$Y_test_true




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("simulation_prediction", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
