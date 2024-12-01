## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(SBMTrees)
library(mitml)
library(lme4)

## ----prediction_sim-----------------------------------------------------------
# Simulate data
data <- simulation_prediction(
  n_subject = 300, 
  seed = 123, 
  nonlinear = TRUE, 
  nonrandeff = TRUE, 
  nonresidual = TRUE
)

# Extract training and testing datasets
X_train <- data$X_train
Y_train <- data$Y_train
Z_train <- data$Z_train
subject_id_train <- data$subject_id_train

X_test <- data$X_test
Z_test <- data$Z_test
subject_id_test <- data$subject_id_test

Y_test_true <- data$Y_test_true


