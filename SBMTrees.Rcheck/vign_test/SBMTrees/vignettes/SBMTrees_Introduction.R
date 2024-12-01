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


## ----prediction---------------------------------------------------------------
# Fit the predictive model
model <- BMTrees_prediction(
  X_train, Y_train, Z_train, subject_id_train, 
  X_test, Z_test, subject_id_test, 
  model = "BMTrees", 
  binary = FALSE, 
  nburn = 3000L, 
  npost = 4000L, 
  skip = 1L, 
  verbose = FALSE, 
  seed = 123
)

# Posterior expectation for the testing dataset
posterior_predictions <- model$post_predictive_y_test
head(colMeans(posterior_predictions))

## ----prediction_evaluation----------------------------------------------------
point_predictions = colMeans(posterior_predictions)

# Compute MAE
mae <- mean(abs(point_predictions - Y_test_true))
cat("Mean Absolute Error (MAE):", mae, "\n")

# Compute MSE
mse <- mean((point_predictions - Y_test_true)^2)
cat("Mean Squared Error (MSE):", mse, "\n")

# Compute 95% credible intervals
lower_bounds <- apply(posterior_predictions, 2, quantile, probs = 0.025)
upper_bounds <- apply(posterior_predictions, 2, quantile, probs = 0.975)

# Check if true values fall within the intervals
coverage <- mean(Y_test_true >= lower_bounds & Y_test_true <= upper_bounds)
cat("95% Posterior Predictive Interval Coverage:", coverage * 100, "%\n")



plot(Y_test_true, point_predictions, 
     xlab = "True Values", 
     ylab = "Predicted Values", 
     main = "True vs Predicted Values")
abline(0, 1, col = "red") # Add a 45-degree reference line


## ----imputation_sim-----------------------------------------------------------
# Simulate data with missing values
data <- simulation_imputation(
  n_subject = 200, 
  seed = 1234, 
  nonrandeff = TRUE, 
  nonresidual = TRUE, 
  alligned = FALSE
)

# Extract components of the dataset
X_mis <- data$X_mis   # Covariates with missing values
Y_mis <- data$Y_mis   # Outcomes with missing values
Z <- data$Z           # Random predictors
subject_id <- data$subject_id


