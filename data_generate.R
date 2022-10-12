library(tidyverse)
#source("BMETrees_sequential.R")
#source("ribart_sequential.R")
# sigma_1 = 1
# tau_1 = 2
# delta_1 = 1
# sigma_2 = 1
# tau_2 = 2
# delta_2 = 1
# n_total = 500
# number_of_time_points = 10
data_generation <- function(sigma_1, tau_1, delta_1, sigma_2, tau_2, delta_2, n_total, number_of_time_points){
  n_unique_subject <- n_total/number_of_time_points
  subject_id <- rep(1:n_unique_subject, number_of_time_points)
  x_0 <- runif(n_total)
  time <- rep(1:number_of_time_points,each  = n_unique_subject)
  
  # Generate x_1 given x_0
  x_1_given_x_0_mean <- sin(x_0) + x_0 ^ 2 + x_0
  x_1_given_x_0_tree  <- rnorm(n_total, x_1_given_x_0_mean, sigma_1)
  alpha_1 <- rnorm(n_unique_subject, 0, tau_1)
  beta_1 <- rnorm(n_unique_subject, 0, delta_1)
  x_1_given_x_0 <- x_1_given_x_0_tree + alpha_1[subject_id] + time* beta_1[subject_id]
  
  # Generate x_2 given x_0 and x_1
  x_2_given_x_1_0_mean <- sin(x_1_given_x_0 * x_0) + x_1_given_x_0 ^ 1/2 + x_0
  x_2_given_x_1_0_tree  <- rnorm(n_total, x_2_given_x_1_0_mean, sigma_2)
  alpha_2 <- rnorm(n_unique_subject, 0, tau_2)
  beta_2 <- rnorm(n_unique_subject, 0, delta_2)
  x_2_given_x_1_0 <- x_2_given_x_1_0_tree + alpha_2[subject_id] + time* beta_2[subject_id]
  
  df <- data.frame(x_0 = x_0, x_1 = x_1_given_x_0, x_2 = x_2_given_x_1_0) %>% cbind(time) %>% cbind(subject_id)
  return(df)
}