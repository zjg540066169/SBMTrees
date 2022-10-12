library(tidyverse)
library(dbarts)
library(lme4)
library(mice)
source("BMETrees_sequential.R")
source("ribart_sequential.R")
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
set.seed(111)
df <- data_generation(sigma_1  = 1, tau_1 = 2, delta_1 = 1, sigma_2  = 1, tau_2 = 2, delta_2 = 1, n_total = 5000, number_of_time_points = 10)
head(df)
# Add missingness 
set.seed(111)
df_with_missing <- df %>% 
  mutate(
    x_1 = case_when(subject_id %in% sample(1:500, size = 188) & time <= 4 ~ NA_real_,
                    TRUE ~ x_1),
    x_2 = case_when(subject_id %in% sample(1:500, size = 250) & time >= 5 ~ NA_real_,
                    TRUE ~ x_2)) %>% 
  mutate(time = time/10)
# mutate(time = time/10)
head(df_with_missing)
sum(is.na(df_with_missing[,2]))/5000
sum(is.na(df_with_missing[,3]))/5000


# Random intercept BART--------------------------------------------------------
rbart_vi_sequential_result <-
  rbart_vi_sequential(
    data_x =  df_with_missing[1:3],
    subject_id = df_with_missing$subject_id,
    time = df_with_missing$time,
    number_of_MCMC = 1000
  )
save(rbart_vi_miss_imputation_result, file = "Rdata/ri_bart_imputation_burnin_100.Rdata")


# Random slope BART--------------------------------------------------------
rbart_vi_miss_imputation_result$x_1
nu_prior_tau = 2
lambda_prior_tau = 1
nu_prior_delta = 2
lambda_prior_delta = 1
initial_tau <- 1
initial_random_intercept <- 0
initial_delta <- 1
random_slope <- df_with_missing$time
BMETrees_sequential_result <-
  BMETrees_sequential(
    data_x =  df_with_missing[1:3],
    subject_id = df_with_missing$subject_id,
    time = df_with_missing$time,
    number_of_MCMC = 2
  )
# linear mixed effect model--------------------------------------------------------
head(df_with_missing)
sum(is.na(df_with_missing[,4]))
result_x_1_imputation <- matrix(NA, nrow = nrow(df_with_missing), ncol = 1000)
result_x_2_imputation <- matrix(NA, nrow = nrow(df_with_missing), ncol = 1000)
for (h in 1:1000){
  # set.seed(h + 111)
  # df <- data_generation(sigma_1  = 1, tau_1 = 2, delta_1 = 1, sigma_2  = 1, tau_2 = 2, delta_2 = 1, n_total = 5000, number_of_time_points = 10)
  # head(df)
  # # Add missingness 
  # set.seed(h + 111)
  df_with_missing <- df %>% 
    mutate(
      x_1 = case_when(subject_id %in% sample(1:500, size = 188) & time <= 4 ~ NA_real_,
                      TRUE ~ x_1),
      x_2 = case_when(subject_id %in% sample(1:500, size = 250) & time >= 5 ~ NA_real_,
                      TRUE ~ x_2)) %>% 
    mutate(time/10)
  lmer_1 <- lmer(x_1 ~ x_0 + (time|subject_id), data = df_with_missing)
  result_x_1_imputation[, h] <- predict(lmer_1, newdata = df_with_missing[,c("x_0", "time", "subject_id")]) 
  df_with_missing[,"x_1"] <- predict(lmer_1, newdata = df_with_missing[,c("x_0", "time", "subject_id")])
  lmer_2 <- lmer(x_2 ~ x_0 + x_1 + (time|subject_id), data = df_with_missing)
  result_x_2_imputation[, h] <- predict(lmer_2, newdata = df_with_missing[,c("x_0", "x_1","time", "subject_id")]) 
  # df_with_missing[,"x_2"] <-predict(lmer_2, newdata = df_with_missing[,c("x_0", "x_1","time", "subject_id")]) 
  print(h)
}
save(result_x_1_imputation, result_x_2_imputation, file = "Rdata/linear_mixed_effect_model_imputation.Rdata")
# multivariate linear mixed model -----------------------------------------

imp_df_with_missing <- mice( as.matrix(df_with_missing)  , maxit=0 )
predM_imp_df_with_missing <- imp_df_with_missing$predictorMatrix
impM_imp_df_with_missing <- imp_df_with_missing$method
predM1_df_with_missing <- predM_imp_df_with_missing
predM1_df_with_missing["x_1",] <- c(1,0,1,2,-2)
predM1_df_with_missing["x_2",] <- c(1,1,0,2,-2)

impM1_df_with_missing <- impM_imp_df_with_missing
impM1_df_with_missing <- c("","2l.pan","2l.pan","","")

imp1 <- mice( as.matrix( df_with_missing ) , predictorMatrix = predM1_df_with_missing , 
              method = impM1_df_with_missing , maxit=1,m = 200)
save(imp1, file = "Rdata/multivariate_linear_mixed_model_imputation.Rdata")
complete(imp1) %>% dim
imp1$imp$x_1 %>% dim
df$x_1[as.numeric(rownames(imp1$imp$x_1))]
# (df$x_1[as.numeric(rownames(imp1$imp$x_1))]- apply(imp1$imp$x_1,1,mean)) %>% mean

apply(imp1$imp$x_1,1,sd) 
imp1$imp$x_2 %>% dim


