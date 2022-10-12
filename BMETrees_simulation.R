library(dbarts)
library(tidyverse)
# Data generation function
data_generation <- function(sigma, tau, delta, n_total, number_of_time_points){
  n_unique_subject <- n_total/number_of_time_points
  x  <- matrix(runif(n_total * 5), n_total, 5)
  time <- rep(1:number_of_time_points,each  = n_unique_subject)/10
  Ey <- 10 * sin(pi * x[,1] * x[,2]) + 20 * (x[,3] - 0.5)^2 + 10 * x[,4] + 5 * x[,5] 
  y_tree  <- rnorm(n_total, Ey, sigma)
  subject_id <- rep(1:n_unique_subject, number_of_time_points)
  alpha <- rnorm(n_unique_subject, 0, tau)
  beta<- rnorm(n_unique_subject, 0, delta)
  y <- y_tree + alpha[subject_id] + time* beta[subject_id]
  df <- as.data.frame(x[,1:5]) %>% cbind(time)
  colnames(df)[1:5] <- paste0("x_", 1:5)
  df$y <- y
  df$subject_id <- subject_id
  return(df)
}

bart_random_intercept_tau <- rep(NA, 200)
bart_random_intercept_sigma <- rep(NA, 200)
bart_random_intercept_rmse <- rep(NA,200)
bart_random_intercept_bias <- rep(NA,200)
bart_random_intercept_bias_time_1 <- rep(NA,200)
bart_random_intercept_bias_time_2 <- rep(NA,200)
bart_random_intercept_bias_time_3 <- rep(NA,200)
bart_random_intercept_bias_time_4 <- rep(NA,200)
bart_random_intercept_bias_time_5 <- rep(NA,200)
bart_random_intercept_bias_time_6 <- rep(NA,200)
bart_random_intercept_bias_time_7 <- rep(NA,200)
bart_random_intercept_bias_time_8 <- rep(NA,200)
bart_random_intercept_bias_time_9 <- rep(NA,200)
bart_random_intercept_bias_time_10 <- rep(NA,200)

bart_random_intercept_rmse_time_1 <- rep(NA,200)
bart_random_intercept_rmse_time_2 <- rep(NA,200)
bart_random_intercept_rmse_time_3 <- rep(NA,200)
bart_random_intercept_rmse_time_4 <- rep(NA,200)
bart_random_intercept_rmse_time_5 <- rep(NA,200)
bart_random_intercept_rmse_time_6 <- rep(NA,200)
bart_random_intercept_rmse_time_7 <- rep(NA,200)
bart_random_intercept_rmse_time_8 <- rep(NA,200)
bart_random_intercept_rmse_time_9 <- rep(NA,200)
bart_random_intercept_rmse_time_10 <- rep(NA,200)

bart_random_slope_tau <- rep(NA, 200)
bart_random_slope_delta <- rep(NA, 200)
bart_random_slope_sigma <- rep(NA, 200)
bart_random_slope_rmse <- rep(NA,200)
bart_random_slope_bias <- rep(NA,200)
bart_random_intercept_bias_time_1_once <- rep(NA,500)
bart_random_intercept_bias_time_2_once <- rep(NA,500)
bart_random_intercept_bias_time_3_once <- rep(NA,500)
bart_random_intercept_bias_time_4_once <- rep(NA,500)
bart_random_intercept_bias_time_5_once <- rep(NA,500)
bart_random_intercept_bias_time_6_once <- rep(NA,500)
bart_random_intercept_bias_time_7_once <- rep(NA,500)
bart_random_intercept_bias_time_8_once <- rep(NA,500)
bart_random_intercept_bias_time_9_once <- rep(NA,500)
bart_random_intercept_bias_time_10_once <- rep(NA,500)


# Simulation --------------------------------------------------------------

for (i in 1:200){
  set.seed(i+999)
  df <- data_generation(sigma = 1, tau = 1.5, delta = 2.5, n_total = 5000, number_of_time_points = 10)
  ## Random intercept: 
  bart_random_intercept <- rbart_vi(y ~ . - subject_id, df, group.by = subject_id,
                                        n.samples = 5000L, n.burn = 1000L, n.chains = 1L, n.threads = 1L)
  bart_random_intercept_tau[i] <- mean(bart_random_intercept$tau)
  bart_random_intercept_sigma[i] <- mean(bart_random_intercept$sigma)
  bart_random_intercept_rmse[i] <- sqrt(mean((
    df$y - apply(bart_random_intercept$yhat.train, 2, mean)
  ) ^ 2))
  bart_random_intercept_bias[i] <- mean(df$y - apply(bart_random_intercept$yhat.train,2,mean))
  bart_random_intercept_rmse_time_1[i] <- sqrt(mean((
    df$y[(1+500*(1-1)):(500+500*(1-1))] - apply(bart_random_intercept$yhat.train, 2, mean)[(1+500*(1-1)):(500+500*(1-1))]
  ) ^ 2))
  bart_random_intercept_bias_time_1[i] <- mean(df$y[(1+500*(1-1)):(500+500*(1-1))] - apply(bart_random_intercept$yhat.train,2,mean)[(1+500*(1-1)):(500+500*(1-1))])
  
  bart_random_intercept_rmse_time_2[i] <- sqrt(mean((
    df$y[(1+500*(2-1)):(500+500*(2-1))] - apply(bart_random_intercept$yhat.train, 2, mean)[(1+500*(2-1)):(500+500*(2-1))]
  ) ^ 2))
  bart_random_intercept_bias_time_2[i] <- mean(df$y[(1+500*(2-1)):(500+500*(2-1))] - apply(bart_random_intercept$yhat.train,2,mean)[(1+500*(2-1)):(500+500*(2-1))])
  
  
  bart_random_intercept_rmse_time_3[i] <- sqrt(mean((
    df$y[(1+500*(3-1)):(500+500*(3-1))] - apply(bart_random_intercept$yhat.train, 2, mean)[(1+500*(3-1)):(500+500*(3-1))]
  ) ^ 2))
  bart_random_intercept_bias_time_3[i] <- mean(df$y[(1+500*(3-1)):(500+500*(3-1))] - apply(bart_random_intercept$yhat.train,2,mean)[(1+500*(3-1)):(500+500*(3-1))])
  
  bart_random_intercept_rmse_time_4[i] <- sqrt(mean((
    df$y[(1+500*(4-1)):(500+500*(4-1))] - apply(bart_random_intercept$yhat.train, 2, mean)[(1+500*(4-1)):(500+500*(4-1))]
  ) ^ 2))
  bart_random_intercept_bias_time_4[i] <- mean(df$y[(1+500*(4-1)):(500+500*(4-1))] - apply(bart_random_intercept$yhat.train,2,mean)[(1+500*(4-1)):(500+500*(4-1))])
  
  bart_random_intercept_rmse_time_5[i] <- sqrt(mean((
    df$y[(1+500*(5-1)):(500+500*(5-1))] - apply(bart_random_intercept$yhat.train, 2, mean)[(1+500*(5-1)):(500+500*(5-1))]
  ) ^ 2))
  bart_random_intercept_bias_time_5[i] <- mean(df$y[(1+500*(5-1)):(500+500*(5-1))] - apply(bart_random_intercept$yhat.train,2,mean)[(1+500*(5-1)):(500+500*(5-1))])
  
  bart_random_intercept_rmse_time_6[i] <- sqrt(mean((
    df$y[(1+500*(6-1)):(500+500*(6-1))] - apply(bart_random_intercept$yhat.train, 2, mean)[(1+500*(6-1)):(500+500*(6-1))]
  ) ^ 2))
  bart_random_intercept_bias_time_6[i] <- mean(df$y[(1+500*(6-1)):(500+500*(6-1))] - apply(bart_random_intercept$yhat.train,2,mean)[(1+500*(6-1)):(500+500*(6-1))])
  
  bart_random_intercept_rmse_time_7[i] <- sqrt(mean((
    df$y[(1+500*(7-1)):(500+500*(7-1))] - apply(bart_random_intercept$yhat.train, 2, mean)[(1+500*(7-1)):(500+500*(7-1))]
  ) ^ 2))
  bart_random_intercept_bias_time_7[i] <- mean(df$y[(1+500*(7-1)):(500+500*(7-1))] - apply(bart_random_intercept$yhat.train,2,mean)[(1+500*(7-1)):(500+500*(7-1))])
  
  bart_random_intercept_rmse_time_8[i] <- sqrt(mean((
    df$y[(1+500*(8-1)):(500+500*(8-1))] - apply(bart_random_intercept$yhat.train, 2, mean)[(1+500*(8-1)):(500+500*(8-1))]
  ) ^ 2))
  bart_random_intercept_bias_time_8[i] <- mean(df$y[(1+500*(8-1)):(500+500*(8-1))] - apply(bart_random_intercept$yhat.train,2,mean)[(1+500*(8-1)):(500+500*(8-1))])
  
  bart_random_intercept_rmse_time_9[i] <- sqrt(mean((
    df$y[(1+500*(9-1)):(500+500*(9-1))] - apply(bart_random_intercept$yhat.train, 2, mean)[(1+500*(9-1)):(500+500*(9-1))]
  ) ^ 2))
  bart_random_intercept_bias_time_9[i] <- mean(df$y[(1+500*(9-1)):(500+500*(9-1))] - apply(bart_random_intercept$yhat.train,2,mean)[(1+500*(9-1)):(500+500*(9-1))])
  
  bart_random_intercept_rmse_time_10[i] <- sqrt(mean((
    df$y[(1+500*(10-1)):(500+500*(10-1))] - apply(bart_random_intercept$yhat.train, 2, mean)[(1+500*(10-1)):(500+500*(10-1))]
  ) ^ 2))
  bart_random_intercept_bias_time_10[i] <- mean(df$y[(1+500*(10-1)):(500+500*(10-1))] - apply(bart_random_intercept$yhat.train,2,mean)[(1+500*(10-1)):(500+500*(10-1))])
  bart_random_intercept_bias_time_1_once <- df$y[(1+500*(1-1)):(500+500*(1-1))] - apply(bart_random_intercept$yhat.train,2,mean)[(1+500*(1-1)):(500+500*(1-1))]
  bart_random_intercept_bias_time_2_once <- df$y[(1+500*(2-1)):(500+500*(2-1))] - apply(bart_random_intercept$yhat.train,2,mean)[(1+500*(2-1)):(500+500*(2-1))]
  bart_random_intercept_bias_time_3_once <- df$y[(1+500*(3-1)):(500+500*(3-1))] - apply(bart_random_intercept$yhat.train,2,mean)[(1+500*(3-1)):(500+500*(3-1))]
  bart_random_intercept_bias_time_4_once <- df$y[(1+500*(4-1)):(500+500*(4-1))] - apply(bart_random_intercept$yhat.train,2,mean)[(1+500*(4-1)):(500+500*(4-1))]
  bart_random_intercept_bias_time_5_once <- df$y[(1+500*(5-1)):(500+500*(5-1))] - apply(bart_random_intercept$yhat.train,2,mean)[(1+500*(5-1)):(500+500*(5-1))]
  bart_random_intercept_bias_time_6_once <- df$y[(1+500*(6-1)):(500+500*(6-1))] - apply(bart_random_intercept$yhat.train,2,mean)[(1+500*(6-1)):(500+500*(6-1))]
  bart_random_intercept_bias_time_7_once <- df$y[(1+500*(7-1)):(500+500*(7-1))] - apply(bart_random_intercept$yhat.train,2,mean)[(1+500*(7-1)):(500+500*(7-1))]
  bart_random_intercept_bias_time_8_once <- df$y[(1+500*(8-1)):(500+500*(8-1))] - apply(bart_random_intercept$yhat.train,2,mean)[(1+500*(8-1)):(500+500*(8-1))]
  bart_random_intercept_bias_time_9_once <- df$y[(1+500*(9-1)):(500+500*(9-1))] - apply(bart_random_intercept$yhat.train,2,mean)[(1+500*(9-1)):(500+500*(9-1))]
  bart_random_intercept_bias_time_10_once <- df$y[(1+500*(10-1)):(500+500*(10-1))] - apply(bart_random_intercept$yhat.train,2,mean)[(1+500*(10-1)):(500+500*(10-1))]
  
  # Random slope:
  random_slope <- df$time
  lambda_prior_delta <- 1
  nu_prior_tau = 2
  lambda_prior_tau = 1
  nu_prior_delta = 2
  lambda_prior_delta = 1
  initial_tau <- 1
  initial_delta <- 1
  initial_sigma <- df %>%
    group_by(subject_id) %>%
    nest() %>%
    mutate(model = map(data,~lm(y ~ ., data = .))) %>%
    mutate(model_info = map(model, broom::augment))%>%
    dplyr::select(-data, -model) %>%
    unnest(cols = c(model_info)) %>%
    ungroup %>%
    pull(.resid) %>%
    sd

  initial_random_intercept <- df %>%
    group_by(subject_id) %>%
    nest() %>%
    mutate(model = map(data,~lm(y ~ ., data = .))) %>%
    mutate(model_info = map(model, broom::augment))%>%
    dplyr::select(-data, -model) %>%
    unnest(cols = c(model_info)) %>%
    ungroup %>%
    pull(.resid) %>%
    mean

  initial_tau <- df %>%
    group_by(subject_id) %>%
    nest() %>%
    mutate(model = map(data,~lm(y ~ ., data = .))) %>%
    mutate(model_info = map(model, broom::tidy)) %>%
    dplyr::select(-data, -model) %>%
    unnest(cols = c(model_info)) %>%
    ungroup %>%
    filter(term == "(Intercept)") %>%
    pull(estimate) %>%
    sd

  initial_delta <- df %>%
    group_by(subject_id) %>%
    nest() %>%
    mutate(model = map(data,~lm(y ~ ., data = .))) %>%
    mutate(model_info = map(model, broom::tidy)) %>%
    dplyr::select(-data, -model) %>%
    unnest(cols = c(model_info)) %>%
    ungroup %>%
    filter(term == "x_5") %>%
    pull(estimate) %>%
    sd
  bart_random_slope <- BMETrees(y ~ . - subject_id, df, group.by = df$subject_id,prior = gamma,
                                           random_slope = x_5,
                                           nu_prior_tau = 3,
                                           lambda_prior_tau = 1,
                                           nu_prior_delta = 2,
                                           lambda_prior_delta = 1,
                                           initial_sigma = initial_sigma,
                                           initial_random_intercept = initial_random_intercept,
                                           initial_tau = initial_tau,
                                           initial_delta  = initial_delta,
                                           n.samples = 5000L, n.burn = 1000L, n.chains = 1L, n.threads = 1L)
  bart_random_slope_tau[i] <- mean(bart_random_slope$tau[1001:5000])
  bart_random_slope_sigma[i] <- mean(bart_random_slope$sigma[1001:5000])
  bart_random_slope_delta[i] <- mean(bart_random_slope$delta[1001:5000])
  bart_random_slope_rmse[i] <- sqrt(mean((
    df$y - apply(bart_random_slope$yhat.train, 2, mean)
  ) ^ 2))
  bart_random_slope_bias[i] <- mean(df$y - apply(bart_random_slope$yhat.train,2,mean))
  bart_random_slope_bias[i] <- mean(df$y - apply(bart_random_slope$yhat.train,2,mean))
  bart_random_slope_rmse_time_1[i] <- sqrt(mean((
    df$y[(1+500*(1-1)):(500+500*(1-1))] - apply(bart_random_slope$yhat.train, 2, mean)[(1+500*(1-1)):(500+500*(1-1))]
  ) ^ 2))
  bart_random_slope_bias_time_1[i] <- mean(df$y[(1+500*(1-1)):(500+500*(1-1))] - apply(bart_random_slope$yhat.train,2,mean)[(1+500*(1-1)):(500+500*(1-1))])
  
  bart_random_slope_rmse_time_2[i] <- sqrt(mean((
    df$y[(1+500*(2-1)):(500+500*(2-1))] - apply(bart_random_slope$yhat.train, 2, mean)[(1+500*(2-1)):(500+500*(2-1))]
  ) ^ 2))
  bart_random_slope_bias_time_2[i] <- mean(df$y[(1+500*(2-1)):(500+500*(2-1))] - apply(bart_random_slope$yhat.train,2,mean)[(1+500*(2-1)):(500+500*(2-1))])
  
  
  bart_random_slope_rmse_time_3[i] <- sqrt(mean((
    df$y[(1+500*(3-1)):(500+500*(3-1))] - apply(bart_random_slope$yhat.train, 2, mean)[(1+500*(3-1)):(500+500*(3-1))]
  ) ^ 2))
  bart_random_slope_bias_time_3[i] <- mean(df$y[(1+500*(3-1)):(500+500*(3-1))] - apply(bart_random_slope$yhat.train,2,mean)[(1+500*(3-1)):(500+500*(3-1))])
  
  bart_random_slope_rmse_time_4[i] <- sqrt(mean((
    df$y[(1+500*(4-1)):(500+500*(4-1))] - apply(bart_random_slope$yhat.train, 2, mean)[(1+500*(4-1)):(500+500*(4-1))]
  ) ^ 2))
  bart_random_slope_bias_time_4[i] <- mean(df$y[(1+500*(4-1)):(500+500*(4-1))] - apply(bart_random_slope$yhat.train,2,mean)[(1+500*(4-1)):(500+500*(4-1))])
  
  bart_random_slope_rmse_time_5[i] <- sqrt(mean((
    df$y[(1+500*(5-1)):(500+500*(5-1))] - apply(bart_random_slope$yhat.train, 2, mean)[(1+500*(5-1)):(500+500*(5-1))]
  ) ^ 2))
  bart_random_slope_bias_time_5[i] <- mean(df$y[(1+500*(5-1)):(500+500*(5-1))] - apply(bart_random_slope$yhat.train,2,mean)[(1+500*(5-1)):(500+500*(5-1))])
  
  bart_random_slope_rmse_time_6[i] <- sqrt(mean((
    df$y[(1+500*(6-1)):(500+500*(6-1))] - apply(bart_random_slope$yhat.train, 2, mean)[(1+500*(6-1)):(500+500*(6-1))]
  ) ^ 2))
  bart_random_slope_bias_time_6[i] <- mean(df$y[(1+500*(6-1)):(500+500*(6-1))] - apply(bart_random_slope$yhat.train,2,mean)[(1+500*(6-1)):(500+500*(6-1))])
  
  bart_random_slope_rmse_time_7[i] <- sqrt(mean((
    df$y[(1+500*(7-1)):(500+500*(7-1))] - apply(bart_random_slope$yhat.train, 2, mean)[(1+500*(7-1)):(500+500*(7-1))]
  ) ^ 2))
  bart_random_slope_bias_time_7[i] <- mean(df$y[(1+500*(7-1)):(500+500*(7-1))] - apply(bart_random_slope$yhat.train,2,mean)[(1+500*(7-1)):(500+500*(7-1))])
  
  bart_random_slope_rmse_time_8[i] <- sqrt(mean((
    df$y[(1+500*(8-1)):(500+500*(8-1))] - apply(bart_random_slope$yhat.train, 2, mean)[(1+500*(8-1)):(500+500*(8-1))]
  ) ^ 2))
  bart_random_slope_bias_time_8[i] <- mean(df$y[(1+500*(8-1)):(500+500*(8-1))] - apply(bart_random_slope$yhat.train,2,mean)[(1+500*(8-1)):(500+500*(8-1))])
  
  bart_random_slope_rmse_time_9[i] <- sqrt(mean((
    df$y[(1+500*(9-1)):(500+500*(9-1))] - apply(bart_random_slope$yhat.train, 2, mean)[(1+500*(9-1)):(500+500*(9-1))]
  ) ^ 2))
  bart_random_slope_bias_time_9[i] <- mean(df$y[(1+500*(9-1)):(500+500*(9-1))] - apply(bart_random_slope$yhat.train,2,mean)[(1+500*(9-1)):(500+500*(9-1))])
  
  bart_random_slope_rmse_time_10[i] <- sqrt(mean((
    df$y[(1+500*(10-1)):(500+500*(10-1))] - apply(bart_random_slope$yhat.train, 2, mean)[(1+500*(10-1)):(500+500*(10-1))]
  ) ^ 2))
  bart_random_slope_bias_time_10[i] <- mean(df$y[(1+500*(10-1)):(500+500*(10-1))] - apply(bart_random_slope$yhat.train,2,mean)[(1+500*(10-1)):(500+500*(10-1))])
  bart_random_slope_bias_time_1_once <- df$y[(1+500*(1-1)):(500+500*(1-1))] - apply(bart_random_slope$yhat.train,2,mean)[(1+500*(1-1)):(500+500*(1-1))]
  bart_random_slope_bias_time_2_once <- df$y[(1+500*(2-1)):(500+500*(2-1))] - apply(bart_random_slope$yhat.train,2,mean)[(1+500*(2-1)):(500+500*(2-1))]
  bart_random_slope_bias_time_3_once <- df$y[(1+500*(3-1)):(500+500*(3-1))] - apply(bart_random_slope$yhat.train,2,mean)[(1+500*(3-1)):(500+500*(3-1))]
  bart_random_slope_bias_time_4_once <- df$y[(1+500*(4-1)):(500+500*(4-1))] - apply(bart_random_slope$yhat.train,2,mean)[(1+500*(4-1)):(500+500*(4-1))]
  bart_random_slope_bias_time_5_once <- df$y[(1+500*(5-1)):(500+500*(5-1))] - apply(bart_random_slope$yhat.train,2,mean)[(1+500*(5-1)):(500+500*(5-1))]
  bart_random_slope_bias_time_6_once <- df$y[(1+500*(6-1)):(500+500*(6-1))] - apply(bart_random_slope$yhat.train,2,mean)[(1+500*(6-1)):(500+500*(6-1))]
  bart_random_slope_bias_time_7_once <- df$y[(1+500*(7-1)):(500+500*(7-1))] - apply(bart_random_slope$yhat.train,2,mean)[(1+500*(7-1)):(500+500*(7-1))]
  bart_random_slope_bias_time_8_once <- df$y[(1+500*(8-1)):(500+500*(8-1))] - apply(bart_random_slope$yhat.train,2,mean)[(1+500*(8-1)):(500+500*(8-1))]
  bart_random_slope_bias_time_9_once <- df$y[(1+500*(9-1)):(500+500*(9-1))] - apply(bart_random_slope$yhat.train,2,mean)[(1+500*(9-1)):(500+500*(9-1))]
  bart_random_slope_bias_time_10_once <- df$y[(1+500*(10-1)):(500+500*(10-1))] - apply(bart_random_slope$yhat.train,2,mean)[(1+500*(10-1)):(500+500*(10-1))]
  print(i)
}


