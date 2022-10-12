library(dbarts)
# data_x <- df_with_missing[1:3]
# subject_id <- df_with_missing$subject_id
# time <- df_with_missing$time
# number_of_MCMC <- 50

rbart_vi_sequential <- function(data_x, subject_id, time, number_of_MCMC){
  df <- data_x
  subject <- subject_id
  time <- time
  ncol_df <- ncol(df)
  nrow_df <- nrow(df)
  # Matrix to store the results
  result_x_1_imputation <- matrix(NA, nrow = nrow(df), ncol = number_of_MCMC)
  result_x_2_imputation <- matrix(NA, nrow = nrow(df), ncol = number_of_MCMC)
  # Calculate the number of NA in each column
  sum_na <- rep(NA, ncol_df)
  for (i in 1:ncol_df) {
    sum_na[i] <- sum(is.na(df[, i]))
  }  
  sum_na
  # Reorder the colum from largest to smallest
  n_variable_need_impute <- ncol_df - length(sum_na[sum_na == 0])
  df_reorder <- df[, order(sum_na)]
  # Rename the columns x_0 to x_k
  colnames(df_reorder) <- paste0("x_", 0:(ncol_df-1))
  head(df_reorder)
  # sample means as the initial values for each column with missingness
  # imputation by column mean
  for (j in 1:n_variable_need_impute){
    df_reorder[, ncol_df + 1 - j][is.na(df_reorder[, ncol_df + 1 - j])] <- mean(df_reorder[,ncol_df + 1 - j], na.rm = T)
  }
  # Initial values for random intercept BART
  df_impute_1 <- cbind(df_reorder[1:2], subject_id = subject, time = time)
  lambda_prior_delta_1 <- 1
  nu_prior_tau = 2
  lambda_prior_tau_1 = 1
  nu_prior_delta = 2
  lambda_prior_delta = 1
  lambda_prior_tau <- 1
  initial_tau <- 1
  initial_delta <- 1
  random_slope <- time
  initial_sigma_1 <- df_impute_1 %>%
    group_by(subject_id) %>% 
    nest() %>%
    mutate(model = map(data,~lm(x_1 ~ ., data = .))) %>% 
    mutate(model_info = map(model, broom::augment))%>% 
    dplyr::select(-data, -model) %>%
    unnest(cols = c(model_info)) %>% 
    ungroup %>% 
    pull(.resid) %>% 
    sd
  
  initial_random_intercept_1 <- df_impute_1 %>%
    group_by(subject_id) %>% 
    nest() %>%
    mutate(model = map(data,~lm(x_1 ~ ., data = .))) %>% 
    mutate(model_info = map(model, broom::augment))%>% 
    dplyr::select(-data, -model) %>%
    unnest(cols = c(model_info)) %>% 
    ungroup %>% 
    pull(.resid) %>% 
    mean
  
  initial_tau_1 <- df_impute_1 %>%
    group_by(subject_id) %>%
    nest() %>%
    mutate(model = map(data,~lm(x_1 ~ ., data = .))) %>%
    mutate(model_info = map(model, broom::tidy)) %>%
    dplyr::select(-data, -model) %>%
    unnest(cols = c(model_info)) %>%
    ungroup %>%
    filter(term == "(Intercept)") %>%
    pull(estimate) %>%
    sd
  
  initial_delta_1 <- df_impute_1 %>%
    group_by(subject_id) %>%
    nest() %>%
    mutate(model = map(data,~lm(x_1 ~ ., data = .))) %>%
    mutate(model_info = map(model, broom::tidy)) %>%
    dplyr::select(-data, -model) %>%
    unnest(cols = c(model_info)) %>%
    ungroup %>%
    filter(term == "time") %>%
    pull(estimate) %>%
    sd
  
  df_impute_2 <- cbind(df_reorder[1:3], subject_id = subject, time = time)
  lambda_prior_delta_2 <- 1
  nu_prior_tau_2 = 2
  lambda_prior_tau_2 = 1
  nu_prior_delta_2 = 2
  lambda_prior_delta_2 = 1
  initial_sigma_2 <- df_impute_2 %>%
    group_by(subject_id) %>% 
    nest() %>%
    mutate(model = map(data,~lm(x_2 ~ ., data = .))) %>% 
    mutate(model_info = map(model, broom::augment))%>% 
    dplyr::select(-data, -model) %>%
    unnest(cols = c(model_info)) %>% 
    ungroup %>% 
    pull(.resid) %>% 
    sd
  
  initial_random_intercept_2 <- df_impute_2 %>%
    group_by(subject_id) %>% 
    nest() %>%
    mutate(model = map(data,~lm(x_2 ~ ., data = .))) %>% 
    mutate(model_info = map(model, broom::augment))%>% 
    dplyr::select(-data, -model) %>%
    unnest(cols = c(model_info)) %>% 
    ungroup %>% 
    pull(.resid) %>% 
    mean
  
  initial_tau_2 <- df_impute_2 %>%
    group_by(subject_id) %>%
    nest() %>%
    mutate(model = map(data,~lm(x_2 ~ ., data = .))) %>%
    mutate(model_info = map(model, broom::tidy)) %>%
    dplyr::select(-data, -model) %>%
    unnest(cols = c(model_info)) %>%
    ungroup %>%
    filter(term == "(Intercept)") %>%
    pull(estimate) %>%
    sd
  
  initial_delta_2 <- df_impute_2 %>%
    group_by(subject_id) %>%
    nest() %>%
    mutate(model = map(data,~lm(x_2 ~ ., data = .))) %>%
    mutate(model_info = map(model, broom::tidy)) %>%
    dplyr::select(-data, -model) %>%
    unnest(cols = c(model_info)) %>%
    ungroup %>%
    filter(term == "time") %>%
    pull(estimate) %>%
    sd
  # matrix for M-H to impute missingness
  x_1_h <- matrix(NA, nrow = number_of_MCMC, ncol = nrow_df)
  
  x_1_h[1,] <- df_reorder[,"x_1"]
  
  
  for (h in 1:number_of_MCMC){
    df_impute_1 <- cbind(df_reorder[1:2], subject_id = subject, time = time)
    
    rbart_vi_1 <-
      rbart_vi(
        x_1 ~ x_0,
        df_impute_1,
        group.by = subject_id,
        n.chains = 1L,
        n.threads = 1L,
        test = df_impute_1,
        group.by.test = subject_id,
        n.samples = 1,
        n.burn = 100L,
        keepTrees = T
      )
    
    # Update thetas
    rbart_vi_1_sigest <- rbart_vi_1$sigest
    rbart_vi_1_tau <- rbart_vi_1$tau
    rbart_vi_1_sigma <- rbart_vi_1$sigma
    rbart_vi_1_fit <- rbart_vi_1$fit
    rbart_vi_1_delta <- rbart_vi_1$delta
    
    df_impute_2 <- cbind(df_reorder[1:3], subject_id = subject, time = time)
    rbart_vi_2 <-
      rbart_vi(
        x_2 ~ x_0 + x_1,
        df_impute_2,
        group.by = subject_id,
        n.chains = 1L,
        n.threads = 1L,
        test = df_impute_2,
        group.by.test = subject_id,
        n.samples = 1,
        n.burn = 100L,
        keepTrees = T
      )
    # Update thetas
    df_impute_2_sigest <- rbart_vi_2$sigest
    df_impute_2_tau <- rbart_vi_2$tau
    df_impute_2_sigma <- rbart_vi_2$sigma
    df_impute_2_fit <- rbart_vi_2$fit
    df_impute_2_delta <- rbart_vi_2$delta
    
    # Impute missing
    df_reorder[,"x_1"] <- rbart_vi_1$yhat.test.mean
    df_reorder[,"x_2"] <- rbart_vi_2$yhat.test.mean
    result_x_2_imputation[,h] <- df_reorder[,"x_2"]
    # For x_1 use M-H sampling 
    df_impute_1_updated <- cbind(df_reorder[1:2], subject_id = subject, time = time)
    x_1_h <- as.numeric(predict(rbart_vi_1, newdata = df_impute_1_updated, group.by = df_impute_1_updated$subject_id))
    x_1_last_h <- rbart_vi_1$yhat.test.mean
    x_1_acccept <- pnorm(x_1_h, rbart_vi_1_sigma)/pnorm(x_1_last_h, rbart_vi_1_sigma)
    df_reorder[,"x_1"] <- ifelse(x_1_acccept > 1, x_1_h, ifelse(x_1_acccept > runif(1), x_1_h, x_1_last_h))
    result_x_1_imputation[,h] <- df_reorder[,"x_1"]
    print(h)
  }
  return(list(x_1 = result_x_1_imputation, x_2 = result_x_2_imputation))
}



