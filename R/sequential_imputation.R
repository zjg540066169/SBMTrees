library(lme4)
library(optimx)
library(dplyr)
library(jomo)
sourceCpp("./src/sequential_imputation.cpp")



apply_locf_nocb <- function(X, subject_id) {
  locf <- function(x) {
    # Fill missing values by carrying forward the last non-missing value
    for (i in 2:length(x)) {
      if (is.na(x[i])) {
        x[i] <- x[i - 1]
      }
    }
    return(x)
  }
  
  nocb <- function(x) {
    # Fill missing values by carrying backward the next non-missing value
    for (i in (length(x) - 1):1) {
      if (is.na(x[i])) {
        x[i] <- x[i + 1]
      }
    }
    return(x)
  }
  
  apply_locf <- function(X, subject_id) {
    unique_ids <- unique(subject_id)
    X_imputed <- as.matrix(X)  # Copy of the matrix to store imputed values
    
    # Loop over each individual
    for (id in unique_ids) {
      # Get the rows corresponding to the current individual
      rows <- which(subject_id == id)
      
      # Apply LOCF for each covariate (column) for this individual
      for (col in 1:ncol(X)) {
        X_imputed[rows, col] <- locf(X[rows, col])
      }
    }
    
    return(X_imputed)
  }
  
  apply_nocb <- function(X, subject_id) {
    unique_ids <- unique(subject_id)
    X_imputed <- as.matrix(X)  # Copy of the matrix to store imputed values
    
    # Loop over each individual
    for (id in unique_ids) {
      # Get the rows corresponding to the current individual
      rows <- which(subject_id == id)
      
      # Apply NOCB for each covariate (column) for this individual
      for (col in 1:ncol(X)) {
        X_imputed[rows, col] <- nocb(X[rows, col])
      }
    }
    
    return(X_imputed)
  }
  
  X_locf <- apply_locf(X, subject_id)  # Apply LOCF first
  X_locf_nocb <- apply_nocb(X_locf, subject_id)  # Then apply NOCB
  return(X_locf_nocb)
}




sequential_imputation <- function(X, Y,  type, Z = NULL, subject_id, model = c("BMTrees", "BMTrees_R", "BMTrees_RE", "mixedBART"), nburn = 0L, npost = 3L, skip = 1L, verbose = TRUE, seed = NULL, tol = 1e-20, resample = 0, ntrees = 200, reordering = T) {
  model = match.arg(model)
  if(is.null(dim(X))){
    stop("More than one covariate is needed!")
  }
  
  
  omit_sub = sapply(unique(subject_id), function(sub){
    t = sum(subject_id == sub)
    X_sub = cbind(X[subject_id == sub,], Y[subject_id == sub])
    if (sum(colSums(is.na(X_sub)) == t) > 0){
      return(sub)
    }
    return(NA)
  })
  omit_sub = omit_sub[!is.na(omit_sub)]
  if(length(omit_sub) > 0){
    warning(paste0("Some variables of ", length(omit_sub), " subjects are missing at all time points, we delete data from these ", length(omit_sub), " subjects:\n", paste(omit_sub, collapse = ", ")))
    X = X[!subject_id %in% omit_sub,]
    Y = Y[!subject_id %in% omit_sub]
    if (!is.null(Z)){
      Z = Z[!subject_id %in% omit_sub,]
    }
    subject_id = subject_id[!subject_id %in% omit_sub]
  }
  trajectory = c(sapply(table(subject_id), function(i) 1:i))
  subject_id = as.factor(subject_id)
  p = dim(X)[2]
  
  mis_num = colSums(is.na(X))
  if(reordering == T){
    mis_num[mis_num ==0] = -Inf
    mis_num[mis_num > 0 & type == 1] = mis_num[mis_num > 0 & type == 1] - base::max(mis_num)
    mis_order = order(mis_num)
    re_order = seq(ncol(is.na(X)))
    type = type[mis_order]
    X[,re_order] = X[,mis_order]
    colnames(X) = colnames(X)[mis_order]
    cat("reordering: new covariates order is ", colnames(X), "\n")
  }
  

  R = is.na(X)
  R = cbind(R, "Y" = is.na(Y))
  X_locf_nocb <- apply_locf_nocb(cbind(X, Y), subject_id)
  Y = X_locf_nocb[,dim(X_locf_nocb)[2]]
  X = X_locf_nocb[,1:(dim(X_locf_nocb)[2] - 1)]

  if(sum(mis_num == -Inf) == 0){
    X = cbind("intercept" = 1, X)
    type = c(0, type)
    R = cbind(0, R)
  }
  
  
 
  cat("Complete initializing imputing missing data\n")
  cat("Start to impute by Sequential Imputation\n")
 
  if(model == "BMTrees_R"){
    imputation_X_DP = sequential_imputation_cpp(as.matrix(X), as.numeric(Y), as.logical(type), as.matrix(Z), as.character(subject_id), as.matrix(R), binary_outcome = FALSE, nburn = nburn, npost = npost, skip = skip, verbose = verbose, CDP_residual = TRUE, CDP_re = FALSE, seed = seed, ncores = 0, ntrees = ntrees, fit_loss = F, resample = resample)
  }
  else if(model == "BMTrees_RE"){
    imputation_X_DP = sequential_imputation_cpp(as.matrix(X), as.numeric(Y), as.logical(type), as.matrix(Z), as.character(subject_id), as.matrix(R), binary_outcome = FALSE, nburn = nburn, npost = npost, skip = skip, verbose = verbose, CDP_residual = FALSE, CDP_re = TRUE, seed = seed, ncores = 0, ntrees = ntrees, fit_loss = F, resample = resample)
  }
  else if(model == "BMTrees"){
    imputation_X_DP = sequential_imputation_cpp(as.matrix(X), as.numeric(Y), as.logical(type), as.matrix(Z), as.character(subject_id), as.matrix(R), binary_outcome = FALSE, nburn = nburn, npost = npost, skip = skip, verbose = verbose, CDP_residual = TRUE, CDP_re = TRUE, seed = seed, ncores = 0, ntrees = ntrees, fit_loss = F, resample = resample)
  }
  else if(model == "mixedBART"){
    imputation_X_DP = sequential_imputation_cpp(as.matrix(X), as.numeric(Y), as.logical(type), as.matrix(Z), as.character(subject_id), as.matrix(R), binary_outcome = FALSE, nburn = nburn, npost = npost, skip = skip, verbose = verbose, CDP_residual = FALSE, CDP_re = FALSE, seed = seed, ncores = 0,  ntrees = ntrees, fit_loss = F, resample = resample)
  }
  else{
    imputation_X_DP = sequential_imputation_cpp(as.matrix(X), as.numeric(Y), as.logical(type), as.matrix(Z), as.character(subject_id), as.matrix(R), binary_outcome = FALSE, nburn = nburn, npost = npost, skip = skip, verbose = verbose, CDP_residual = TRUE, CDP_re = TRUE, seed = seed, ncores = 0, ntrees = ntrees, fit_loss = F, resample = resample)
  }
  
  imputation_Y = t(do.call(cbind, imputation_X_DP$imputation_Y_DP))
  if(reordering == T){
    imputation_X_DP = lapply(imputation_X_DP$imputation_X_DP, function(x){
      if(sum(mis_num == -Inf) == 0){
        x = x[,-1]
      }
      x[,mis_order] = x[,re_order]
      colnames(x)[mis_order] = colnames(x)[re_order]
      x
    })
  }else{
    imputation_X_DP = lapply(imputation_X_DP$imputation_X_DP, function(x){
      if(sum(mis_num == -Inf) == 0){
        x = x[,-1]
      }
      x
    })
  }
  
  imputation_X = array(NA, dim = c(length(imputation_X_DP), dim(imputation_X_DP[[1]])[1], dim(imputation_X_DP[[1]])[2] + 1))
  for (i in 1:length(imputation_X_DP)) {
    imputation_X[i,,] = cbind(imputation_X_DP[[i]], imputation_Y[i,])
  }
  return(list(imputed_data = imputation_X))
}
