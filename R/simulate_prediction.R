library(invgamma)
library(mvtnorm)
library("stringr")
library(dplyr)
library(tidyr)
library(mnonr)
library(mice)
library(mixtools)
library(extraDistr)
library(sn)

sigmoid = function(x){
  1 / (1 + exp(-x))
}

simulation_prediction = function(n_subject = 800, n_obs_per_sub = 6, seed = 123, nonlinear = FALSE, nonrandeff = FALSE, nonresidual = FALSE){
  set.seed(123)
  n_obs_per_sub = sapply(1:n_subject, function(x) n_obs_per_sub)
  subject_id = c(unlist(sapply(1:n_subject, function(x) rep(x, n_obs_per_sub[x]))))
  n_obs = length(subject_id)
  Z = c(unlist(sapply(1:n_subject, function(x) 1:n_obs_per_sub[x])))
  trajectory = cbind(Z)
  Z = cbind(Z, Z^2)
  Z = cbind(rep(1, length(subject_id)), Z)
  Z = apply(Z[,-1], 2, scale)
  Z = cbind(1, Z)
  Z_O = cbind(Z)
  n = dim(Z)[1]

  
  X = rmvnorm(n_obs, c(0, 2, 5, -10, -3, -8, 9))
  
  X[,1] = apply(cbind(X[,1]), 2, function(X_v){
    Bi = rmvnorm(n_subject, rep(0, dim(Z)[2]))
    re = sapply(1:length(subject_id), function(x){
      Z[x,] %*% Bi[subject_id[x],]
    })
    X_v + re
  })
  
  
  X[,2] = apply(cbind(X[,2]), 2, function(X_v){
    Bi = rmvnorm(n_subject, rep(0, dim(Z)[2]))
    re = sapply(1:length(subject_id), function(x){
      Z[x,] %*% Bi[subject_id[x],]
    })
    X_v + re
  })
  
  X[,3] = apply(cbind(X[,3]), 2, function(X_v){
    Bi = rmvnorm(n_subject, rep(0, dim(Z)[2]))
    re = sapply(1:length(subject_id), function(x){
      Z[x,] %*% Bi[subject_id[x],]
    })
    X_v + re
  })
  
  X[,4] = apply(cbind(X[,4]), 2, function(X_v){
    Bi = rmvnorm(n_subject, rep(0, dim(Z)[2]))
    re = sapply(1:length(subject_id), function(x){
      Z[x,] %*% Bi[subject_id[x],]
    })
    X_v + re
  })
  
  X[,5] = apply(cbind(X[,5]), 2, function(X_v){
    Bi = rmvnorm(n_subject, rep(0, dim(Z)[2]))
    re = sapply(1:length(subject_id), function(x){
      Z[x,] %*% Bi[subject_id[x],]
    })
    X_v + re
  })
  
  X[,6] = apply(cbind(X[,6]), 2, function(X_v){
    Bi = rmvnorm(n_subject, rep(0, dim(Z)[2]))
    re = sapply(1:length(subject_id), function(x){
      Z[x,] %*% Bi[subject_id[x],]
    })
    X_v + re
  })
  
  X[,7] = apply(cbind(X[,7]), 2, function(X_v){
    Bi = rmvnorm(n_subject, rep(0, dim(Z)[2]))
    re = sapply(1:length(subject_id), function(x){
      Z[x,] %*% Bi[subject_id[x],]
    })
    X_v + re
  })
  
  Y = 5 + X[,1] - X[,2] + 3 * X[,3] - 2 * X[,4] + 0.5 * X[,5] - X[,6] + 0.5 * X[,7] 
  if(nonlinear)
    Y = 1 * (X[,3] * X[,4]) +        # Interaction1 (scaled down)
    1 * (X[,7] * X[,1]) +        # Interaction2 (scaled down)
    1 * log(3 * (X[,5] + X[,6])^2 + 0.5) +        # Simple addition interaction (scaled down)
    0.5 * X[,2]^2
  
  if(nonresidual){
    Y = Y + (-1)^rbinom(n_obs, 1, prob = 0.5) * rchisq(n_obs, 10)#rnormmix(n_obs, mu = c(-8, -3, 2, 4, 4.5), sigma = c(1, 1, 1, 1, 1), lambda = c(1/5, 1/5, 0.15, 0.25, 1/5))
  }else{
    Y = Y + rnorm(n_obs,sd = 5)
  }

  if (nonrandeff){
    
    # Define enhanced parameters
    xi <- c(0, 0, 0)  # Location vector
    
    
    Omega <- matrix(c(1.0, 0.5, 0.3,
                      0.5, 1.0, 0.4,
                      0.3, 0.4, 1.0),
                    nrow = 3, byrow = TRUE)
    alpha <- c(1, 1, 1) 
    Bi <- rmst(n = n_subject, xi = xi, Omega = Omega, alpha = alpha, nu = 5)
  }else{
    Bi = rmvnorm(n_subject, rep(0, dim(Z)[2]))
  }
  re = sapply(1:length(subject_id), function(x){
    Z[x,] %*% Bi[subject_id[x],]
  })
  Y = Y + re

  
  Y_copy = as.numeric(Y)
  
  
  data_O = tibble(
    subject_id = subject_id,
    time = c(unlist(sapply(1:n_subject, function(x) 1:n_obs_per_sub[x]))),
    X1 = X[,1],
    X2 = X[,2],
    X3 = X[,3],
    X4 = X[,4],
    X5 = X[,5],
    X6 = X[,6],
    X7 = X[,7],
    Y = Y_copy
  )
  
  data = data_O %>% 
    pivot_wider(id_cols = subject_id, names_from = time, values_from = c(X1, X2, X3, X4, X5, X6, X7, Y))
  
  
  
  
  patterns = matrix(1, nrow = n_obs_per_sub[1], ncol = n_obs_per_sub[1] * 8)
  patterns[1, 43] = 0
  patterns[2, 44] = 0
  patterns[3, 45] = 0
  patterns[4, 46] = 0
  patterns[5, 47] = 0
  patterns[6, 48] = 0
  
  for (i in 43:47) {
    for (j in (i + 1):48) {
      a = rep(1, n_obs_per_sub[1] * 8)
      a[i] = 0
      a[j] = 0
      #print(a)
      patterns = rbind(patterns, a)
    }
  }
  
  for (i in 43:46) {
    for (j in (i + 1):47) {
      for (k in (j + 1):48) {
        a = rep(1, n_obs_per_sub[1] * 8)
        a[i] = 0
        a[j] = 0
        a[k] = 0
        #print(a)
        patterns = rbind(patterns, a)
      }
    }
  }
  
  freq = rep(0.1, dim(patterns)[1])
  #freq = c(rep(0.3 / (n_obs_per_sub[1]), n_obs_per_sub[1]), rep(0.3 / 15, 15), rep(0.3 / 20, 20))
  set.seed(seed)
  a = mice::ampute(data[,-1], mech = "MCAR", patterns = patterns, prop = 0.9999, freq = freq)
 
  
  
  

  
  b = as_tibble(a$amp) %>%
    mutate(subject_id = data$subject_id) %>% 
    pivot_longer(cols = c(starts_with("X"), starts_with("Y")), names_prefix = "_", values_to = "value", names_to = "Variable") %>%
    separate(col = Variable, into = c("Variable", "Time"), sep = "_") %>% 
    pivot_wider(id_cols = c(subject_id, Time), names_from = "Variable", values_from = "value") %>% 
    mutate(followup = ifelse(
      is.na(X1) & is.na(X2) & is.na(X3) & is.na(X4) & is.na(X5) & is.na(X6) & is.na(X7) & is.na(Y),
      0,
      1
    ))
  
  
  subject_id = b$subject_id
  trajectory = b$Time
  Y = unlist(b[,10])
  #mean(is.na(Y))
  

  
  Y_O = Y_copy

  subject_id = data_O$subject_id
  
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
    #print(omit_sub)
    #print((unique(subject_id) %in% omit_sub))
    #print(Bi[! (unique(subject_id) %in% omit_sub)])
    warning(paste0("Some variables of ", length(omit_sub), " subjects are missing at all time points, we delete data from these ", length(omit_sub), " subjects:\n", paste(omit_sub, collapse = ", ")))
    
    X = X[!subject_id %in% omit_sub,]
    Y = Y[!subject_id %in% omit_sub]
    Y_copy = Y_copy[!subject_id %in% omit_sub]
    
    Z = Z[!subject_id %in% omit_sub,]
    Z_O = Z_O[!subject_id %in% omit_sub,]
    
    trajectory = trajectory[!subject_id %in% omit_sub]
    #print(!(unique(subject_id) %in% omit_sub))
    Bi = Bi[-omit_sub,]
    re = re[!subject_id %in% omit_sub]
    subject_id = subject_id[!subject_id %in% omit_sub]
    #print(Bi)
  }
  
  return(list(re = re, Y_O = Y_copy, Bi = Bi, Y = Y, X = X, Z = Z_O, subject_id = subject_id, trajectory = trajectory))
}



