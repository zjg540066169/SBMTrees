library(invgamma)
library(mvtnorm)
library("stringr")
library(dplyr)
library(tidyr)
library(mnonr)
library(mice)
library(mixtools)
library(extraDistr)

sigmoid = function(x){
  1 / (1 + exp(-x))
}



# library(parallel)
# library(doParallel) a
# registerDoParallel(4) 
# 
# sim = foreach(i = 1:1000, .combine = cbind) %dopar% {
#   d = generate_longitudinal_MAR(seed = i)
#   colMeans(is.na(d$X_mis))[7:10]
# }

# simulation scenario
# 1: nonlinear normal-random-effect nonnormal-residual
# 2: linear normal-random-effect nonnormal-residual

simulation = function(n_subject = 1000, n_obs_per_sub = 5, seed = 123, nonlinear = FALSE, nonrandeff = FALSE, nonresidual = FALSE, alligned = F, survival = FALSE ){
  set.seed(seed)
  #n_obs_per_sub = sapply(1:n_subject, function(x) sample(n_obs_per_sub, 1))
  n_obs_per_sub = sapply(1:n_subject, function(x) n_obs_per_sub)
  subject_id = c(unlist(sapply(1:n_subject, function(x) rep(x, n_obs_per_sub[x]))))
  n_obs = length(subject_id)
  Z = c(unlist(sapply(1:n_subject, function(x) 1:n_obs_per_sub[x])))
  trajectory = cbind(Z)
  Z = cbind(Z, Z^2)
  Z = cbind(rep(1, length(subject_id)), Z)
  Z_O = cbind(Z)
  #Z_O = cbind(Z)
  #Z = apply(Z[,-1], 2, function(z){(z - min(z))/(max(z) - min(z))})
  Z = apply(Z[,-1], 2, scale)
  #Z = cbind(1, Z)
  #print(Z)
  #Z = cbind(rep(1, length(subject_id)))
  n = dim(Z)[1]

  
  X_1 = cbind(rnorm(n_obs, rep(0, 2)))
  X_1 = apply(X_1, 2, function(X_v){
    Bi = rmvnorm(n_subject, rep(0, dim(Z)[2]))
    re = sapply(1:length(subject_id), function(x){
      Z[x,] %*% Bi[subject_id[x],]
    })
    X_v + re
  })
  X_1 = cbind(rbinom(n, 1, sigmoid(X_1)))
  
  X_2 = X_1 * 2 
  if(nonlinear)
    X_2 = exp(X_1 * -2)
  
  X_2 = apply(X_2, 2, function(X_v){
    if (nonrandeff){
      Bi = rmvt(n_subject, sigma = diag(dim(Z)[2]), df = 3)
    }else{
      Bi = rmvnorm(n_subject, rep(0, dim(Z)[2]))
    }
    re = sapply(1:length(subject_id), function(x){
      Z[x,] %*% Bi[subject_id[x],]
    })
    X_v + re
  })
  X_2 = cbind(rbinom(n, 1, sigmoid(X_2)))
  
  X_3 = X_1 * -1 + 0.5 * X_2
  if(nonlinear)
    X_3 = (X_1 - 2 * X_2)^2
    #X_3 = X_3 + (X_1 - 0.2 * X_2)^3 * 0.1
  
  X_3 = apply(X_3, 2, function(X_v){
    if (nonrandeff){
      Bi = rmvt(n_subject, sigma = diag(dim(Z)[2]), df = 10)
    }else{
      Bi = rmvnorm(n_subject, rep(0, dim(Z)[2]))
    }
    re = sapply(1:length(subject_id), function(x){
      Z[x,] %*% Bi[subject_id[x],]
    })
    X_v + re
  })
  
  
  X_4 = X_1 * 0.2 + 2 * X_2 + 3 * X_3 
  if(nonlinear)
    X_4 = 2 * (X_1 + X_2 - 2 * X_1 * X_3 + X_3 - 1 * log(X_3^2 + 0.5))
    #X_4 = 2 * (X_1 + X_2 - 0.2 * X_3 )^2
    #X_4 = 5 * sigmoid((X_1 + X_2 + 4 * X_3 + 20) / 8)
  
  
  if(nonresidual){
    X_4 = X_4 + runif(n_obs, 0, 2)
    #X_4 = X_4 + rnormmix(n_obs, mu = c(-10, -5, 0, 5, 10), sigma = c(1, 1, 1, 1, 1), lambda = c(1/5, 1/5, 1/5, 1/5, 1/5))#(-1)^rbinom(n_obs, 1, prob = 0.5) * rgamma(n_obs, 2, 2) * X_4/1
  }else{
    X_4 = X_4 + rnorm(n_obs)
    #X_4 = X_4 + runif(n_obs, -1, 1)
    #X_4 = X_4 + rnormmix(n_obs, mu = c(-10, 0, 10), sigma = c(0.5, 0.5, 0.5), lambda = c(1/3, 1/3, 1/3))
  }
  X_4 = apply(X_4, 2, function(X_v){
    if (nonrandeff){
      Bi = sapply(1:(dim(Z)[2]), function(i){
        (-1)^rbinom(n_subject, 1, prob = 0.5) * rchisq(n_subject, 1)
      })
    }else{
      Bi = rmvnorm(n_subject, rep(0, dim(Z)[2]))
    }
    re = sapply(1:length(subject_id), function(x){
      Z[x,] %*% Bi[subject_id[x],]
    })
    X_v + re
  })
  
  X_5 = X_1 - 5 * X_2 - X_3 + X_4
  if(nonlinear)
    X_5 =  -(X_1 * X_2 - 0.3 * X_3^2 - X_2 * X_4)#-(X_1 + X_2 - 2 * X_3 - X_2 * X_4)
  if(nonresidual){
    X_5 = X_5 + rgamma(n_obs, 2, 2)
  }else{
    X_5 = X_5 + rnorm(n_obs)
  }
  X_5 = apply(X_5, 2, function(X_v){
    if (nonrandeff){
      Bi = sapply(1:(dim(Z)[2]), function(i){
        rgamma(n_subject, 2, 1)
      })
    }else{
      Bi = rmvnorm(n_subject, rep(0, dim(Z)[2]))
    }
    re = sapply(1:length(subject_id), function(x){
      Z[x,] %*% Bi[subject_id[x],]
    })
    X_v + re
  })
  
  
  
  
  
  X_6 = X_1 * -4 + 5 * X_2 - X_3 + X_4  + X_5
  if(nonlinear)
    #X_6 = X_6 - log(abs(X_1 - 3 * X_2 + 2 * X_3 - 0.3 * X_4 - X_5) + 1)
    #X_6 = - 5 * log(abs(X_1 - 3 * X_2 + 3 * X_3 + X_4 * X_2 - X_5) + 1) + 10
    X_6 = - 10 * log((X_1 - 3 * X_2 + 3 * X_3 + X_4 * X_2 - X_5)^2 + 10) + 20
  if(nonresidual){
    X_6 = X_6 + rexp(n_obs, 1)
    #X_6 = X_6 + (-1)^rbinom(n_obs, 1, prob = 0.5) * rchisq(n_obs, 8)
  }else{
    X_6 = X_6 + rnorm(n_obs)
  }
  X_6 = apply(X_6, 2, function(X_v){
    if (nonrandeff){
      Bi = sapply(1:(dim(Z)[2]), function(i){
        (-1)^rbinom(n_subject, 1, prob = 0.5) * rgamma(n_subject, 10, 10)
      })
    }else{
      Bi = rmvnorm(n_subject, rep(0, dim(Z)[2]))
    }
    re = sapply(1:length(subject_id), function(x){
      Z[x,] %*% Bi[subject_id[x],]
    })
    X_v + re
  })

  X_7 = 0.5 * X_1 - 3 * X_2 + X_3 - 1.5 * X_4 - X_5 + 0.5 * X_6 
  if(nonlinear)
    X_7 = (X_1 + X_2) * X_6 + 5 * log(X_5^2 + 1) + (X_4 - X_3) * X_1 - 20
    #X_7 = (X_1 + X_2) * X_4 + 5 * log(X_5^2 + 1) - sigmoid(X_6) +  abs(X_4 - X_3)
  if(nonresidual){
    #X_7 = X_7 + (-1)^rbinom(n_obs, 1, prob = 0.5) * rgamma(n_obs, 2, 2)
    X_7 = X_7 + (-1)^rbinom(n_obs, 1, prob = 0.2) * rgamma(n_obs, 8, 3)
  }else{
    X_7 = X_7 + rnorm(n_obs)
  }
  X_7 = apply(X_7, 2, function(X_v){
    if (nonrandeff){
      Bi = sapply(1:(dim(Z)[2]), function(i){
        rexp(n_subject, 1)
        #rnorm(n_subject) * rgamma(n_subject, 2, 2) + runif(n_subject, -2, 2)
      })
    }else{
      Bi = rmvnorm(n_subject, rep(0, dim(Z)[2]))
    }
    re = sapply(1:length(subject_id), function(x){
      Z[x,] %*% Bi[subject_id[x],]
    })
    X_v + re
  })
  
  Y = X_1 - X_2 + 3 * X_3 - 2 * X_4 + 0.5 * X_5 - X_6 + 0.5 * X_7 
  if(nonlinear)
    Y = X_1 - X_2 + 3 * X_3 - 2 * X_4 + 0.5 * X_5 - X_6 + 0.5 * X_7 
    #Y = Y + sqrt(abs(X_1 - X_2 + X_3 - X_4 - X_5 + X_6 + X_7) + 1)
    # Y = 0.5 * (X_2 * X_4) +        # Interaction1 (scaled down)
    # 0.5 * (X_7 * X_1) +        # Interaction2 (scaled down)
    # 0.5 * (X_5 + X_6) +        # Simple addition interaction (scaled down)
    # 0.5 * X_3^2  
   # Y = 3 * sin(X_1 + X_3*pi + X_2) +  tanh(X_4) * sqrt(abs(X_6) + 10) +  log(X_5^2) - log(abs(X_7))#sqrt(abs(4 * X_1 - 3 * X_2 + X_3 * X_4 + X_6 + X_7) + 1)# - atan(X_5 * X_7) * sigmoid(X_4) * X_3 + X_4  * cos(X_4 / 10) + abs(X_6)
  Y_copy = as.numeric(Y)
  if(nonresidual){
    #noise = rexp(n_obs)#//(-1)^rbinom(n_obs, 1, prob = 0.5) * rgamma(n_obs, 7, 2) * exp(-abs(Y)/50)
    #hist(noise, breaks = n_obs / 3)
    #Y = Y + noise
    #Y = Y + rnormmix(n_obs, mu = c(-8, -4, 0, 4, 8), sigma = c(1, 1, 1, 1, 1), lambda = c(1/5, 1/5, 1/5, 1/5, 1/5))
    Y = Y + rnorm(n_obs, -3, 1) * rgamma(n_obs, 5, 3)
  }else{
    Y = Y + rnorm(n_obs,sd = 1)
  }
  
  
  
  if (nonrandeff){
    
    Bi = sapply(1:(dim(Z)[2]), function(i){
      #rnorm(n_subject) * rgamma(n_subject, 2, 2) + runif(n_subject, -2, 2)
      rnorm(n_subject, -1, 1) * runif(n_subject, 0, 2)
    })
    
    # Bi = sapply(1:(dim(Z)[2]), function(i){
    #   #runif(n_subject, -5, 5)
    #   #(-1)^rbinom(n_subject, 1, prob = 0.5) * rgamma(n_subject, 7, 2)
    #   rexp(n_subject, 1)
    #   #rnormmix(n_subject, mu = c(-8, -4, 0, 4, 8), sigma = c(1, 1, 1, 1, 1), lambda = c(1/5, 1/5, 1/5, 1/5, 1/5))
    # })
    # lambda <- rep(1, 4)/4
    # mu <- matrix(2*(1:12)* rbinom(12, 1, 0.5), 4, 3) 
    # sigma <- matrix(1,4,3)
    # Bi = rmvnormmix(n_subject, lambda, mu, sigma)
   
    # Bi = sapply(1:(dim(Z)[2]), function(i){
    #   #runif(n_subject, -5, 5)
    #   (-1)^rbinom(n_subject, 1, prob = 0.5) * rgamma(n_subject, 7, 2)
    #   #rnormmix(n_subject, mu = c(-8, -4, 0, 4, 8)+5, sigma = c(1, 1, 1, 1, 1), lambda = c(1/5, 1/5, 1/5, 1/5, 1/5))
    # })
  }else{
    Bi = rmvnorm(n_subject, rep(0, dim(Z)[2]))
  }
  re = sapply(1:length(subject_id), function(x){
    Z[x,] %*% Bi[subject_id[x],]
  })
  #print(re)
  Y = Y + re
 
  
  data_O = tibble(
    subject_id = subject_id,
    time = c(unlist(sapply(1:n_subject, function(x) 1:n_obs_per_sub[x]))),
    X1 = X_1[,1],
    X2 = X_2[,1],
    X3 = X_3[,1],
    X4 = X_4[,1],
    X5 = X_5[,1],
    X6 = X_6[,1],
    X7 = X_7[,1],
    Y = Y[,1]
  )
  
  if(!alligned){
    data_O = data_O %>% 
      dplyr::select(subject_id,  time, X2,   X1, X6, X5, X4, X7, X3,Y)
  }
  
  data = data_O %>% 
    pivot_wider(id_cols = subject_id, names_from = time, values_from = c(X1, X2, X3, X4, X5, X6, X7, Y))
  
  
  
  patterns = matrix(1, nrow = n_obs_per_sub[1] * 8, ncol = n_obs_per_sub[1] * 8)
  diag(patterns) = 0
  patterns = rbind(patterns, t(sapply(1:n_obs_per_sub[1], function(x){
    d = rep(rep(1,8), n_obs_per_sub[1])
    d[5 * (1:8 - 1) + x] = 0
    d
  })))
  
  patterns = rbind(patterns,t(sapply(4:(n_obs_per_sub[1] - 1), function(x){
    d = rep(rep(1,8), n_obs_per_sub[1])
    for (i in x:n_obs_per_sub[1]){
      d[5 * (1:8 - 1) + i] = 0
    }
    
    d
  })))
  
  
  patterns = rbind(patterns, t(sapply(1:100, function(i){
    rbinom(dim(patterns)[2], 1, 0.7)
  })))         
  
  #freq = c(rep(0.2 / 40, 40), rep(0.1 / 4, 4), rep(0.5 / 4, 4), rep(0.5 / 40, 40))
  
  #freq = c(rep(0.0 / 40, 40), rep(0.0 / 4, 4), rep(0.0 / 4, 4), rep(1 / 100, 100))
  freq = c(rep(0.2 / 40, 40), rep(0.1 / 4, 4), rep(0.2 / 2, 2), rep(0.2 / 100, 100))
  if(nonrandeff){
    weights = matrix(0, nrow = dim(patterns)[1], ncol = dim(patterns)[2])
    for (i in 1:5) {
      for (p in 1:8) {
        weights[, (p - 1) * 5 + i] = abs(rnorm(dim(patterns)[1], (10 + i + p) * 1000, 10000))
      }
    }
    a = mice::ampute(data[,-1], mech = "MAR", patterns = patterns, prop = 0.8, freq = freq, type = "TAIL", weights = weights)
  }else{
    weights = matrix(0, nrow = dim(patterns)[1], ncol = dim(patterns)[2])
    for (i in 1:5) {
      for (p in 1:8) {
        weights[, (p - 1) * 5 + i] = abs(rnorm(dim(patterns)[1], (10 + i + p) * 100, 10000))
      }
    }
    a = mice::ampute(data[,-1], mech = "MAR", patterns = patterns, prop = 0.8, freq = freq, type = "LEFT", weights = weights)
  }
  
  
  if(nonrandeff == F & F == nonlinear & F == nonresidual){
    weights = matrix(0, nrow = dim(patterns)[1], ncol = dim(patterns)[2])
    for (i in 1:5) {
      for (p in 1:8) {
        weights[, (p - 1) * 5 + i] = abs(rnorm(dim(patterns)[1], (10 + i + p) * 10000, 10000))
      }
    }
    a = mice::ampute(data[,-1], mech = "MAR", patterns = patterns, prop = 0.8, freq = freq, type = "LEFT", weights = weights, std = F)
  }
  
  
  
  
  
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
  
  
  
  if(nonresidual){
    X_3 = X_3 + rnormmix(n_obs, mu = c(-10, -5, 0, 5, 10), sigma = c(1, 1, 1, 1, 1), lambda = c(1/5, 1/5, 1/5, 1/5, 1/5))#(-1)^rbinom(n_obs, 1, prob = 0.5) * rgamma(n_obs, 2, 2) * X_3/1
  }else{
    X_3 = X_3 + rnorm(n_obs)
  }
  
  
  X_O = data_O[,3:9]
  X_mis = b[,3:9]
  if(!alligned){
    X_mis = X_mis %>% 
      dplyr::select( X2,   X1, X6, X5, X4, X7, X3)
  }
  
  
  Y_O = data_O$Y
  Y = b$Y
  subject_id = data_O$subject_id
  
  omit_sub = sapply(unique(subject_id), function(sub){
    t = sum(subject_id == sub)
    X_sub = cbind(X_mis[subject_id == sub,], Y[subject_id == sub])
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
    X_mis = X_mis[!subject_id %in% omit_sub,]
    X_O = X_O[!subject_id %in% omit_sub,]
    Y = Y[!subject_id %in% omit_sub]
    Y_O = Y_O[!subject_id %in% omit_sub]
    Y_copy = Y_copy[!subject_id %in% omit_sub]
    
    Z_O = Z_O[!subject_id %in% omit_sub,]
    
    
    trajectory = trajectory[!subject_id %in% omit_sub]
    #print(!(unique(subject_id) %in% omit_sub))
    Bi = Bi[-omit_sub,]
    re = re[!subject_id %in% omit_sub]
    subject_id = subject_id[!subject_id %in% omit_sub]
    #print(Bi)
  }
  
  
  if(survival)
    return(list(Y = Y, X_O = X_O, X_mis = X_mis, Z = Z_O, subject_id = subject_id, trajectory = trajectory, censor = censor, censor_right_proportion = sum(censor == 1) / length(unique(subject_id)), censor_end_proportion = sum(censor == 2) / length(unique(subject_id)),  keep_proportion = mean(a$keep)))
  return(list(re = re, Y_copy = Y_copy, Bi = Bi, Y = Y, Y_O = Y_O, X_O = X_O, X_mis = X_mis, Z = Z_O, subject_id = subject_id, trajectory = trajectory))
}

