library(MASS)
library(BART, quietly = TRUE)
library(invgamma)
library(MCMCpack)
library(mvtnorm)
library(svMisc)
library(parallel) # one of the core R packages
library(doParallel)
library(DescTools)
suppressMessages(library(BART))




DP = function(G0, M, n = 1000, dp = 1) {
  b <- rbeta(n, 1, M)
  p <- numeric(n)
  p[1] <- b[1]
  p[2:n] <- sapply(2:n, function(i) b[i] * prod(1 - b[1:(i-1)]))
  #print(sum(p))
  y <- G0(n, dp)
  if (dp == 1){
    y <- y - sum(y * p)
    names(y) = seq(length(p))
    theta <- sample(y, prob = p, replace = TRUE)
  }
  else{
    y <- y - colSums(y * p)
    rownames(y) = seq(length(p))
    theta <- sample(rownames(y), prob = p, replace = TRUE)
    theta = y[theta,]
  }
  return(list(theta = theta, pi = p, h = y))
}

CDP = function(G0, M, n = 1000, dp = 1) {
  b <- rbeta(n, 1, M)
  p <- numeric(n)
  p[1] <- b[1]
  p[2:n] <- sapply(2:n, function(i) b[i] * prod(1 - b[1:(i-1)]))
  y <- G0(n, dp)
  if (dp == 1){
    y <- y - sum(y * p)
    names(y) = seq(length(p))
    theta <- sample(y, prob = p, replace = TRUE)
  }
  else{
    y <- y - colSums(y * p)
    rownames(y) = seq(length(p))
    theta <- sample(rownames(y), prob = p, replace = TRUE)
    theta = y[theta,]
  }
  return(list(theta = theta, pi = p, h = y))
}

G_0 <- function(n, p) {
  if (p == 1)
    rnorm(n, 0, 1)
  else{
    rmvnorm(n, rep(0, p))
  }
}
hist(CDP(G_0, 1000, dp = 1)$theta)



source("data_generate.R")

set.seed(111)
data <- data_generation(sigma_1  = 1, tau_1 = 2, delta_1 = 1, sigma_2  = 1, tau_2 = 2, delta_2 = 1, n_total = 5000, number_of_time_points = 10)



v = 3
k = 1

phi_1 = 2
phi_2 = 0.1


y = data$x_2
X = data[,c(1,2)]
z = data[, "time"]
p = dim(X)[2]
n = dim(X)[1]
subject_id = data[,"subject_id"]
subject_id = subject_id + 1
df = NULL
inverse_wishart_matrix = NULL
Covariance_tau = NULL
Mu_tau = NULL
tol = 1e-20



DP_LMM_BART = function(nburn, npost, y, X, z, subject_id, v = 3, k, phi_1 = 2, phi_2 = 0.1, sigma_tau = 1, df = NULL, inverse_wishart_matrix = NULL, Covariance_tau = NULL, Mu_tau = NULL, tol = 1e-20){
  nCores <- detectCores()  # to set manually 
  registerDoParallel(nCores) 
  if(is.null(dim(z))){d = 1} else{d = dim(z)[2]}
  n = length(y)
  p = dim(X)[2]
  n_subject = length(unique(subject_id))
  Z = matrix(1, nrow = n, ncol = d + 1)
  Z[, 2:(d+1)] = as.matrix(z)
  subject_to_B = seq(n_subject)
  names(subject_to_B) = unique(subject_id)
  
  
  # initial parameters
  if (is.null(inverse_wishart_matrix))
    inverse_wishart_matrix = diag(rep(1, d + 1))
  initial_M = rgamma(1, phi_1, phi_2)
  if (is.null(Covariance_tau))
    Covariance_tau = diag(d + 1)
  if (is.null(Mu_tau))
    Mu_tau = rep(0, d + 1)
  if (is.null(df))
    df = d + 1
  dp = CDP(function(n, p) {
    if (p == 1)
      rnorm(n, Mu_tau, Covariance_tau)
    else{
      rmvnorm(n, Mu_tau, Covariance_tau)
    }
  }, initial_M, n = n_subject, dp = d + 1)
  initial_location = dp$theta
  initial_pi = dp$pi
  initial_cluster = as.numeric(rownames(initial_location))
  initial_unique_location = initial_location[unique(rownames(initial_location)),]


  # posterior sampling
  post_M = array(NA, dim = npost)
  post_sigma = array(NA, dim = npost)
  post_Sigma = array(NA, dim = c(npost, (d + 1) * (d + 1)))
  post_B = array(NA, dim = c(npost, n_subject * (d + 1)))
  post_location = array(NA, dim = c(npost, n_subject * (d + 1)))
  post_unique_location = matrix(0, nrow = (npost), ncol = length(initial_unique_location))
  post_tree_hat = array(NA, dim = c(npost, length(y)))
  post_y_predict = array(NA, dim = c(npost, length(y)))
  post_tree = list()

  
  # initial values
  M = initial_M
  Covariance = riwish(df, inverse_wishart_matrix)
  B = t(apply(initial_location, 1, function (x) rmvnorm(1, x, Covariance)))
  location = initial_location
  unique_location = initial_unique_location
  cluster = initial_cluster
  pi = initial_pi
  names(pi) = cluster
  pi = t(rowsum(pi, group = names(pi)))[1,]
  
  for(kl in 1:(nburn + npost)){
    # burn in
    
    
    mixed_effect <- foreach(j = seq(length(y)), .combine = rbind) %dopar% {
      sub = subject_id[j]
      b_sub = subject_to_B[as.character(sub)]
      Z[j,] %*% B[b_sub,]
    }
    
    
    tree = wbart(X, y - mixed_effect, ndpost = 1, nskip = 10) # get tree
    
    y_hat = tree$yhat.train.mean
    sigma = tree$sigma[11]
    
    
    
    # prob = t(sapply(unique(subject_id), function(x) {
    #   b_sub = subject_to_B[as.character(x)]
    #   b_i = B[b_sub,]
    #   Z_i = Z[subject_id == x,]
    #   R_i = y[subject_id == x] - y_hat[subject_id == x]
    #   mean_i = t(Z_i %*% t(unique_location))
    #   var_i = Z_i %*% Covariance %*% t(Z_i)
    #   apply(mean_i, 1, function(y) pmvnorm(upper=R_i, mean = as.numeric(y), sigma=var_i)[1])
    #   }))
    
    
    
    prob = foreach(j = unique(subject_id), .combine = rbind) %dopar% {
      b_sub = subject_to_B[as.character(j)]
      b_i = B[b_sub,]
      proba = foreach(j2 = 1:(dim(unique_location)[1]), .combine = c) %dopar% {
        pmvnorm(upper=b_i, mean = as.numeric(unique_location[j2,]), sigma=Covariance)[1]
      }
      names(proba) = rownames(unique_location)
      proba
    }
    rownames(prob) = NULL
    prob[prob == 0] = 1e-20
    
    
    
    sample_frame = prob[,names(pi)] * pi
    cluster = apply(sample_frame, 1, function(x) sample(names(x), size = 1, replace = TRUE, prob = x)) 
    location = unique_location[cluster, ]
    
    
    
    # update B
    B = foreach(j = unique(subject_id), .combine = rbind) %dopar% {
      b_sub = subject_to_B[as.character(j)]
      Zi = Z[subject_id == as.numeric(j), ]
      R_i = y[subject_id == j] - y_hat[subject_id == j]
      var = solve(solve(Covariance) + t(Zi) %*% (Zi) / (sigma^2))
      mu = (t(t(Zi) %*% R_i) / (sigma^2) + location[b_sub,] %*% solve(Covariance))  %*% var 
      mvrnorm(1, mu, var)
    }
    rownames(B) = NULL
    
    nk = table(cluster)
    nk[setdiff(rownames(unique_location), names(nk))] = 0
    
   
    
    
    
    
    
    beta = matrix(0, length(nk) - 1, 2)
    beta[,1] = nk[1:(length(nk) - 1)] + 1
    beta[,2] = M + (rev(cumsum(rev(nk))) - nk)[1:(length(nk) - 1)]
    Vh = c(rbeta(length(nk) - 1, beta[,1], beta[,2]), 1)
    pi = numeric(length(nk))
    pi[1] = Vh[1]
    if(length(nk) >= 2)
      pi[2:length(nk)] = sapply(2:length(nk), function(i) Vh[i] * prod(1 - Vh[1:(i-1)]))
    names(pi) = names(nk)
    
    # update M
    phi_1_new = phi_1 + length(nk) - 1
    if(length(nk) >= 2)
      phi_2_new = phi_2 - sum(log(1 - Vh[1:(length(nk) - 1)]))
    M = rgamma(1, phi_1_new, phi_2_new)
    
    
    unique_location[,] = NA
    unique_location = foreach(j = names(nk), .combine = rbind) %dopar% {
      if (nk[j] == 0){
        new = rmvnorm(1, Mu_tau, Covariance_tau)#, Covariance)
        rownames(new) = j
        new
      }
        
      else{
        var_k = solve(solve(Covariance_tau) + nk[j] * solve(Covariance))
        if (nk[j] == 1)
          mean_k = (B[rownames(location) == j,] %*% solve(Covariance) + Mu_tau %*% solve(Covariance_tau)) %*% var_k
        
        
        #(B[unname(subject_to_B[as.character(unique(subject_id[cluster == as.numeric(j)]))]),] %*% solve(Covariance)) %*% var_k
        if (nk[j] > 1)
          mean_k = (colSums(B[rownames(location) == j,] %*% solve(Covariance)) + Mu_tau %*% solve(Covariance_tau)) %*% var_k
        
        
        #(colSums(B[unname(subject_to_B[as.character(unique(subject_id[cluster == as.numeric(j)]))]),]) %*% solve(Covariance)) %*% var_k
        new = rmvnorm(1, mean_k, var_k)#, Covariance)
        rownames(new) = j
        new
      }
    }
    rownames(unique_location) = names(nk)
    unique_location = unique_location - colSums(unique_location[rownames(unique_location), ] * pi[rownames(unique_location)])
    location_ = t(sapply(rownames(location), function(x) unique_location[x,]))
    rownames(location_) = rownames(location)
    location = location_
    
    
    # update Covariance matrix
    Covariance = riwish(df + n_subject, inverse_wishart_matrix + t(B - location) %*% (B - location))
    

    
    

    
    # bi =  B[i,]
    # bi ~ iid MVN(mu_i, Covariance)
    # mu_i ~ DP(M, G0)
    # G0 ~ MVN(mu_tau, Covariance_tau)
    
    
    
    
    if(kl > nburn){
      post_Sigma[kl - nburn,] = Covariance
      post_sigma[kl - nburn] = sigma
      post_B[kl - nburn,] = B
      post_unique_location[kl - nburn,] = unique_location
      post_location[kl - nburn,] = location
      post_tree_hat[kl - nburn, ] = y_hat
      post_y_predict[kl - nburn, ] = pwbart(X, tree$treedraws)[1,] + mixed_effect
      post_tree[[kl - nburn]] = tree$treedraws
    }
    progress(kl, npost + nburn)
  }
  return(list(tree_hat = post_tree_hat, y_predict = post_y_predict, tree = post_tree, y_hat = post_y_hat, Sigma = post_Sigma, sigma = post_sigma, B = post_B, location = post_location, unique_location = post_unique_location))
}









chain1 = DP_LMM_BART(5, 2, y, X, z, subject_id, v = 3, k, phi_1 = 2, phi_2 = 0.1, sigma_tau = 1, tol = 1e-20)
chain2 = DP_LMM_BART(2000, 3000, y, X, z, subject_id, v = 3, k, phi_1 = 2, phi_2 = 0.1, sigma_tau = 1, tol = 1e-20)
chain3 = DP_LMM_BART(2000, 3000, y, X, z, subject_id, v = 3, k, phi_1 = 2, phi_2 = 0.1, sigma_tau = 1, tol = 1e-20)

tree = chain1$tree[[length(chain1$tree)]]
B = matrix(chain1$B[length(chain1$tree),], nrow = dim(X), byrow = TRUE)
sigma = chain1$sigma[length(chain1$tree)]




proposed_distribution = function(predict_y, X, z, subject_id, tree, sigma, B, log = TRUE){
  if(is.null(dim(z))){d = 1} else{d = dim(z)[2]}
  Z = matrix(1, nrow = n, ncol = d + 1)
  Z[, 2:(d+1)] = as.matrix(z)
  subject_to_B = seq(n_subject)
  names(subject_to_B) = unique(subject_id)
  
  
  mixed_effect <- foreach(j = seq(length(y)), .combine = rbind) %dopar% {
    sub = subject_id[j]
    b_sub = subject_to_B[as.character(sub)]
    Z[j,] %*% B[b_sub,]
  }
  mean = pwbart(X, tree)[1,] + mixed_effect
  
  if(log)
      sum(dnorm(predict_y, mean, sigma, log))
  else
    prod(dnorm(predict_y, mean, sigma, FALSE))
}

X_new = cbind(X, y)
R_miss = matrix(rbernoulli((dim(X_new)[1] * (dim(X_new)[2] - 1)), 0.2), nrow = dim(X_new)[1])
R_miss = cbind(FALSE, R_miss)
X_new[R_miss] = NA



X_train = matrix(X_new[,1], ncol = 1)
y_train = matrix(X_new[,2], ncol = 1)
miss_marker = is.na(y_train)
y_train = LOCF(y_train)



chain1 = DP_LMM_BART(5, 1, y_train, X_train, z, subject_id, v = 3, k, phi_1 = 2, phi_2 = 0.1, sigma_tau = 1, tol = 1e-20)
new_y_train = chain1$y_predict[1,]
new_y_train[!miss_marker] = y_train[!miss_marker]



# check sigma^2
combinedchains = mcmc.list(mcmc(chain1$sigma), mcmc(chain2$sigma))
gelman.diag(combinedchains)



# check Covariance matrix
for (i in 1:dim(chain1$Sigma)[2]) {
  combinedchains = mcmc.list(mcmc(chain1$Sigma[, i]), mcmc(chain2$Sigma[, i]))
  if(gelman.diag(combinedchains)[[1]][1] > 1.1)
    print(paste(i, gelman.diag(combinedchains)[[1]][1]))
}

# check mean from CDP
for (i in 1:dim(chain1$location)[2]) {
  combinedchains = mcmc.list(mcmc(chain1$location[,i]), mcmc(chain2$location[,i]))
  if(gelman.diag(combinedchains)[[1]][1] > 1.1)
    print(paste(i, gelman.diag(combinedchains)[[1]][1]))
}



# check random effect
for (i in 1:dim(chain1$B)[2]) {
  combinedchains = mcmc.list(mcmc(chain1$B[, i]), mcmc(chain2$B[, i]))
  if(gelman.diag(combinedchains)[[1]][1] > 1.1)
    print(paste(i, gelman.diag(combinedchains)[[1]][1]))
}
