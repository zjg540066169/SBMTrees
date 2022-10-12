library(MASS)
library(BART, quietly = TRUE)
library(invgamma)
library(MCMCpack)
library(mvtnorm)
library(svMisc)
library(parallel) # one of the core R packages
library(doParallel)
suppressMessages(library(BART))
source("data_generate.R")


sink("./myfile.log", append=TRUE, split=TRUE)

nburn = 2
npost = 2

nCores <- detectCores()  # to set manually 
registerDoParallel(nCores) 

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

set.seed(111)
data <- data_generation(sigma_1  = 1, tau_1 = 2, delta_1 = 1, sigma_2  = 1, tau_2 = 2, delta_2 = 1, n_total = 5000, number_of_time_points = 10)


v = 3
k = 1

phi_1 = 2
phi_2 = 0.1
phi_1_2 = 2
phi_2_2 = 0.1

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


DP_MIX_BART = function(nburn, npost, y, X, z, subject_id, v = 3, k, phi_1 = 2, phi_2 = 0.1, phi_1_2 = 2, phi_2_2 = 0.1, sigma_tau = 1, df = NULL, inverse_wishart_matrix = NULL, Covariance_tau = NULL, Mu_tau = NULL, tol = 1e-20){
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
  initial_M_2 = rgamma(1, phi_1_2, phi_2_2)
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
  
  
  M_2 = initial_M_2
  dp_2 = CDP(function(n, p) {
    if (p == 1)
      rnorm(n, 0, 1)
    else{
      rmvnorm(n, rep(0, p))
    }
  }, initial_M_2, n = n, dp = 1)
  initial_location_2 = dp_2$theta
  unique_location_2 = initial_location_2[unique(names(initial_location_2))]
  
  
  
  
  # posterior sampling
  post_y_hat = array(NA, dim = c(npost, length(y)))
  post_M = array(NA, dim = npost)
  post_sigma = array(NA, dim = npost)
  post_Sigma = array(NA, dim = c(npost, (d + 1) * (d + 1)))
  post_B = array(NA, dim = c(npost, n_subject * (d + 1)))
  post_location = array(NA, dim = c(npost, n_subject * (d + 1)))
  post_unique_location = matrix(0, nrow = (npost), ncol = length(initial_unique_location))
  
  
  post_M_2 = array(NA, dim = npost)
  post_location_2 = array(NA, dim = c(npost, n))
  post_unique_location_2 = matrix(0, nrow = (npost), ncol = length(unique_location_2))
  
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
  
  
  M_2 = initial_M_2
  initial_pi_2 = dp_2$pi
  initial_h_2 = dp_2$h
  initial_sigma = rinvgamma(1, v / 2, k * v / 2)
  initial_cluster_2 = as.numeric(names(initial_location_2))
  
  cluster_map_2 = as.numeric(names(unique_location_2))
  names(cluster_map_2) = unique_location_2
  colnames(post_unique_location_2) = names(unique_location_2)
  
  
  sigma = initial_sigma
  location_2 = initial_location_2
  pi_2 = initial_pi_2
  names(pi_2) = names(initial_location_2)
  pi_2 = t(rowsum(pi_2, group = names(pi_2)))[1,]
  
  
  
  for(kl in 1:(nburn + npost)){
    # burn in
    
    
    mixed_effect <- foreach(j = seq(length(y)), .combine = rbind) %dopar% {
      sub = subject_id[j]
      b_sub = subject_to_B[as.character(sub)]
      Z[j,] %*% B[b_sub,]
    }
    
    
    tree = wbart(X, y - mixed_effect - location_2, ndpost = 1, printevery = 1000, nskip = 10) # get tree
    
    y_hat = tree$yhat.train.mean
    sigma = tree$sigma[11]
    
    
    
    
    unique_location_2 = unique_location_2[names(pi_2)] 
    y_residue = y - y_hat - mixed_effect
    y_red = t(sapply(y_residue, function(x) x - unique_location_2))
    prob_2 = pnorm(y_red / sigma)
    prob_2[prob_2 == 0] = tol
    sample_frame_2 =  t(apply(prob_2, 1, function(x) x * pi_2[colnames(prob_2)]))
    
    cluster_2 = apply(sample_frame_2, 1, function(x) sample(as.numeric(names(x)), size = 1, replace = TRUE, prob = x)) 
    location_2 = unique_location_2[as.character(cluster_2)]
    nk_2 = table(location_2)
    names(nk_2) = cluster_map_2[names(nk_2)]
    nk_2[setdiff(names(unique_location_2), names(nk_2))] = 0
    beta_2 = matrix(0, length(nk_2) - 1, 2)
    beta_2[,1] = nk_2[1:(length(nk_2) - 1)] + 1
    beta_2[,2] = M_2 + (rev(cumsum(rev(nk_2))) - nk_2)[1:(length(nk_2) - 1)]
    Vh_2 = c(rbeta(length(nk_2) - 1, beta_2[,1], beta_2[,2]), 1)
    pi_2 = numeric(length(nk_2))
    pi_2[1] = Vh_2[1]
    if(length(nk_2) >= 2)
      pi_2[2:length(nk_2)] = sapply(2:length(nk_2), function(i) Vh_2[i] * prod(1 - Vh_2[1:(i-1)]))
    names(pi_2) = names(nk_2)
    
    unique_location_2[] = NA
    for (j in names(nk_2)) {
      if (nk_2[j] == 0)
        unique_location_2[j] = rnorm(1, 0,  sigma_tau)
      else
        unique_location_2[j] = rnorm(1, sigma_tau^2 / (nk_2[j]*sigma_tau^2 + sigma^2) * sum((y - y_hat - mixed_effect)[cluster_2 == as.numeric(j)]),  sqrt(sigma_tau^2 * sigma^2 / (nk_2[j]*sigma_tau^2 + sigma^2)))
    }
    unique_location_2 = unique_location_2 - sum(unique_location_2[names(unique_location_2)] * pi_2[names(unique_location_2)])
    cluster_map_2 = as.numeric(names(unique_location_2))
    names(cluster_map_2) = unique_location_2
    location__2= sapply(names(location_2), function(x) unique_location_2[x])
    names(location__2) = names(location_2)
    location_2 = location__2
    
    
    
    
    
    
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
      post_y_hat[kl - nburn, ] = y_hat
      #post_unique_location[kl - nburn,] = unique_location
      post_location[kl - nburn,] = location
      
      #post_unique_location_2[kl - nburn,] = unique_location_2[colnames(post_location_2)]
      #post_M_2[kl - nburn] = M_2
      post_location_2[kl - nburn, ] = location_2
    }
    #progress(kl, npost + nburn)
    cat(kl,"/", npost + nburn, "\n")
  }
  return(list(Sigma = post_Sigma, sigma = post_sigma, B = post_B, location = post_location, location_2 = post_location_2))
}


chain1 = DP_MIX_BART(nburn, npost, y, X, z, subject_id, v = 3, k, phi_1 = 2, phi_2 = 0.1, sigma_tau = 1, tol = 1e-20)
chain2 = DP_MIX_BART(nburn, npost, y, X, z, subject_id, v = 3, k, phi_1 = 2, phi_2 = 0.1, sigma_tau = 1, tol = 1e-20)
chain3 = DP_MIX_BART(nburn, npost, y, X, z, subject_id, v = 3, k, phi_1 = 2, phi_2 = 0.1, sigma_tau = 1, tol = 1e-20)


# check sigma^2
cat("check convergence for sigma^2 \n")
combinedchains = mcmc.list(mcmc(chain1$sigma), mcmc(chain2$sigma), mcmc(chain3$sigma))
cat(gelman.diag(combinedchains)[[1]], "\n")



# check Covariance matrix
cat("check convergence for Covariance matrix \n")
for (i in 1:dim(chain1$Sigma)[2]) {
  combinedchains = mcmc.list(mcmc(chain1$Sigma[, i]), mcmc(chain2$Sigma[, i]), mcmc(chain3$Sigma[, i]))
  if(gelman.diag(combinedchains)[[1]][1] > 1.1)
    cat(paste(i, gelman.diag(combinedchains)[[1]][1]), "\n")
}

# check mean from CDP
cat("check convergence for CDP location for random effect \n")
for (i in 1:dim(chain1$location)[2]) {
  combinedchains = mcmc.list(mcmc(chain1$location[,i]), mcmc(chain2$location[,i]), mcmc(chain3$location[,i]))
  if(gelman.diag(combinedchains)[[1]][1] > 1.1)
    cat(paste(i, gelman.diag(combinedchains)[[1]][1]), "\n")
    #print(paste(i, gelman.diag(combinedchains)[[1]][1]))
}



# check mean from CDP
cat("check convergence for CDP location for random noise")
for (i in 1:dim(chain1$location_2)[2]) {
  combinedchains = mcmc.list(mcmc(chain1$location_2[,i]), mcmc(chain2$location_2[,i]), mcmc(chain3$location_2[,i]))
  if(gelman.diag(combinedchains)[[1]][1] > 1.1)
    cat(paste(i, gelman.diag(combinedchains)[[1]][1]), "\n")
}




# check random effect
cat("check convergence for random effect")
for (i in 1:dim(chain1$B)[2]) {
  combinedchains = mcmc.list(mcmc(chain1$B[, i]), mcmc(chain2$B[, i]), mcmc(chain3$B[, i]))
  if(gelman.diag(combinedchains)[[1]][1] > 1.1)
    cat(paste(i, gelman.diag(combinedchains)[[1]][1]), "\n")
}
