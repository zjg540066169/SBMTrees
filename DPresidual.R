library(MASS)
library(BART, quietly = TRUE)
library(invgamma)
library(MCMCpack)
library(mvtnorm)
library(svMisc)
suppressMessages(library(BART))
source("data_generate.R")

set.seed(111)
data <- data_generation(sigma_1  = 1, tau_1 = 2, delta_1 = 1, sigma_2  = 1, tau_2 = 2, delta_2 = 1, n_total = 5000, number_of_time_points = 10)



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

hist(CDP(G_0, 1000))



y = data$x_2
X = data[,c(1,2)]
z = data[, "time"]
p = dim(X)[2]
n = dim(X)[1]
subject_id = data[,"subject_id"]
subject_id = subject_id + 1

v = 3
k = 1

phi_1 = 2
phi_2 = 0.1
sigma_tau = 1


DP_BART = function(nburn, npost, y, X, z, subject_id, v = 3, k, phi_1 = 2, phi_2 = 0.1, sigma_tau = 1, tol = 1e-20){
    if(is.null(dim(z))){d = 1} else{d = dim(z)[2]}

    df = d + 1
    inverse_wishart_matrix = diag(rep(1, d + 1))
    n = length(y)
    p = dim(X)[2]
    n_subject = length(unique(subject_id))
    
    Z = matrix(1, nrow = n, ncol = d + 1)
    Z[, 2:(d+1)] = as.matrix(z)
    
    post_M = array(NA, dim = npost)
    post_sigma = array(NA, dim = npost)
    post_Sigma = array(NA, dim = c(npost, (d + 1) * (d + 1)))
    post_B = array(NA, dim = c(npost, n * (d + 1)))
    post_tau = array(NA, dim = c(npost, n))
    post_location = matrix(0, nrow = (npost), ncol = length(unique_location))
    post_y_hat = array(NA, dim = c(npost, length(y)))
    
    initial_M = rgamma(1, phi_1, phi_2)
    M = initial_M
    dp = CDP(function(n, p) {
      if (p == 1)
        rnorm(n, 0, 1)
      else{
        rmvnorm(n, rep(0, p))
      }
    }, M, n = n)
    initial_location = dp$theta
    unique_location = initial_location[unique(names(initial_location))]
    
    
    
    initial_pi = dp$pi
    initial_h = dp$h
    initial_sigma = rinvgamma(1, v / 2, k * v / 2)
    initial_cluster = as.numeric(names(initial_location))
    
    cluster_map = as.numeric(names(unique_location))
    names(cluster_map) = unique_location
    colnames(post_location) = names(unique_location)
    
    
    sigma = initial_sigma
    location = initial_location
    pi = initial_pi
    names(pi) = names(initial_location)
    pi = t(rowsum(pi, group = names(pi)))[1,]
    
    Covariance = riwish(df, inverse_wishart_matrix)
    a = rinvgamma(n_subject, 1, 1)
    B = mvrnorm(n_subject, rep(0, d + 1), Covariance)
    
    subject_to_B = seq(n_subject)
    names(subject_to_B) = unique(subject_id)
    
    mixed_effect = sapply(seq(dim(X)[1]), function(x) {
      sub = subject_id[x]
      b_sub = subject_to_B[as.character(sub)]
      Z[x,] %*% B[b_sub,]
    })
    
    print("Burning:")
    for(i in 1:(nburn)){
      #print(i)
      tree = bart(X, y - location - mixed_effect, ndpost = 1, verbose = F, nskip = 10) # get tree
      y_hat = tree$yhat.train.mean
      unique_location = unique_location[names(pi)] 
      y_residue = y - y_hat - mixed_effect
      y_red = t(sapply(y_residue, function(x) x - unique_location))
      prob = pnorm(y_red / sigma)
      prob[prob == 0] = tol
      sample_frame =  t(apply(prob,1, function(x) x * pi[colnames(prob)]))
      
      cluster = apply(sample_frame, 1, function(x) sample(as.numeric(names(x)), size = 1, replace = TRUE, prob = x)) 
      location = unique_location[as.character(cluster)]
      nk = table(location)
      names(nk) = cluster_map[names(nk)]
      nk[setdiff(names(unique_location), names(nk))] = 0
      beta = matrix(0, length(nk) - 1, 2)
      beta[,1] = nk[1:(length(nk) - 1)] + 1
      beta[,2] = M + (rev(cumsum(rev(nk))) - nk)[1:(length(nk) - 1)]
      Vh = c(rbeta(length(nk) - 1, beta[,1], beta[,2]), 1)
      pi = numeric(length(nk))
      pi[1] = Vh[1]
      if(length(nk) >= 2)
          pi[2:length(nk)] = sapply(2:length(nk), function(i) Vh[i] * prod(1 - Vh[1:(i-1)]))
      names(pi) = names(nk)
      
      unique_location[] = NA
      for (j in names(nk)) {
        if (nk[j] == 0)
          unique_location[j] = rnorm(1, 0,  sigma_tau)
        else
          unique_location[j] = rnorm(1, sigma_tau^2 / (nk[j]*sigma_tau^2 + sigma^2) * sum((y - y_hat - mixed_effect)[cluster == as.numeric(j)]),  sqrt(sigma_tau^2 * sigma^2 / (nk[j]*sigma_tau^2 + sigma^2)))
      }
      unique_location = unique_location - sum(unique_location[names(unique_location)] * pi[names(unique_location)])
      cluster_map = as.numeric(names(unique_location))
      names(cluster_map) = unique_location
      location_= sapply(names(location), function(x) unique_location[x])
      names(location_) = names(location)
      location = location_
      
      
      # update B
      
      B = t(sapply(names(subject_to_B), function(x) {
        b_sub = subject_to_B[as.character(x)]
        Zi = Z[subject_id == x, ]
        Ri = (y - y_hat - location)[subject_id == x]
        var = solve(solve(Covariance) + t(Zi) %*% (Zi) / sigma)
        mu = (t(t(Zi) %*% Ri) / sigma)  %*% var 
        mvrnorm(1, mu, var)
      }))
      
      
      # update Covariance matrix
      Covariance = riwish(df + n_subject, inverse_wishart_matrix + t(B) %*% B)
      
      
      mixed_effect = sapply(seq(dim(X)[1]), function(x) {
        sub = subject_id[x]
        b_sub = subject_to_B[as.character(sub)]
        Z[x,] %*% B[b_sub,]
      })
      
      
      # update M
      phi_1_new = phi_1 + length(nk) - 1
      if(length(nk) >= 2)
          phi_2_new = phi_2 - sum(log(1 - Vh[1:(length(nk) - 1)]))
      M = rgamma(1, phi_1_new, phi_2_new)
      
      # update sigma2
      s2 = sum((y - y_hat - location - mixed_effect)^2)
      sigma = rinvgamma(1, (v + n) / 2, (s2 + k * v) / 2)
      
      progress(i, nburn)
    }
    
    print("Sampling:")
    for(i in 1:(npost)){
      #print(i)
      tree = bart(X, y - location - mixed_effect, ndpost = 1, verbose = F, nskip = 10) # get tree
      y_hat = tree$yhat.train.mean
      unique_location = unique_location[names(pi)] 
      y_residue = y - y_hat - mixed_effect
      y_red = t(sapply(y_residue, function(x) x - unique_location))
      prob = pnorm(y_red / sigma)
      prob[prob == 0] = tol
      sample_frame =  t(apply(prob,1, function(x) x * pi[colnames(prob)]))
      
      cluster = apply(sample_frame, 1, function(x) sample(as.numeric(names(x)), size = 1, replace = TRUE, prob = x)) 
      location = unique_location[as.character(cluster)]
      nk = table(location)
      names(nk) = cluster_map[names(nk)]
      nk[setdiff(names(unique_location), names(nk))] = 0
      beta = matrix(0, length(nk) - 1, 2)
      beta[,1] = nk[1:(length(nk) - 1)] + 1
      beta[,2] = M + (rev(cumsum(rev(nk))) - nk)[1:(length(nk) - 1)]
      Vh = c(rbeta(length(nk) - 1, beta[,1], beta[,2]), 1)
      pi = numeric(length(nk))
      pi[1] = Vh[1]
      if(length(nk) >= 2)
        pi[2:length(nk)] = sapply(2:length(nk), function(i) Vh[i] * prod(1 - Vh[1:(i-1)]))
      names(pi) = names(nk)
      
      unique_location = numeric()
      for (j in names(nk)) {
        if (nk[j] == 0)
          unique_location[j] = rnorm(1, 0,  sigma_tau)
        else
          unique_location[j] = rnorm(1, sigma_tau^2 / (nk[j]*sigma_tau^2 + sigma^2) * sum((y - y_hat - mixed_effect)[cluster == as.numeric(j)]),  sqrt(sigma_tau^2 * sigma^2 / (nk[j]*sigma_tau^2 + sigma^2)))
      }
      unique_location = unique_location - sum(unique_location[names(unique_location)] * pi[names(unique_location)])
      cluster_map = as.numeric(names(unique_location))
      names(cluster_map) = unique_location
      location_= sapply(names(location), function(x) unique_location[x])
      names(location_) = names(location)
      location = location_
      
      
      # update M
      phi_1_new = phi_1 + length(nk) - 1
      if(length(nk) >= 2)
        phi_2_new = phi_2 - sum(log(1 - Vh[1:(length(nk) - 1)]))
      M = rgamma(1, phi_1_new, phi_2_new)
      
      
      # update B
      
      B = t(sapply(names(subject_to_B), function(x) {
        b_sub = subject_to_B[as.character(x)]
        Zi = Z[subject_id == x, ]
        Ri = (y - y_hat - location)[subject_id == x]
        var = solve(solve(Covariance) + t(Zi) %*% (Zi) / sigma)
        mu = (t(t(Zi) %*% Ri) / sigma)  %*% var 
        mvrnorm(1, mu, var)
      }))
      
      
      
      
      
      # update Covariance matrix
      Covariance = riwish(df + n_subject, inverse_wishart_matrix + t(B) %*% B)
      
      
      mixed_effect = sapply(seq(dim(X)[1]), function(x) {
        sub = subject_id[x]
        b_sub = subject_to_B[as.character(sub)]
        Z[x,] %*% B[b_sub,]
      })
      
      

      
      # update sigma2
      s2 = sum((y - y_hat - location - mixed_effect)^2)
      sigma = rinvgamma(1, (v + n) / 2, (s2 + k * v) / 2)
      
      
      # colllect
      post_location[i,] = unique_location[colnames(post_location)]
      post_M[i] = M
      post_sigma[i] = sigma
      post_Sigma[i,] = Covariance
      post_B[i,] = B
      post_tau[i, ] = location
      post_y_hat[i, ] = y_hat
      progress(i, npost)
    }   
    return(list(post_tau = post_tau, post_Sigma = post_Sigma, post_B = post_B, post_location =  post_location, post_M = post_M, post_sigma = post_sigma))
}

npost = 10
nburn = 10

chain1 = DP_BART(nburn, npost, y, X, z, subject_id, v, k, phi_1 = 2, phi_2 = 0.1, sigma_tau = 1)
chain2 = DP_BART(nburn, npost, y, X, z, subject_id, v, k, phi_1 = 2, phi_2 = 0.1, sigma_tau = 1)
chain3 = DP_BART(nburn, npost, y, X, z, subject_id, v, k, phi_1 = 2, phi_2 = 0.1, sigma_tau = 1)

# check the convergence of B
for (i in 1:dim(chain1$post_B)[2]) {
  combinedchains = mcmc.list(mcmc(chain1$post_B[, i]), mcmc(chain2$post_B[, i]), mcmc(chain3$post_B[, i]))
  if(gelman.diag(combinedchains)[[1]][1] > 1.1)
    print(paste(i, gelman.diag(combinedchains)[[1]][1]))
}
matrix(colMeans(chain1$post_B), ncol = dim(z)[2] + 1, byrow = TRUE)

# check the convergence of location
for (i in 1:dim(chain1$post_tau)[2]) {
  combinedchains = mcmc.list(mcmc(chain1$post_tau[, i]), mcmc(chain2$post_tau[, i]), mcmc(chain3$post_tau[, i]))
  if(gelman.diag(combinedchains)[[1]][1] > 1.1)
    print(paste(i, gelman.diag(combinedchains)[[1]][1]))
}
colMeans(chain1$post_tau)

# check the convergence of sigma
combinedchains = mcmc.list(mcmc(chain1$post_sigma), mcmc(chain2$post_sigma), mcmc(chain3$post_sigma))
gelman.diag(combinedchains)
mean(chain1$post_sigma)

# check the convergence of Sigma, covariance matrix
for (i in 1:dim(chain1$post_Sigma)[2]) {
  combinedchains = mcmc.list(mcmc(chain1$post_Sigma[, i]), mcmc(chain2$post_Sigma[, i]), mcmc(chain3$post_Sigma[, i]))
  if(gelman.diag(combinedchains)[[1]][1] > 1.1)
    print(paste(i, gelman.diag(combinedchains)[[1]][1]))
}
matrix(colMeans(chain1$post_Sigma), nrow = dim(z)[2] + 1)

# check the convergence of M
combinedchains = mcmc.list(mcmc(chain1$post_M), mcmc(chain2$post_M), mcmc(chain3$post_M))
gelman.diag(combinedchains)

# check the convergence of location
combinedchains = mcmc.list(mcmc(chain2$post_location), mcmc(chain3$post_location))
gelman.diag(combinedchains)



alcohol1 <- read.table("https://stats.idre.ucla.edu/stat/r/examples/alda/data/alcohol1_pp.txt", header=T, sep=",")
attach(alcohol1)
y = alcohol1$alcuse
X = alcohol1[,c(3, 4, 7, 8, 9)]
z = alcohol1[, c("age_14", "age")]
p = dim(X)[2]
n = dim(X)[1]
subject_id = alcohol1[,"id"]


npost = 3000
nburn = 3000

chain1 = DP_BART(nburn, npost, y, X, z, subject_id, v, k, phi_1 = 2, phi_2 = 0.1, sigma_tau = 1)
chain2 = DP_BART(nburn, npost, y, X, z, subject_id, v, k, phi_1 = 2, phi_2 = 0.1, sigma_tau = 1)
chain3 = DP_BART(nburn, npost, y, X, z, subject_id, v, k, phi_1 = 2, phi_2 = 0.1, sigma_tau = 1)
