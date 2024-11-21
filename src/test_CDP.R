library(Rcpp)
library(RcppParallel)
library(RcppArmadillo)
library(RcppDist)
library(tidyverse)
sourceCpp("./src/BMTrees_MCMC.cpp")
#sourceCpp("./src/bmtrees2.cpp")
sourceCpp("./src/bmtrees.cpp")
a = test(Y, X, Z, subject_id, F, T, F)

sourceCpp("./src/utils.cpp")

create_row_id_to_row(seq(1, 10, 2))





source("./R/nCDP.R")

sourceCpp("./src/bart_model.cpp")

a = bart_train(X, Y, nburn = 100, npost = 1, T)


sourceCpp("./src/DP2.cpp")

a = 1:3
a = matrix_multiply(as.matrix(a), t(as.matrix(a)))

a = rtgamma(100, 2, 1 / 10, 1, 2)
hbmem::rtgamma(100, 2, 2, 1, 2)

DP_sampler(100, list("p" = 3, "sigma" = diag(3) * 0.1, "mu" = c(3, 2, -1)))

n <- 100
mat1 <- matrix(runif(n^2), nrow = n, ncol = n)
mat2 <- matrix(runif(n^2), nrow = n, ncol = n)

mul = function(mat1, mat2){
  mat1 %*% mat2
}


mb <- microbenchmark(
  matrix_multiply(mat1, mat2),
  mul(mat1, mat2),
  times = 100  # Number of repetitions
)

ID = sort(rep(-5:5, 5))
X = cbind(c(sapply(ID, function(x) rnorm(20, x, 0.1))), c(sapply(ID, function(x) rnorm(20, x, 0.1))), c(sapply(ID, function(x) rnorm(20, x, 0.1))))
ID = as.character(as.integer(as.factor(ID)))


#X = c(rnorm(100, 20, 1), rnorm(100, -20, 1), rnorm(100, 10, 1) , rnorm(100, -10, 1))
#ID = c(rep("1", 100), rep("2", 100), rep("3", 100), rep("4", 100))


#tau = DP(parameters = list("p" = 1), M = 2,  N_truncated = 300, N_sample = length(X), CDP = FALSE)
 

#tau = update_nCDP(as.matrix(X), tau, 0.01, tol = 1e-80, CDP = T)

#tau2["sigma"] = 0.01
#tau2 = update_nDP_normal(as.matrix(X), tau2, tol = 1e-80)


#cbind(do.call(rbind, lapply(tau$y, function(x) x$y)), do.call(rbind, lapply(tau2$y, function(x) x$y)))
start = Sys.time()
tau2 = DP(parameters = list("p" = 3, "Sigma" = cov(X)), M = 2, N_truncated = sqrt(dim(X)[1]), N_sample = dim(X)[1], CDP = T)
tau2 = DP(parameters = list("p" = 1, "sigma" = sd(X)), M = 2, N_truncated = sqrt(dim(X)[1]), N_sample = dim(X)[1], CDP = F)
sample_post2 = matrix(NA, nrow = 3000, ncol = length(tau2$samples))
cluster_post2 = matrix(NA, nrow = 3000, ncol = length(tau2$y))
for(i in 1:6000){
  print(paste0(i, ":", tau2$M, " ", length(table(tau2$cluster))))
  #tau2[["sigma"]] = 0.1#^2 * diag(3)
  tau2[["Sigma"]] = 0.1^2 * diag(3)
  tau2 = update_DP_normal(as.matrix(X), tau2)
  if(i > 3000){
    sample_post2[i - 3000,] = tau2$samples
    cluster_post2[i - 3000,] = tau2$y
  }
  
  print(table(tau2$cluster))
}


tau = DP(parameters = list("p" = 3, "Sigma" = cov(X)), M = 2, N_truncated = sqrt(dim(X)[1]), N_sample = dim(X)[1], CDP = T)
tau = DP(parameters = list("p" = 1, "sigma" = sd(X)), M = 2, N_truncated = sqrt(dim(X)[1]), N_sample = dim(X)[1], CDP = F)
sample_post = matrix(NA, nrow = 3000, ncol = length(tau$samples))
cluster_post = matrix(NA, nrow = 3000, ncol = length(tau$y))
for(i in 1:6000){
  print(paste0(i, ":", tau$M, " ", length(table(tau$cluster))))
  tau[["Sigma"]] = 0.1^2 * diag(3)
  #tau[["sigma"]] = 0.1#^2 * diag(3)
  tau = update_DP_normal(as.matrix(X), tau)
  if(i > 3000){
    sample_post[i - 3000,] = tau$samples
    cluster_post[i - 3000,] = tau$y
  }
  print(table(tau$cluster))
}

end = Sys.time()


combinedchains = mcmc.list(mcmc(sample_post), mcmc(sample_post2)) 

a = gelman.diag(combinedchains, multivariate = T)

library(dirichletprocess)
dp <- DirichletProcessGaussian(X, g0Priors = c(0, 1, 0.01, 1))
dp<-Fit(dp, 1000)

samples <- list()
for(s in 1:1000){
  print(paste0(i, ":", dp$alpha, " ", length(table(dp$clusterLabels))))
  dp = ClusterComponentUpdate(dp)
  dp = ClusterParameterUpdate(dp)
  dp = UpdateAlpha(dp)
  print(table(dp$clusterLabels))
}


start2 = Sys.time()
tau = nestedCDP(G_0, parameters = list("p" = 1), cluster_id = ID, M_alpha = 2, M_beta = 2, L = 2, K = 51, dp = 1,  CDP = F)
for(i in 1:300){
  print(paste0(i, ":", tau$M_alpha, " ", tau$M_beta))
  #tau[["sigma"]] = 0.01
  tau = update_nCDP(as.matrix(X), tau, 0.01, tol = 1e-80, CDP = F)
  print(table(tau$cluster))
}
end2 = Sys.time()




library(BART3) 

X1 = rnorm(1000)
X2 = X1 + X1^2 + rgamma(1000, 1,1) * (-1)^rbinom(1000,1, prob = 0.5)
X3 = X1 + X2 + 10 - X1 * X2 + X2^2+ rgamma(1000, 10,10) * (-1)^rbinom(1000,1, prob = 0.5)

sourceCpp("./src/BMTrees_MCMC.cpp")

library(mvtnorm)
dat = simulation(100, nonlinear = T, nonrandeff = T, nonresidual = T, alligned = T)
subject_id = dat$subject_id

nburn = 3000
npost = 3000
Y = unlist(dat$Y_O)
X = as.matrix(dat$X_O)

#X = as.matrix(rep(1, length(Y)))
Z = dat$Z

sourceCpp("./src/update_Covariance.cpp")
subject_to_B = test(Y, X = X, Z = Z, subject_id = subject_id, F)

mu = matrix(0, nrow = length(unique(subject_id)), ncol = ncol(Z) + 1)

covariance = bdiag(VarCorr(slmm))
covariance = as.matrix(covariance)


coe = as.matrix(coef(slmm)$subject_id)[,c(3,1,2)]
coe[,1] = coe[,1] - mean(Y)

a = sapply(1:1000, function(i)
  as.numeric(update_Covariance(coe, mu, as.matrix(cov(coe)) , 0, 100))
  )
rowMeans(a)
covariance




update_Covariance(as.matrix(coef(slmm)$subject_id)=, mu, diag(0), 4, 100)

a = t(a)
a[,1] = a[,1] + mean(Y)
summary(slmm)

a = update_B_CPP(Y - mean(Y), as.matrix(cbind(1, Z)), subject_id, mu, subject_to_B[[1]], covariance, 8.32234)

sourceCpp("./src/bart_model.cpp")
m = bart_train(X, Y, 100, 1000, T)

sourceCpp("./src/bmtrees2.cpp")
chain1 = test(Y, X, Z, subject_id, F, T, T)
sourceCpp("./src/BMTrees_MCMC.cpp")

chain1 = foreach(1:2, .combine = list) %dopar%{
  BMTrees_MCMC(6000, 6000, Y, as.matrix(X), as.matrix(Z), subject_id, FALSE, verbose =  F, NULL, tol = 1e-100, F, F)
}
  
chain5 = foreach(1:2, .combine = list) %dopar%{
  BMTrees_MCMC(6000, 6000, Y, as.matrix(X), as.matrix(Z), subject_id, FALSE, verbose =  F, NULL, tol = 1e-100, T, F)
}


chain3 = foreach(1:2, .combine = list) %dopar%{
  BMTrees_MCMC(6000, 6000, Y, as.matrix(X), as.matrix(Z), subject_id, FALSE, verbose =  F, NULL, tol = 1e-100, F, T)
}

chain4 = foreach(1:2, .combine = list) %dopar%{
  BMTrees_MCMC(6000, 6000, Y, as.matrix(X), as.matrix(Z), subject_id, FALSE, verbose =  F, NULL, tol = 1e-100, T, T)
}

  
library(doParallel)
library(foreach)
registerDoParallel(2)




chain4 = BMTrees_MCMC(6000, 3000, Y, X, Z, subject_id, FALSE, verbose = T, NULL, tol = 1e-100, F,T)
chain3 = BMTrees_MCMC2(nburn, npost, Y, X, Z, subject_id, FALSE, verbose =  F, 123123, tol = 1e-100, T, T)
chain4 = BMTrees_MCMC(nburn, npost, Y, X, Z, subject_id, FALSE, verbose =  F, 121233, tol = 1e-100, T, T)
sourceCpp("./src/DP2.cpp")

update_DP_normal(chain3$B, chain3$B_tau)

combinedchains = mcmc.list(mcmc(mse_bart), mcmc(mse_bart2))
a = gelman.diag(combinedchains, multivariate = F)



chain5 = BMTrees_MCMC2(nburn, npost, Y, X, Z, subject_id, FALSE, verbose =  F, 121212333, tol = 1e-100, T, T)
chain6 = BMTrees_MCMC2(nburn, npost, Y, X, Z, subject_id, train = rbinom(length(Y),1 , prob = 0.9), FALSE, verbose =  F, 121, tol = 1e-100, T, T)

chain2 = BMTrees_MCMC(nburn, npost, Y, X, Z, subject_id, FALSE, verbose =  F, 1212334, tol = 1e-100, F, T)

chain3 = BMTrees_MCMC(nburn, npost, Y, X, Z, subject_id, FALSE, verbose =  F, 123485, tol = 1e-100, F, T)

library(doParallel)
library(foreach)
registerDoParallel(3)

chain = foreach(i = 19:21, .combine = list, .multicombine = T) %dopar%{
  BMTrees_MCMC(nburn, npost, Y, X, Z, subject_id, FALSE, verbose =  F, NULL, tol = 1e-100, T, T)
}

Z_scale = scale(Z)
X_scale = scale(X)
Y_scale = scale(Y)
library(lme4)
slmm = lmer(Y ~ X_scale + (Z|subject_id), REML = F)
bdiag(VarCorr(slmm))


slmm_Y = predict(slmm)

mean((slmm_Y - Y))


mf = MixRF(Y = Y, X = X, random='(Z1+Z2|subject_id)',data = as.data.frame(cbind(X, Y, Z1 = Z[,1], Z2 = Z[,2], subject_id)),initialRandomEffects=0, ErrorTolerance=0.01, MaxIterations=100)
mf_Y = predict.MixRF(mf, newdata = as.data.frame(cbind(X = X, Z1 = Z[,1], Z2 = Z[,2], subject_id = subject_id)), EstimateRE=F)
mean((mf_Y - Y)^2)
mean(mf_Y - Y)

mse_bart = apply(chain1[[1]]$post_y_predict, 1, function(x) mean((x - Y)^2))
mse_bart2 = apply(chain1[[2]]$post_y_predict, 1, function(x) mean((x - Y)^2))

mse_residual_cdp = apply(chain2[[1]]$post_y_predict_new, 1, function(x) mean((x - Y)^2))
mse_residual_cdp2 = apply(chain2[[2]]$post_y_predict_new, 1, function(x) mean((x - Y)^2))

mse_residual_cdp3 = mean((colMeans(chain5[[1]]$post_x_hat) - Y))
mse_residual_cdp4 = mean((colMeans(chain5[[2]]$post_x_hat + chain5[[2]]$post_tau_samples) - Y)^2)



mse_ranef_cdp = apply(chain3[[1]]$post_y_predict, 1, function(x) mean((x - Y)^2))
mse_ranef_cdp2 = apply(chain3[[2]]$post_y_predict, 1, function(x) mean((x - Y)^2))

mse_dp_double = apply(chain4[[1]]$post_y_predict, 1, function(x) mean((x - Y)^2))
mse_dp_double2 = apply(chain4[[2]]$post_y_predict, 1, function(x) mean((x - Y)^2))

summary(cbind(mse_bart, mse_bart2, mse_residual_cdp, mse_residual_cdp2, mse_ranef_cdp, mse_ranef_cdp2, mse_dp_double, mse_dp_double2))


mse_cdp_double = apply(chain5$post_y_predict, 1, function(x) mean((x - Y)^2))
mse_cdp_non = apply(chain3$post_y_predict, 1, function(x) mean((x - Y)^2))

saveRDS(chain, "double_DP.rds")
chain = readRDS("double_DP.rds")

combinedchains = mcmc.list(mcmc(chain[[1]]$post_tau_samples), mcmc(chain[[2]]$post_tau_samples), mcmc(chain[[3]]$post_tau_samples))
a = gelman.diag(combinedchains, multivariate = T)
print(gelman.diag(combinedchains, multivariate = F))



combinedchains = mcmc.list(mcmc(colMeans(apply(chain[[1]]$post_y_predict, 1, function(y) (y - Y)^2))), mcmc(colMeans(apply(chain[[2]]$post_y_predict, 1, function(y) (y - Y)^2))), mcmc(colMeans(apply(chain[[3]]$post_y_predict, 1, function(y) (y - Y)^2))))







chain4 = BMTrees_MCMC(nburn, npost, Y, X, Z, subject_id, FALSE, verbose =  F, 123, tol = 1e-100, T, T)

chain5 = BMTrees_MCMC(nburn, npost, Y, X, Z, subject_id, FALSE, verbose =  F, 1212334, tol = 1e-100, T, T)

chain6 = BMTrees_MCMC(nburn, npost, Y, X, Z, subject_id, FALSE, verbose =  F, 123485, tol = 1e-100, T, T)






sourceCpp("./src/DP2.cpp")
sourceCpp("./src/nDP2.cpp")

n_subject = length(unique(subject_id))
M_re = n_subject^(runif(1, 0, 0.5))
B_tau = DP(list("p" = dim(Z)[2]), M_re, sqrt(n_subject), n_subject, T)
B_tau_samples = B_tau["samples"]


M_alpha = n_subject^(runif(1, 0, 0.5))
M_beta = 5^(runif(1, -1, 2))
tau = nDP(list("p" = 1), subject_id, M_alpha, M_beta, 5, 30, T)
tau_samples = tau["samples"]
Z = dat$Z
Z = apply(Z, 2, function(x) (x - min(x)) / (max(x) - min(x)))
Z = cbind(1, Z)
dz = dim(Z)[2]

covariance = (Y - mean(Y)) %*% t(Y - mean(Y))
belongs = matrix(NA, length(Y), length(Y))
Z_2 = matrix(NA, length(Y), length(Y))
for(i in 1:length(Y)){
  for(j in 1:i){
    Z_2[i, j] = Z[i,] %*% Z[j,]
    if(B_tau$cluster[subject_id[i]] == B_tau$cluster[subject_id[j]]){
      belongs[i, j] = 0
    }else{
      belongs[i, j] = 2
    }
    
    if(tau$cluster[i] == tau$cluster[j] & tau$within_cluster[i] == tau$within_cluster[j]){
      belongs[i, j] = belongs[i, j] + 1
    }
    if(tau$cluster[i] == tau$cluster[j] & tau$within_cluster[i] != tau$within_cluster[j]){
      belongs[i, j] = belongs[i, j] + 2
    }
    if(tau$cluster[i] != tau$cluster[j]){
      belongs[i, j] = belongs[i, j] + 3
    }
  }
}

Y_C = covariance[!is.na(belongs)]
X_b = belongs[!is.na(belongs)]
X_Z2 = Z_2[!is.na(belongs)]
X = cbind(X_b, X_Z2) %>% 
  as_tibble() %>% 
  mutate(
    tau = case_when(
      X_b == 1 ~ 1,
      X_b == 2 ~ 1,
      X_b == 3 ~ 1,
      X_b == 4 ~ 1 / (1 + M_re),
      X_b == 5 ~ 1 / (1 + M_re),
      X_b == 6 ~ 1 / (1 + M_re)
    ),
    tau = tau * X_Z2,
    k = case_when(
        X_b == 1 ~ 1,
        X_b == 2 ~ 1,
        X_b == 3 ~ 1,
        X_b == 4 ~ 0,
        X_b == 5 ~ 0,
        X_b == 6 ~ 0
      ),
    k = k * X_Z2,
    kv = case_when(
      X_b == 1 ~ 1,
      X_b == 2 ~ 0,
      X_b == 3 ~ 0,
      X_b == 4 ~ 1,
      X_b == 5 ~ 0,
      X_b == 6 ~ 0
    ),
    sigma2T = case_when(
      X_b == 1 ~ 1,
      X_b == 2 ~ 1 / (1 + M_beta),
      X_b == 3 ~ 1 / (1 + M_beta) * 1 / (1 + M_alpha),
      X_b == 4 ~ 1,
      X_b == 5 ~ 1 / (1 + M_beta),
      X_b == 6 ~ 1 / (1 + M_beta) * 1 / (1 + M_alpha),
    )
  )
sigma2T = X$sigma2T
tau = X$tau

library(colf)

colf_nlxb(Y_C - sigma2T - tau ~ 0 + k + kv, data = X, lower = c(0, 0))








#tree = gbart(cbind(X1, X2), X3, ndpost = 10, nskip = 10)
#tau = DP(parameters = list("p" = 1), M = 2,  N_truncated = sqrt(length(X3)), N_sample = length(X3), CDP = TRUE)
tree = gbart(as.matrix(X), unlist(Y), ndpost = 10, nskip = 10)
#tau = update_nDP_normal(G_0, parameters = list("p" = 1), cluster_id = ID, M_alpha = 2, M_beta = 2, L = 5, K = 10, dp = 1,  CDP = T)
tau = nDP(parameters = list("p" = 1, "sigma" = 10), ID, M_alpha = 2, M_beta = 2, L = 5, K = 10, CDP = T, seed = 1234)
for(i in 1:6000){
  print(paste0(i, ":", tau$M_alpha, " ", tau$M_beta))
  print(table(tau$cluster))
  print(table(tau$within_cluster))
  tree = gbart(as.matrix(X), Y - tau$samples, ndpost = 1, nskip = 1, verbose = F)
  tau[["sigma"]] = mean(tree$sigma)
  tau = update_nDP_normal(as.matrix(Y - tree$yhat.train.mean), tau, tol = 1e-80)
}

tree_predict = array(NA, dim = c(6000, 1000))
residual = array(NA, dim = c(6000, 1000))
sigma = array(NA, dim = c(6000, 1))
for(i in 1:6000){
  print(paste0(i, ":", tau$M_alpha, " ", tau$M_beta))
  tree = gbart(as.matrix(X), Y - tau$samples, ndpost = 1, nskip = 1, verbose = F)
  tree_predict[i, ] = tree$yhat.train.mean
  tau[["sigma"]] = mean(tree$sigma)
  sigma[i,1] = mean(tree$sigma)
  tau = update_nDP_normal(as.matrix(Y - tree$yhat.train.mean), tau, tol = 1e-80)
  residual[i,] = tau$samples
}


tree_predict2 = array(NA, dim = c(6000, 1000))
residual2 = array(NA, dim = c(6000, 1000))
sigma2 = array(NA, dim = c(6000, 1))
tree = gbart(as.matrix(X), unlist(Y), ndpost = 10, nskip = 10)
set.seed(12345)
#tau = update_nDP_normal(G_0, parameters = list("p" = 1), cluster_id = ID, M_alpha = 2, M_beta = 2, L = 5, K = 10, dp = 1,  CDP = T)
tau = nDP(parameters = list("p" = 1, "sigma" = 10), ID, M_alpha = 2, M_beta = 2, L = 5, K = 10, CDP = T, seed = 12345)
for(i in 1:6000){
  print(paste0(i, ":", tau$M_alpha, " ", tau$M_beta))
  print(table(tau$cluster))
  print(table(tau$within_cluster))
  tree = gbart(as.matrix(X), Y - tau$samples, ndpost = 1, nskip = 1, verbose = F)
  tau[["sigma"]] = mean(tree$sigma)
  tau = update_nDP_normal(as.matrix(Y - tree$yhat.train.mean), tau, tol = 1e-80)
}

for(i in 1:6000){
  print(paste0(i, ":", tau$M_alpha, " ", tau$M_beta))
  tree = gbart(as.matrix(X), Y - tau$samples, ndpost = 1, nskip = 1, verbose = F)
  tree_predict2[i, ] = tree$yhat.train.mean
  tau[["sigma"]] = mean(tree$sigma)
  sigma2[i,1] = mean(tree$sigma)
  tau = update_nDP_normal(as.matrix(Y - tree$yhat.train.mean), tau, tol = 1e-80)
  residual2[i,] = tau$samples
}

combinedchains = mcmc.list(mcmc(sigma[1:2000]), mcmc(sigma2[1:2000]), mcmc(chain1$post_sigma))
a = gelman.diag(combinedchains, multivariate = FALSE)
print(gelman.diag(combinedchains, multivariate = FALSE))


var(Y) = E(var(Y|mu)) + var(E(Y|mu))
E(sigma) + 


a = test(data$Y, data$X_O, as.matrix(Z), subject_id, binary = FALSE, nCDP_residual = TRUE, CDP_re = TRUE)




X = data$Y
subject_id = data$subject_id
tau = nDP(parameters = list("p" = 1), cluster_id = subject_id, M_alpha = 2, M_beta = 2, L = 55, K = 35, CDP = TRUE)
tau$sigma  = 2
tau = update_nDP_normal(as.matrix(X), tau)

a1 = BMTrees_MCMC(10, 10, data$Y, data$X_O, as.matrix(Z), subject_id, TRUE, TRUE, 12, 1e-40, FALSE, FALSE)


for (i in 1:dim(a1$post_y_predict)[2]) {
  combinedchains = mcmc.list(mcmc(a1$post_y_predict[,i]), mcmc(a2$post_y_predict[,i]), mcmc(a3$post_y_predict[,i]))
  print(gelman.diag(combinedchains, multivariate = FALSE))
}

combinedchains = mcmc.list(mcmc(a1$post_sigma), mcmc(a2$post_sigma), mcmc(a3$post_sigma))


sourceCpp("./src/DP2.cpp")

sourceCpp("./Rcpp/update_BART.cpp")
sourceCpp("./Rcpp/update_B.cpp")
sourceCpp("./Rcpp/update_alpha.cpp")
sourceCpp("./Rcpp/update_Covariance.cpp")
sourceCpp("./Rcpp/create_subject_to_B.cpp")
sourceCpp("./Rcpp/random_effects.cpp")

sourceCpp("./src/DP_LMM_BART2.cpp")

sourceCpp("./Rcpp/sequential_imputation.cpp")

update_alpha2 = function(alpha, B, mu_s, covariance, subject_to_B, seed = 123){
  set.seed(seed)
  for (x in names(subject_to_B)) {
    b_sub = unlist(subject_to_B[as.character(x)])
    #print(b_sub)
    Bi = B[b_sub, ]
    musi = mu_s[b_sub, ]
    alpha[b_sub] = MCMCpack::rinvgamma(1, 1 + (Bi - musi) %*% solve(covariance) %*% (Bi - musi) / 2 )
  }
  return(alpha)
}



path_ = "./simulation_continuous/"  
path_dir = list.dirs(path_, recursive = FALSE)
i = 1
scenario = list.dirs(path_dir[i], recursive = FALSE)
#scenario = scenario[stringr::str_ends(scenario, "FALSE")]
cat(i, "\n")
j = 1
path = paste(scenario[j], "/", sep = "")
if("SLMM1.csv" %in% list.files(path, recursive = FALSE)){
  return()
}
Original_X = read.csv(paste(path, "X_O.csv", sep = ""))
X = read.csv(paste(path, "X_mis.csv", sep = ""))
R = is.na(read.csv(paste(path, "X_mis.csv", sep = "")))
subject_id = read.csv(paste(path, "subject_id.csv", sep = ""))[,1]
Z = read.csv(paste(path, "Z.csv", sep = ""))
colnames(Z) = paste("Z", seq(1,2),sep="")
Y = read.csv(paste(path, "Y.csv", sep = ""))

R = is.na(X)
type = c(0, 0, 0, 0, 0, 1, 1, 0, 0, 0)
Z = as.matrix(Z)
X = as.matrix(X)


Z = scale(Z)

library(lme4)


for(i in 2:dim(R)[2]){
  if(sum(R[,i]) != 0){
    #print(i)
    x_train = as.matrix(X[,1:(i-1)])
    y_train = X[,i]
    if(type[i] == 0){
      lmer_1 <- lmer(y_train ~ x_train + (Z|subject_id))
      y_predict = predict(lmer_1, allow.new.levels = TRUE, newdata = as.data.frame(cbind(x_train, Z, subject_id)))
    }else{
      lmer_1 <- lmer(y_train ~ x_train + (Z|subject_id))
      y_predict = 1 * (predict(lmer_1, allow.new.levels = TRUE, newdata = as.data.frame(cbind(x_train, Z, subject_id))) >= 0.5)
    }
    X[R[,i],i] = y_predict[R[,i]]
  }
}







X = as.matrix(Original_X[, 1:6])
Y = as.numeric(Original_X[,7])
Z = as.matrix(Z)
Z = scale(Z)
subject_id = as.character(subject_id)



b = update_BART(X, unlist(Y))




tau = DP_LMM_BART(chain = 1, nburn = 0, npost = 1, Y, X, Z, subject_id, verbose = FALSE, binary = TRUE, seed =123, last_states = tau)

sourceCpp("./Rcpp/DP_LMM_BART.cpp")


start_time <- Sys.time()
tau2 = DP_LMM_BART(chain = 2, nburn = 10, npost = 10, Y, X, Z, subject_id, seed =1234, nDP_residual = TRUE, DP_re = FALSE)
end_time <- Sys.time()
end_time - start_time

tau3 = DP_LMM_BART(chain = 3, nburn = 6000, npost = 6000, Y, X, Z, subject_id, seed =12345)

#-3.57385 3.86291 -6.81601

subject_to_B = create_subject_to_B(subject_id)

alpha = rep(1, length(unique(subject_id)))
alpha2 = rep(1, length(unique(subject_id)))
B_dp = DP(list(p = 2), M = 2, N_truncated = 100, N_sample = length(unique(subject_id)))







for(i in 1:100){

  B_dp$Sigma = Covariance
  b = update_BART(as.matrix(Original_X) - re, unlist(Y))
  residual = Y - b$tree_pre
  sigma = b$sigma
  
  
  B = update_B(residual, Z, subject_id, B_mu, subject_to_B, Covariance, sigma, alpha, 123)
  
  
  B_dp = update_DP_normal(B, B_dp, seed = 1234)
  B_mu = B_dp$samples
  
  
  Covariance = update_Covariance(B, alpha, B_mu, subject_to_B, diag(dim(Z)[2]), 3, 1000)
  alpha = update_alpha(alpha, B, B_mu, Covariance, subject_to_B)
  
  re = random_effects(as.matrix(Z), subject_id, B, subject_to_B)
  
  print(Covariance)
}












random_effects(as.matrix(Z), subject_id, B, subject_to_B)


subject_to_B = create_subject_to_B(subject_id)
b = update_BART(as.matrix(Original_X), unlist(Y))
re =  unlist(Y) - b$tree_pre
re = unname(re)
sigma = b$sigma

alpha = rep(1, length(unique(subject_id)))

mu = DP(list(p = 2), M = 2, N_truncated = 100, N_sample = length(unique(subject_id)))$samples
inv_covariance = diag(2)

B = update_B(re, as.matrix(Z), subject_id, mu, subject_to_B, inv_covariance, sigma, alpha)
summary(B)

#alpha = update_alpha(alpha, B, mu, inv_covariance, subject_to_B)
#Covariance = update_Covariance(B, alpha, mu, subject_to_B, diag(dim(Z)[2]), 3, 1000)
inv_covariance = solve(Covariance)
Covariance

Covariance1 = diag(2)
subject_to_B1 = sapply(names(subject_to_B), function(l){
  unname(unlist(subject_to_B[l])) + 1
}, simplify = TRUE)


alpha1 = rep(1, length(unique(subject_id)))
B1 = update_B(re, mu, Z, subject_to_B1, Covariance1, sigma, alpha1)
Covariance1 = update_Covariance(B1, alpha1, mu, subject_to_B, diag(dim(Z)[2]), 3, 2)
alpha1 = update_alpha(alpha1, B1, mu_s, Covariance1, subject_to_B1)
Covariance1







      
    
  starts = time

  
  sourceCpp("./src/nDP2.cpp")
  your_code = function(){
    X = as.matrix(c(rnorm(1000, 2), rnorm(4000, -5), rnorm(3000, -2)))
    cluster_id = as.character(c(rep(2, 1000), rep(-5, 4000), rep(-2, 3000)))
    a = nDP(list(p = 1), cluster_id, M_alpha = 100, M_beta = 100, L = 10, K = 35, TRUE)
    
  }
  

  runtime <- system.time(your_code())
  
  # Print the elapsed time
  print(paste("Elapsed time: ", runtime[3], "seconds"))








library(tidyverse)
library(BART3)
Y = BostonHousing$crim
X = dplyr::select(BostonHousing, zn, indus)
tree = gbart(X, Y, verbose = 0, ndpost=10, nskip = 1)

tree = update_BART(as.matrix(X), Y)





#X = rbind(as.matrix(rmvnorm(1000, mean = rep(0, 10))), as.matrix(rmvnorm(1000, mean = rep(-10, 10))))
X = as.matrix(c(rnorm(1000), rnorm(1000, 10)))
subject_id = c(rep(1, 1000), rep(2, 1000))




DP_mcmc = function(X, nburn = 1000, npost = 1000, seed = 1){
  set.seed(seed)
  post_samples1 = array(NA, dim = c(npost, dim(X)[1]))
  
  a = nDP(list(p = 1), subject_id, CDP = FALSE, L = 10, K = 5)
  # a = DP(list(p = 1, mu = 1), N_truncated = 10, N_sample = dim(X)[1], CDP = TRUE)
  a$sigma = 2
  
  for (i in 1:nburn) {
    cat(i, "\n")
    a = update_nDP_normal(X, a, L_alpha = -1, U_alpha = 2, L_beta = -1, U_beta = 2, tol = 1e-20, seed = seed + i)
    #a = update_DP_normal(X, a, L = -1, U = 2, seed = seed + i, tol = 1e-20)
  }
  
  for (i in 1:npost) {
    cat(i, "\n")
    a = update_nDP_normal(X, a, L_alpha = -1, U_alpha = 2, L_beta = -1, U_beta = 2, tol = 1e-20, seed = seed + i + nburn)
    #a = update_DP_normal(X, a, L = -1, U = 2, seed = seed + i + nburn, tol = 1e-20)
    post_samples1[i,] = a$samples
  }
  return(post_samples1)
}


chain1 = DP_mcmc(X, 12000, 12000, 1)
chain2 = DP_mcmc(X, 12000, 12000, 2)
chain3 = DP_mcmc(X, 12000, 12000, 3)






combinedchains = mcmc.list(mcmc(chain1), mcmc(chain2), mcmc(chain3))
a = gelman.diag(combinedchains, multivariate = FALSE)


j = 0
for (i in 1:dim(chain1)[2]) {
  combinedchains = mcmc.list(mcmc(chain1[,i]), mcmc(chain2[,i]), mcmc(chain3[,i]))
  if(gelman.diag(combinedchains)[[1]][1] > 1.1){
    print(paste(i, gelman.diag(combinedchains)[[1]][1]))
    j = j + 1
  }
}


dp = CDP(function(n, p) {
  if (p == 1)
    as.matrix(rnorm(n, 0, 1))
  else{
    rmvnorm(n, rep(0, p))
  }
}, M = 1, n = 1000, n_cluster = 100)







sourceCpp("./src/DP_LMM_BART2.cpp")
X = rmvnorm(100, c(-1,2,-3,4))
subject_id = rep(1:10, 10)
Y = (2* exp(X[,1]) + -1* X[,2] + 3* X[,3]^2 + -1.5*sin(X[,4]) - 25) / 7
Y = sigmoid(Y) >= 0.5

a = DP_LMM_BART(1, 100, 100, Y, as.matrix(X), NULL, subject_id, binary = TRUE, verbose = TRUE, DP_re = TRUE)




# chain: BART 1, pure BART
sourceCpp("./src/BMTrees_MCMC.cpp")
chain1_pure = foreach(1:2, .combine = list) %dopar%{
  BMTrees_MCMC(6000, 6000, Y, as.matrix(X), as.matrix(Z), subject_id, FALSE, verbose =  F, NULL, tol = 1e-100, F, F)
}

chain1_cdp_residual = foreach(1:2, .combine = list) %dopar%{
  BMTrees_MCMC(6000, 6000, Y, as.matrix(X), as.matrix(Z), subject_id, FALSE, verbose =  F, NULL, tol = 1e-100, T, F)
}

chain1_cdp_double = foreach(1:2, .combine = list) %dopar%{
  BMTrees_MCMC(6000, 6000, Y, as.matrix(X), as.matrix(Z), subject_id, FALSE, verbose =  F, NULL, tol = 1e-100, T, T)
}


# chain: BART 5, pure BART
sourceCpp("./src/BMTrees_MCMC.cpp")
chain5_pure = foreach(1:2, .combine = list) %dopar%{
  BMTrees_MCMC(6000, 6000, Y, as.matrix(X), as.matrix(Z), subject_id, FALSE, verbose =  F, NULL, tol = 1e-100, F, F)
}

chain5_cdp_residual = foreach(1:2, .combine = list) %dopar%{
  BMTrees_MCMC(6000, 6000, Y, as.matrix(X), as.matrix(Z), subject_id, FALSE, verbose =  F, NULL, tol = 1e-100, T, F)
}

chain5_cdp_double = foreach(1:2, .combine = list) %dopar%{
  BMTrees_MCMC(6000, 6000, Y, as.matrix(X), as.matrix(Z), subject_id, FALSE, verbose =  F, NULL, tol = 1e-100, T, T)
}


# chain: BART 10, pure BART
sourceCpp("./src/BMTrees_MCMC.cpp")
chain10_pure = foreach(1:2, .combine = list) %dopar%{
  BMTrees_MCMC(6000, 6000, Y, as.matrix(X), as.matrix(Z), subject_id, FALSE, verbose =  F, NULL, tol = 1e-100, F, F)
}

chain10_cdp_residual = foreach(1:2, .combine = list) %dopar%{
  BMTrees_MCMC(6000, 6000, Y, as.matrix(X), as.matrix(Z), subject_id, FALSE, verbose =  F, NULL, tol = 1e-100, T, F)
}

chain10_cdp_double = foreach(1:2, .combine = list) %dopar%{
  BMTrees_MCMC(6000, 6000, Y, as.matrix(X), as.matrix(Z), subject_id, FALSE, verbose =  F, NULL, tol = 1e-100, T, T)
}

