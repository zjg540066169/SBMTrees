#ifndef ARMADILLO_H_
#define ARMADILLO_H_
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#endif

#ifndef UTILS_H_
#define UTILS_H_
#include "utils.h"
#endif


#ifndef BMTREES_H_
#define BMTREES_H_
#include "bmtrees.h"
#endif

#ifndef UPDATE_H_
#define UPDATE_H_
//#include "DP_parallel.h"
#include "nDP.h"
#include "cal_random_effects.h"
#include "update_alpha.h"
#include "DP.h"
#include "update_B.h"
#include "update_Covariance.h"
#include <cmath>
#endif


#ifndef RCPP_H_
#define RCPP_H_
#include <Rcpp.h>
#endif

// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <progress_bar.hpp>


using namespace Rcpp;



/*** R
library(stats)
library(BART3)
library(stats)
library(MCMCpack)
library(mvtnorm)
library(cascsim)
library(truncnorm)
library(hbmem)
library(foreach)
library(doParallel)
library(Matrix)

library(lme4)

makePositiveDefinite <- function(A, epsilon = 1e-8) {
  A_reg <- A
  diag(A_reg) <- diag(A_reg) + epsilon
  return(A_reg)
}



get_inverse_wishart_matrix2 = function(X, Y, Z, subject_id, subject_to_B, binary = F){
  
  constant_cols <- apply(X, 2, function(x) var(x) == 0)
  non_constant_cols <- !constant_cols
  X[, constant_cols] = scale(X[, constant_cols], center = TRUE, scale = FALSE)
  X[, non_constant_cols] <- scale(X[, non_constant_cols])
  
  if(!binary){
    suppressMessages(lmm <- lmer(Y ~ 0 + X + (0 + Z|subject_id), REML = T))
    coe = as.matrix(ranef(lmm)[[1]])
    coe = (coe[names(subject_to_B),])
  }else{
    suppressMessages(lmm <- glmer(as.factor(Y) ~ 0 + X + (0 + Z|subject_id), family = binomial(link = "logit")))
    coe = as.matrix(ranef(lmm)[[1]])
    coe = (coe[names(subject_to_B),])
  }
  co = as.matrix(bdiag(VarCorr(lmm)))
  svd_A <- svd(co)
  tolerance = 1e-10 * max(svd_A$d)
  while(kappa(co) > 500 | !isPositiveDefinite(co)){
    if(tolerance > max(svd_A$d)){
      print("overflow")
      if(!isPositiveDefinite(co)){
        eigenvalues <- eigen(co)$values
        negative_eigenvalues <- eigenvalues[eigenvalues < 0]
        
        # Check if there are any negative eigenvalues
        if (length(negative_eigenvalues) > 0) {
          # Find the negative eigenvalue with the largest absolute value
          largest_abs_negative_eigenvalue <- negative_eigenvalues[which.max(abs(negative_eigenvalues))]
          diag(co) = diag(co) - largest_abs_negative_eigenvalue + 1e-10
          svd_A <- svd(co)
          tolerance = 1e-10 * max(svd_A$d)
          next
        }else{
          break
        }
      }else{
        break
      }
    }
    tolerance = tolerance * 2
    svd_A <- svd(co)
    D_truncated <- diag(svd_A$d)
    diag(D_truncated) = pmax(diag(D_truncated), tolerance)
    co <- svd_A$u %*% D_truncated %*% t(svd_A$v)
    co = makePositiveDefinite(co)
    co = (co + t(co)) / 2
    
  }
  
  # print(A_reconstructed)
  #print(solve(cov(as.matrix(coe))))
  #print(subject_to_B)
  #print(coe)
  #print(chol2inv(chol((as.matrix(bdiag(VarCorr(lmm)))))))
  #print(solve(as.matrix(bdiag(VarCorr(lmm)))))
  #return(list(coe = as.matrix(coe), sigma = sigma(lmm), covariance = as.matrix(bdiag(VarCorr(lmm)))))
  #return(list())
  return(list(coe = as.matrix(coe), sigma = sigma(lmm), covariance = as.matrix(co)))
}

nCDP_residual_update = function(tau_train_data_list, ncores, seed = NULL, step = 0){
  #print(123)
  if(!is.null(seed))
    set.seed(seed + step)
  registerDoParallel(ncores)
  #print(456)
  tau_list = foreach(i = 1:(length(tau_train_data_list)), .multicombine = T, .combine = list) %dopar% {
    residual = tau_train_data_list[[i]]$residual
    tau = tau_train_data_list[[i]]$tau
    tau = update_DP_normal(residual, tau, 0, 0.5, 1e-40)
    return(tau)
  }
 # print(789)
  stopImplicitCluster()
  return(tau_list)
}

CDP_re_update = function(B_tau_train_data_list, ncores, seed = NULL, step = 0){
  
  if(!is.null(seed))
    set.seed(seed + step)
  #print(123)
  registerDoParallel(ncores)
  #print(123)
  B_tau_list = foreach(i = 1:(length(B_tau_train_data_list)), .multicombine = T, .combine = list) %dopar% {
    B = B_tau_train_data_list[[i]]$B
    B_tau = B_tau_train_data_list[[i]]$B_tau
    B_tau = update_DP_normal(B, B_tau, 0, 0.5, 1e-40)
    return(B_tau)
  }
  
  stopImplicitCluster()
  return(B_tau_list)
}

update_B = function(subject_id, subject_to_B, Z, R, Covariance, sigma, B_location){
  B = t(sapply(unique(subject_id), function(j){
    b_sub = subject_to_B[[as.character(j)]] + 1
    Zi = Z[subject_id == as.numeric(j), ]
    R_i = R[subject_id == j]
    var = solve(solve(Covariance) + t(Zi) %*% (Zi) / (sigma))
    mu = (t(t(Zi) %*% R_i) / (sigma) + B_location[b_sub,] %*% solve(Covariance))  %*% var 
    mvrnorm(1, mu, var)
  }))
  # while(!is.positive.definite(t(B - B_location) %*% (B - B_location))){
  #   seed = seed + 1
  #   set.seed(seed)
  #   B = t(sapply(unique(subject_id), function(j){
  #     b_sub = subject_to_B[as.character(j)]
  #     Zi = Z[subject_id == as.numeric(j), ]
  #     R_i = R[subject_id == j]
  #     var = solve(solve(Covariance) + t(Zi) %*% (Zi) / (sigma))
  #     mu = (t(t(Zi) %*% R_i) / (sigma) + B_location[b_sub,] %*% solve(Covariance))  %*% var 
  #     mvrnorm(1, mu, var)
  #   }))
  #   
  # }
  rownames(B) = NULL
  return(B)
}

update_Covariance = function(B, B_location, df, n_subject, inverse_wishart_matrix){
  Covariance = riwish(df + n_subject, inverse_wishart_matrix + t(B - B_location) %*% (B - B_location))
  return(Covariance)
}

predict_ncdp = function(n, subject_id_test, cluster_map, cdp_information, ncores){
  #print(cdp_information)
  registerDoParallel(ncores)
  e = foreach(i = 1:n, .multicombine = T, .combine = c) %dopar% {
    subject = subject_id_test[i]
    ncdp = cluster_map[[subject]] + 1
    cdp_select = cdp_information[[ncdp]]
    values = cdp_select[["y"]]
    pi = cdp_select[["pi"]]
    sample(values, 1, FALSE, pi)
  } 
  stopImplicitCluster()
  return(e)
}


bartModelMatrix=function(X, numcut=0L, usequants=FALSE, type=7,
                         rm.const=FALSE, cont=FALSE, xinfo=NULL) {
  #print(xinfo)
  X.class = class(X)[1]
  
  if(X.class=='factor') {
    X.class='data.frame'
    X=data.frame(X=X)
  }
  grp=NULL
  if(X.class=='data.frame') {
    print(X)
    p=dim(X)[2]
    
    xnm = names(X)
    for(i in 1:p) {
      print(i)
      if(is.factor(X[[i]])) {
        #print(i)
        Xtemp = nnet::class.ind(X[[i]])
        colnames(Xtemp) = paste(xnm[i],1:ncol(Xtemp),sep='')
        X[[i]]=Xtemp
        m=ncol(Xtemp)
        grp=c(grp, rep(m, m))
      } else {
        X[[i]]=cbind(X[[i]])
        colnames(X[[i]])=xnm[i]
        grp=c(grp, 1)
        ##grp=c(grp, i)
      }
    }
    Xtemp=cbind(X[[1]])
    if(p>1) for(i in 2:p) Xtemp=cbind(Xtemp, X[[i]])
    X=Xtemp
  }
  else if(X.class=='numeric' | X.class=='integer') {
    X=cbind(as.numeric(X))
    ##grp=1
  }
  else if(X.class=='NULL') return(X)
  else if(X.class!='matrix')
    stop('Expecting either a factor, a vector, a matrix or a data.frame')
  
  #print(789)
  N <- nrow(X)
  p <- ncol(X)
  
  xinfo. <- matrix(nrow=p, ncol=numcut)
  nc <- numcut
  rm.vars <- c()
  #print(789)
  if(length(xinfo)==0 & N>0 & p>0 & (rm.const | numcut[1]>0)) {
    for(j in 1:p) {
      X.class <- class(X[1, j])[1]
      
      if(X.class=='numeric' | X.class=='integer') {
        xs <- unique(sort(X[ , j]))
        k <- length(xs)
        nc[j] <- numcut
        
        if(k %in% 0:1) { # deal with constant variables
          rm.vars <- c(rm.vars, -j)
          nc[j] <- 1
          if(k==0) xs <- NA
        }
        else if(cont)
          xs <- seq(xs[1], xs[k], length.out=numcut+2)[-c(1, numcut+2)]
        else if(k<numcut) {
          xs <- 0.5*(xs[1:(k-1)]+xs[2:k]) #  if k < numcut, use middle point between values to split
          nc[j] <- k-1
        }
        else if(usequants) { 
          xs <- quantile(X[ , j], type=type,
                         probs=(0:(numcut+1))/(numcut+1))[-c(1, numcut+2)]
          names(xs) <- NULL
        }
        else xs <-
          seq(xs[1], xs[k], length.out=numcut+2)[-c(1, numcut+2)]
      }
      else
        stop(paste0('Variables of type ', X.class, ' are not supported'))
      
      xinfo.[j, 1:nc[j] ] <- xs
    }
  }
  
  X <- data.matrix(X)
  #print(xinfo)
  #print(dim(xinfo.))
  if(length(xinfo)>0) {
    if(is.list(xinfo)) for(j in 1:p){
      
      xinfo.[j, 1:length(xinfo[[j]])] <- xinfo[[j]]
    } 
    else if(is.matrix(xinfo)) xinfo. <- xinfo
    else stop('Only a list or a matrix can be provided for xinfo')
    
    for(j in 1:p) nc[j] <- sum(!is.na(xinfo.[j, ]))
  }
  
  xinfo <- xinfo.
  
  if(rm.const && length(rm.vars)>0 &&
     !(length(rm.vars)==p && all((1:p)==(-rm.vars)))) {
    X <- X[ , rm.vars]
    nc <- nc[rm.vars]
    xinfo <- xinfo[rm.vars, ]
    grp <- grp[rm.vars]
  }
  else if(length(rm.vars)==0 || (length(rm.vars)==p && all((1:p)==(-rm.vars))))
    rm.vars <- 1:p
  
  dimnames(xinfo) <- list(dimnames(X)[[2]], NULL)
  
  if(all(numcut==0)) return(X)
  else return(list(X=X, numcut=as.integer(nc), rm.const=rm.vars,
                   xinfo=xinfo, grp=grp))
}
*/



// [[Rcpp::export]]
List BMTrees_MCMC(LogicalVector obs_ind, long nburn, long npost, NumericVector Y, NumericMatrix X, Nullable<NumericMatrix> Z, CharacterVector subject_id, bool binary = false, bool verbose = true, Nullable<long> seed = R_NilValue, double tol = 1e-40, bool CDP_residual = false, bool CDP_re = false, int pattern = 0, double pi_CDP = 0.99){
  // if(CDP_residual){
  //   bmtrees model = bmtrees(Y, X, Z, subject_id, binary, nCDP_residual, true, CDP_re, tol, seed);
  // }else{
  //   bmtrees model = bmtrees(Y, X, Z, subject_id, binary, nCDP_residual, false, CDP_re,  tol, seed);
  // }
  //Rcout << 123;
  NumericMatrix Z_obs;
  NumericMatrix Z_test;
  NumericVector Y_obs = Y[obs_ind];
  NumericMatrix X_obs = row_matrix(X, obs_ind);
  if(!Z.isNull()){
    NumericMatrix z = as<NumericMatrix>(Z);
    //Rcout << z.rows() << " " << z.cols()<<std::endl;
    //Rcout << obs_ind.length()<<std::endl;
    //Rcout << z << std::endl;
    Z_obs = row_matrix(z, obs_ind);
    Z_test = row_matrix(z, !obs_ind);
  }
  //Rcout << Z_obs << std::endl;
    
  CharacterVector subject_id_obs = subject_id[obs_ind];
  IntegerVector row_id_obs = seqC(1, Y.length())[obs_ind];
  //return List::create();
  bmtrees model = bmtrees(clone(Y_obs), clone(X_obs), clone(Z_obs), subject_id_obs, row_id_obs, binary, CDP_residual, CDP_re, tol, 50,pattern, pi_CDP);
  //return List::create();
  //
  //
  NumericVector Y_test = Y[!obs_ind];
  NumericMatrix X_test = row_matrix(X, !obs_ind);

  CharacterVector subject_id_test = subject_id[!obs_ind];
  IntegerVector row_id_test = seqC(1, Y.length())[!obs_ind];

  int d;
  if(Z.isNull()){
    d = 1;
  }else{
    d = as<NumericMatrix>(Z).ncol();
  }

  long N = Y_obs.length();
  long N_test = Y_test.length();
  int p = X_obs.ncol();
  int n_subject = unique(subject_id_obs).length();
  
  
  List post_trees;
  NumericMatrix post_tree_pre_mean(npost, 1);
  NumericMatrix post_M(npost, 1);
  NumericMatrix post_M_re(npost, 1);
  NumericMatrix post_alpha(npost, n_subject);
  NumericMatrix post_x_hat(npost, N);
  NumericMatrix post_sigma(npost, 1);
  NumericMatrix post_lambda(npost, 1);
  NumericMatrix post_Sigma(npost, d * d);
  NumericMatrix post_B(npost, n_subject * d);
  NumericMatrix post_tau_samples(npost, N);
  NumericMatrix post_B_tau_samples(npost, n_subject * d);
  NumericMatrix post_random_effect(npost, N);
  NumericMatrix post_y_predict(npost, N);
  NumericMatrix post_y_predict_new(npost, N);
  NumericMatrix post_y_expectation(npost, N);
  NumericMatrix post_y_sample(npost, N);
  NumericMatrix post_tau_position(npost, (int)sqrt(N));
  NumericMatrix post_tau_pi(npost, (int)sqrt(N));
  NumericMatrix post_B_tau_pi(npost, (int)sqrt(n_subject));
  NumericMatrix post_B_tau_position(npost, (int)sqrt(n_subject) * d);

  NumericMatrix post_y_expectation_test(npost, N_test);
  NumericMatrix post_y_sample_test(npost, N_test);

 //Progress probar(nburn + npost, true);
  List tau;
  List B_tau;
  Progress progr(nburn + npost, !verbose);
  for(int i = 0 ; i < nburn + npost; ++i){

    if (Progress::check_abort() )
      return -1.0;
    progr.increment();

    model.update_all(verbose);

    if(verbose)
    Rcout << i << " " << nburn + npost << std::endl;




    // dmodel.update();
    List post_sample = model.posterior_sampling();
    if(i >= nburn){
      //return model.predict(clone(X), clone(Z), subject_id);
      // post_trees = post_sample["tree"];
      //
      // post_M(i - nburn, 0) = post_sample["M"];
      // post_M_re(i - nburn, 0) = post_sample["M_re"];
      post_tree_pre_mean(i - nburn, 0) = post_sample["tree_pre_mean"];
      post_sigma(i - nburn, 0) = post_sample["sigma"];
      //
      // post_alpha(i - nburn, _) = as<NumericVector>(post_sample["alpha"]);
      post_x_hat(i - nburn, _) = as<NumericVector>(post_sample["tree_pre"]);
      post_Sigma(i - nburn, _) = as<NumericVector>(post_sample["Sigma"]);
      // //Rcout << post_Sigma(i - nburn, _)  << std::endl;
      post_B(i - nburn, _) = as<NumericVector>(post_sample["B"]);
      // post_tau_samples(i - nburn, _) = as<NumericVector>(post_sample["tau_samples"]);
      // post_B_tau_samples(i - nburn, _) = as<NumericVector>(post_sample["B_tau_samples"]);
      post_random_effect(i - nburn, _) = as<NumericVector>(post_sample["re"]);
      // //post_y_predict(i - nburn, _) = as<NumericVector>(post_sample["y_predict"]);
      // //post_y_predict_new(i - nburn, _) = model.predict(clone(X), clone(Z), subject_id, seq(1, Y.length()));
      post_y_expectation(i - nburn, _) = model.predict_expectation(clone(X_obs), clone(Z_obs), subject_id_obs, row_id_obs);
      post_y_sample(i - nburn, _) = model.predict_sample(clone(X_obs), clone(Z_obs), subject_id_obs, row_id_obs);
     

      post_y_expectation_test(i - nburn, _) = model.predict_expectation(clone(X_test), clone(Z_test), subject_id_test, row_id_test);
      post_y_sample_test(i - nburn, _) = model.predict_sample(clone(X_test), clone(Z_test), subject_id_test, row_id_test);

      if(CDP_residual){
        tau = post_sample["tau"];
        //B_tau = post_sample["B_tau"];
        post_tau_position(i - nburn, _) = as<NumericVector>(tau["y"]);
        post_tau_pi(i - nburn, _) = as<NumericVector>(tau["pi"]);
      }
      if (CDP_re){
        B_tau = post_sample["B_tau"];
        post_B_tau_pi(i - nburn, _) = as<NumericVector>(B_tau["pi"]);
        post_B_tau_position(i - nburn, _) = as<NumericVector>(B_tau["y"]);
        post_lambda(i - nburn, 0) = (double)B_tau["lambda"];
      }

    }
 //probar.increment();
  }
 return List::create(
   Named("post_tree_pre_mean") = post_tree_pre_mean,
   // Named("post_trees") = post_trees,
   // Named("post_tau") = tau,
   // Named("post_B_tau") = B_tau,
   // Named("post_M") = post_M,
   // Named("post_M_re") = post_M_re,
   Named("post_sigma") = post_sigma,
   // Named("post_alpha") = post_alpha,
   Named("post_x_hat") = post_x_hat,
   Named("post_Sigma") = post_Sigma,
   Named("post_lambda") = post_lambda,
   Named("post_B") = post_B,
   // Named("post_tau_samples") = post_tau_samples,
   // Named("post_B_tau_samples") = post_B_tau_samples,
   Named("post_random_effect") = post_random_effect,
   // Named("post_y_predict") = post_y_predict,
   // Named("post_y_predict_new") = post_y_predict_new,
   Named("post_y_expectation") = post_y_expectation,
   Named("post_y_sample") = post_y_sample,
   Named("post_tau_position") = post_tau_position,
   Named("post_tau_pi") = post_tau_pi,
   Named("post_B_tau_position") = post_B_tau_position,
   Named("post_B_tau_pi") = post_B_tau_pi,
   Named("post_y_expectation_test") = post_y_expectation_test,
   Named("post_y_sample_test") = post_y_sample_test
 );

}



// 
// 
// 
// // [[Rcpp::export]]
// List BMTrees_MCMC2(long nburn, long npost, NumericVector Y, NumericMatrix X, NumericMatrix Z, CharacterVector subject_id,  LogicalVector train, bool binary = false, bool verbose = true, Nullable<long> seed = R_NilValue, double tol = 1e-40, bool nCDP_residual = false, bool CDP_re = false){
//   // if(CDP_residual){
//   //   bmtrees model = bmtrees(Y, X, Z, subject_id, binary, nCDP_residual, true, CDP_re, tol, seed);
//   // }else{
//   //   bmtrees model = bmtrees(Y, X, Z, subject_id, binary, nCDP_residual, false, CDP_re,  tol, seed);
//   // }
//   
//   NumericMatrix X_train = row_matrix(X, train);
//   NumericMatrix Z_train = row_matrix(Z, train);
//   NumericVector Y_train = Y[train];
//   CharacterVector subject_train = subject_id[train];
//   Rcpp::Environment G = Rcpp::Environment::global_env();
//   Rcpp::Function nCDP_residual_update = G["nCDP_residual_update"];
//   Rcpp::Function CDP_re_update = G["CDP_re_update"];
//   bmtrees model = bmtrees(clone(Y_train), clone(X_train), clone(Z_train), subject_train, binary, nCDP_residual, CDP_re, tol);
//   //Rcout << 123;
//   //bmtrees(clone(y_train), clone(X_train), clone(Z), subject_id, type[i+1], nCDP_residual, CDP_re, tol, seed)
//   Environment base("package:base");
//   
//   if(seed.isNotNull()){
//     Rcpp::Function set_seed_r = base["set.seed"];
//     set_seed_r(seed);
//     //long seed_ = (long)seed;
//   }
//   
//   int d = as<NumericMatrix>(Z).ncol() + 1;
//   // if(Z.isNull()){
//   //   d = 1;
//   // }else{
//   //   d = as<NumericMatrix>(Z).ncol() + 1;
//   // }
//   
//   long N = Y.length();
//   int p = X.ncol();
//   int n_subject = unique(subject_id).length();
//   
//   
//   List post_trees;
//   NumericMatrix post_M_alpha(npost, 1);
//   NumericMatrix post_M_re(npost, 1);
//   NumericMatrix post_alpha(npost, n_subject);
//   NumericMatrix post_x_hat(npost, N);
//   NumericMatrix post_sigma(npost, 1);
//   NumericMatrix post_Sigma(npost, d * d);
//   NumericMatrix post_B(npost, n_subject * d);
//   NumericMatrix post_tau_samples(npost, N);
//   NumericMatrix post_B_tau_samples(npost, n_subject * d);
//   NumericMatrix post_random_effect(npost, N);
//   NumericMatrix post_y_predict(npost, N);
//   NumericMatrix post_y_predict_function(npost, N);
//   
//   //Progress probar(nburn + npost, true);
//   List tau;
//   List B_tau;
//   for(int i = 0 ; i < nburn + npost; ++i){
//     Rcout << i << " " << nburn + npost << std::endl;
//     model.update_tree();
//     if(nCDP_residual){
//       //Rcout << "nCDP_residual" << std::endl;
//       List tau_train_data_list;
//       tau_train_data_list.push_back(model.get_nCDP_residual_data());
//       List tau_list = nCDP_residual_update(tau_train_data_list, 1, seed, i);
//       //Rcout << tau_list.length();
//       model.set_nCDP_residual_data(tau_list);
//     }
//     
//     if(CDP_re){
//       //Rcout << "CDP_re" << std::endl;
//       List B_tau_train_data_list;
//       B_tau_train_data_list.push_back(model.get_CDP_re_data());
//       List B_tau_list = CDP_re_update(B_tau_train_data_list, 1, seed, i);
//       model.set_CDP_re_data(B_tau_list);
//     }
//     model.update(verbose);
//     List post_sample = model.posterior_sampling();
//     if(i >= nburn){
//       post_trees = post_sample["tree"];
//       post_M_alpha(i - nburn, 0) = post_sample["M_alpha"];
//       post_M_re(i - nburn, 0) = post_sample["M_re"];
//       post_sigma(i - nburn, 0) = post_sample["sigma"];
//       
//       post_alpha(i - nburn, _) = as<NumericVector>(post_sample["alpha"]);
//       post_x_hat(i - nburn, _) = as<NumericVector>(post_sample["tree_pre"]);
//       post_Sigma(i - nburn, _) = as<NumericVector>(post_sample["Sigma"]);
//       post_B(i - nburn, _) = as<NumericVector>(post_sample["B"]);
//       post_tau_samples(i - nburn, _) = as<NumericVector>(post_sample["tau_samples"]);
//       post_B_tau_samples(i - nburn, _) = as<NumericVector>(post_sample["B_tau_samples"]);
//       post_random_effect(i - nburn, _) = as<NumericVector>(post_sample["re"]);
//       post_y_predict(i - nburn, _) = as<NumericVector>(post_sample["y_predict"]);
//       post_y_predict_function(i - nburn, _) = as<NumericVector>(model.predict(clone(X), clone(Z), subject_id));
//       tau = post_sample["tau"];
//       B_tau = post_sample["B_tau"];
//     }
//     //probar.increment();
//   }
//   return List::create(
//     Named("post_trees") = post_trees,
//     Named("post_tau") = tau,
//     Named("post_B_tau") = B_tau,
//     Named("post_M_alpha") = post_M_alpha,
//     Named("post_M_re") = post_M_re,
//     Named("post_sigma") = post_sigma,
//     Named("post_alpha") = post_alpha,
//     Named("post_x_hat") = post_x_hat,
//     Named("post_Sigma") = post_Sigma,
//     Named("post_B") = post_B,
//     Named("post_tau_samples") = post_tau_samples,
//     Named("post_B_tau_samples") = post_B_tau_samples,
//     Named("post_random_effect") = post_random_effect,
//     Named("post_y_predict") = post_y_predict,
//     Named("post_y_predict_function") = post_y_predict_function
//   );
//   
// }
// 
