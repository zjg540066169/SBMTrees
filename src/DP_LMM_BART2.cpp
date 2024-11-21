#ifndef UTILS_H_
#define UTILS_H_
#include "utils.h"
#endif



#include <Rcpp.h>
#include "create_subject_to_B.h"
#include "DP.h"
#include "nDP.h"
#include "cal_random_effects.h"
#include "update_alpha.h"
//#include "update_BART.h"
#include "update_B.h"
#include "update_Covariance.h"
#include <cmath>

using namespace Rcpp;


/*** R
library(stats)
library(BART3)
library(stats)
library(MCMCpack)
library(mvtnorm)
library(cascsim)
library(truncnorm)



update_tree = function(X, y){
  #print(y)
  tree = gbart(X, y, verbose = 0, ndpost=1, nskip = 1)
  return(tree)
  #tree_pre = tree$yhat.train.mean#bart_fit[[1]]
}
*/


// [[Rcpp::export]]
List DP_LMM_BART(long chain, long nburn, long npost, NumericVector Y, NumericMatrix X, Nullable<NumericMatrix> Z, CharacterVector subject_id, bool binary = false, bool verbose = true, Nullable<long> seed = R_NilValue, double tol = 1e-40, bool nDP_residual = false, bool DP_re = false, Nullable<List> last_states = R_NilValue){
  
  
  Environment base("package:base");
  Function solve = base["solve"];
  Rcpp::Environment G = Rcpp::Environment::global_env();
  Rcpp::Environment stats("package:stats");
  Rcpp::Function dnorm_cpp = stats["dnorm"];
  Rcpp::Function pnorm_cpp = stats["pnorm"];
  Environment MCMCpack("package:MCMCpack");
  Function rinvgamma = MCMCpack["rinvgamma"];
  Function riwish = MCMCpack["riwish"];
  Environment mvtnorm("package:mvtnorm");
  Function rmvnorm = mvtnorm["rmvnorm"];
  Rcpp::Environment truncnorm("package:truncnorm");
  Rcpp::Function rtruncnorm = truncnorm["rtruncnorm"];
  Function update_tree = G["update_tree"];
  
  if(seed.isNotNull()){
    Rcpp::Function set_seed_r = base["set.seed"];
    set_seed_r(seed);
    //long seed_ = (long)seed;
  }
  
  int d;
  if(Z.isNull()){
    d = 1;
  }else{
    d = as<NumericMatrix>(Z).ncol() + 1;
  }
  
  long N = Y.length();
  int p = X.ncol();
  int n_subject = unique(subject_id).length();
  NumericMatrix z(N, d);
  z(_, 0) = z(_, 0) + 1;
  if(!Z.isNull()){
    NumericMatrix z0 = as<NumericMatrix>(Z);
    for(int i = 1; i < d; ++i){
      z(_, i) = z0(_, i - 1);
    }
  }
  
  
  NumericVector indicator;
  if(binary){
    indicator = clone(Y);
    NumericVector y_ = clone(Y);
    y_ = y_ * 2 - 1;
    Y = clone(y_);
  }
  
  
  
  List post_trees;
  NumericMatrix post_M_alpha(npost, 1);
  NumericMatrix post_M_beta(npost, 1);
  NumericMatrix post_M_re(npost, 1);
  NumericMatrix post_alpha(npost, n_subject);
  NumericMatrix post_x_hat(npost, N);
  NumericMatrix post_sigma(npost, 1);
  NumericMatrix post_Sigma(npost, d * d);
  NumericMatrix post_B(npost, n_subject * d);
  NumericMatrix post_tau_samples(npost, N);
  NumericMatrix post_B_tau_samples(npost, n_subject * d);
  NumericMatrix post_random_effect(npost, N);
  NumericMatrix post_y_predict(npost, N);
  
  /*
   
   
   
   NumericMatrix post_y_predict(npost, n);
   NumericMatrix post_x_hat(npost, n);
   NumericMatrix post_mixed_effect(npost, n);
   NumericMatrix post_alpha(npost, n_subject);
   List post_tree_para(npost);
   */
  
  
  
  List subject_to_B;
  NumericMatrix inverse_wishart_matrix(d, d);
  NumericMatrix Covariance;
  NumericVector alpha(n_subject);
  double M_re = 0;
  double M_alpha = 0;
  double M_beta = 0;
  
  List tau;
  List B_tau; 
  
  NumericVector tau_samples(N);
  NumericMatrix B(n_subject, d);
  NumericMatrix B_tau_samples(n_subject, d);
  
  
  
  
  if(last_states.isNotNull()){
    if(verbose)
      Rcout << "read from last states" << std::endl;
    List last_state (last_states);
    subject_to_B = last_state["subject_to_B"];
    inverse_wishart_matrix = NumericMatrix(d, d, as<NumericMatrix>(last_state["inverse_wishart_matrix"]).begin()); 
    int last_npost = as<NumericMatrix>(last_state["post_Sigma"]).nrow() - 1;
    Covariance = NumericMatrix(d, d, as<NumericMatrix>(last_state["post_Sigma"])(last_npost, _).begin());
    M_re = as<NumericMatrix>(last_state["post_M_re"])(last_npost, 0);
    if(nDP_residual){
      M_alpha = as<NumericMatrix>(last_state["post_M_alpha"])(last_npost, 0);
      M_beta = as<NumericMatrix>(last_state["post_M_beta"])(last_npost, 0);
      tau = last_state["post_tau"];
      tau_samples = tau["samples"];
    }
    
    
    if(DP_re){
      B_tau = last_state["post_B_tau"];
      B_tau_samples = as<NumericMatrix>(B_tau["samples"]);
    }
    
    alpha = as<NumericMatrix>(last_state["post_alpha"])(last_npost, _);
    
    
    
    B = NumericMatrix(n_subject, d, as<NumericMatrix>(last_state["post_B"])(last_npost, _).begin());
  }else{
    
    subject_to_B = create_subject_to_B(subject_id);
    
    
    inverse_wishart_matrix.fill_diag(1.0);
    
    
    Covariance = riwish(d, inverse_wishart_matrix);
    
    //alpha = 
    //  rinvgamma(n_subject, 1, 1);
    alpha = alpha + 1;
    if(DP_re){
      M_re = pow(n_subject, (double)(runif(1, -1, 2)[0]));
      B_tau = DP(List::create(Named("p") = d), M_re, sqrt(n_subject), n_subject, true);
      B_tau_samples = as<NumericMatrix>(B_tau["samples"]);
    }
    if(nDP_residual){
      M_alpha = pow(sqrt(n_subject), (double)(runif(1, -1, 2)[0]));
      M_beta = pow(sqrt(10), (double)(runif(1, -1, 2)[0]));
      tau = nDP(List::create(Named("p") = 1), subject_id, M_alpha, M_beta, 10, 15, true);
      tau_samples = tau["samples"];
    }  
    
    for(int i = 0 ; i < n_subject; ++i){
      B(i,_) = as<NumericVector>(rmvnorm(1, B_tau_samples(i,_), alpha[i] * Covariance));
    }
  }
  
  
  
  
  //Rcout << B << std::endl;
  //Rcout << chain << " DP2" << std::endl;
  NumericVector re = cal_random_effects(z, subject_id, B, subject_to_B);
  
  
  //Rcout << chain << " DP" << std::endl;
  
  
  
  //return List::create(Named("re") = re);
  //Rcout << 123;
  
  for(int i = 0 ; i < nburn + npost; ++i){
    if(verbose)
      Rcout <<  "Chain:" << chain << " "  << i + 1 << " " << nburn + npost<< std::endl;
    if(verbose)
      Rcout << "update trees" << std::endl;
    //List trees_bart = update_BART(X, );
    
    //Rcout << 123 << std::endl;
    List tree = update_tree(X, Y - re - tau_samples);
    //Rcout << 123 << std::endl;
    
    //Rcout << 456 << std::endl;
    //NumericMatrix tree_pre(X.nrow(), 1, as<NumericVector>(tree["yhat.train.mean"]).begin());
    //return List::create(Named("sigma") = sigma, Named("tree_pre") = tree_pre, Named("trees") = tree);
    
    
    
    
    NumericVector tree_pre = tree["yhat.train.mean"];
    double sigma = 1;
    if(!binary){
      double sigma = tree["sigma.mean"];
    }
    
    if(verbose)
      Rcout << "update residual" << std::endl;
    NumericVector residual_tem = Y - re - tree_pre;
    NumericMatrix residual(N, 1, residual_tem.begin());
    
    
    
    if(nDP_residual){
      if(verbose)
        Rcout << "update nDP" << std::endl;
      tau["sigma"] = sigma;
      tau = update_nDP_normal(residual, tau, -1, 2, -1, 2, tol);
      NumericVector tau_samples = tau["samples"];
      M_alpha = tau["M_alpha"];
      M_beta = tau["M_beta"];
    }
    
    //Rcout << 3 << std::endl;
    if(verbose)
      Rcout << "update Covariance" << std::endl;
    Covariance = update_Covariance(B, alpha, B_tau_samples, subject_to_B, inverse_wishart_matrix, d);
    //Rcout << Covariance << std::endl;
    
    if(DP_re){
      if(verbose)
        Rcout << "update DP" << std::endl;
      B_tau["Sigma"] = Covariance;
      if(d == 1){
        B_tau["sigma"] = Covariance;
      }
      B_tau = update_DP_normal(B, B_tau, -1, 2, tol);
      NumericMatrix B_tau_samples = B_tau["samples"];
      M_re = B_tau["M"];
    }
    //return B_tau_samples;
    //Rcout << Y - as<NumericVector>(tau_samples) - tree_pre << std::endl;
    if(verbose)
      Rcout << "update B" << std::endl;
    B = update_B(Y - tau_samples - tree_pre, z, subject_id, B_tau_samples, subject_to_B, Covariance, sigma, alpha);
    
    if(verbose)
      Rcout << "update alpha" << std::endl;
    //alpha = update_alpha(alpha, B, B_tau_samples, Covariance, subject_to_B);
    //Rcout << alpha << std::endl;
    //update_B(Y - tree_pre - tau_samples, z, subject_id, B_tau_samples, subject_to_B, NumericMatrix Covariance, double sigma, NumericVector alpha)
    if(verbose)
      Rcout << "update random effects" << std::endl;
    re = cal_random_effects(z, subject_id, B, subject_to_B);
    
    if(verbose)
      Rcout << "M_re:" << M_re << "  " << "M_alpha:" << M_alpha << "  " << "M_beta:" << M_beta << std::endl;
    if(binary){
      //Rcout << "sigma:"<<sqrt(sigma_2) <<std::endl;
      for(int i = 0; i < N; ++i){
        if(indicator[i] == 0){
          NumericVector mean_y = rtruncnorm(1, R_NegInf, 0, tree_pre[i] + re[i] + tau_samples[i], sigma);
          Y[i] = mean_y[0];
        }else{
          NumericVector mean_y = rtruncnorm(1, 0, R_PosInf, tree_pre[i] + re [i] + tau_samples[i], sigma);
          Y[i] = mean_y[0];
        }
      }
      if(verbose)
        Rcout << "update probit outcome" << std::endl;
    }
    
    
    if(i >= nburn){
      post_trees = tree;
      post_M_alpha(i - nburn, 0) = M_alpha;
      post_M_beta(i - nburn, 0) = M_beta;
      post_M_re(i - nburn, 0) = M_re;
      post_sigma(i - nburn, 0) = sigma;
      
      post_alpha(i - nburn, _) = alpha;
      post_x_hat(i - nburn, _) = tree_pre;
      post_Sigma(i - nburn, _) = Covariance;
      post_B(i - nburn, _) = B;
      post_tau_samples(i - nburn, _) = tau_samples;
      post_B_tau_samples(i - nburn, _) = B_tau_samples;
      post_random_effect(i - nburn, _) = re;
      post_y_predict(i - nburn, _) = tree_pre + re + tau_samples;
      //NumericVector B__ = B(0, _);
      //Rcout << B__ << std::endl;
    }
    
  }
  
  
  
  
  
  
  
  
  
  //Rcout << "complete LMM" << std::endl;
  
  //(nburn, npost, y, X, z, subject_id, v = 3, k = 1, sigma_tau = 1, tol = 1e-20, verbose = TRUE, seed = NULL, DP = FALSE, N_truncated = 40, REN_truncated = 10, nCDP = TRUE, CDP_randeff = TRUE, normal_invgamma_residual = FALSE
  
  return List::create(
    Named("post_trees") = post_trees,
    Named("post_tau") = tau,
    Named("post_B_tau") = B_tau,
    Named("post_M_alpha") = post_M_alpha,
    Named("post_M_beta") = post_M_beta,
    Named("post_M_re") = post_M_re,
    Named("post_sigma") = post_sigma,
    Named("post_alpha") = post_alpha,
    Named("post_x_hat") = post_x_hat,
    Named("post_Sigma") = post_Sigma,
    Named("post_B") = post_B,
    Named("post_tau_samples") = post_tau_samples,
    Named("post_B_tau_samples") = post_B_tau_samples,
    Named("post_random_effect") = post_random_effect,
    Named("post_y_predict") = post_y_predict,
    Named("subject_to_B") = subject_to_B,
    Named("inverse_wishart_matrix") = inverse_wishart_matrix
  );
  
}