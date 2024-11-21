#ifndef ARMADILLO_H_
#define ARMADILLO_H_
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#endif

#ifndef ARMADILLO_EX_H_
#define ARMADILLO_EX_H_
#include <RcppArmadilloExtensions/sample.h>
// [[Rcpp::depends(RcppArmadillo)]]
#endif


#ifndef DIST_H_
#define DIST_H_
#include <RcppDist.h>
// [[Rcpp::depends(RcppArmadillo, RcppDist)]]
// [[Rcpp::plugins(openmp)]]
#endif

#include <Rcpp.h>
#ifndef UTILS_H_
#define UTILS_H_
#include "utils.h"
#endif
#include <cmath>
#ifndef DP_SAMPLER_H_
#define DP_SAMPLER_H_
#include "DP_sampler_arma.h"
#endif

#ifndef Parallel_H_
#define Parallel_H_
#include <RcppParallel.h>
// [[Rcpp::depends(RcppParallel)]]
#endif


#ifdef _OPENMP
#include <omp.h>
#endif


using namespace Rcpp;
using namespace RcppParallel;
using namespace arma;

// [[Rcpp::export]]
SEXP DP_arma(List parameters, double M, long N_truncated, long N_sample, bool CDP = true){
  
  vec p(N_truncated, fill::zeros);
  IntegerVector index = seq(0, N_truncated - 1);
  IntegerVector cluster;
  
  if(N_truncated > 1){
    vec b(N_truncated, fill::zeros);
    for(int k = 0 ; k < N_truncated; ++k)
      b(k) = R::rbeta(1, M);
    p(0) = b(0);
    for(int i = 1; i <= N_truncated-1; i++){
      p(i) = b(i) * arma::prod(1 - b.subvec(0, i - 1));
    }
  }else{
    p(0) = 1;
  }
  
  // NumericVector p_ = wrap(p);
  // p_.attr("names") = index;
  // 
  
  List samples_parameters = DP_sampler_arma(N_truncated, parameters);
  mat y = samples_parameters["y"];
  
  int dp = parameters["p"];
  parameters["CDP"] = CDP;
  //Rcout << dp << std::endl;
  
  
  if(CDP == true){
    for(size_t i = 0; i < dp; ++i){
      double inner_product = arma::dot(y.col(i), p);
      y.col(i) -= inner_product;  // Subtract the inner product from each element in the column
    }
  }
  //CharacterVector rowname(index.begin(), index.end());
  //return List::create(Named("y") = y);
  cluster = Rcpp::RcppArmadillo::sample(index, N_sample, true, p);
  //return List::create(Named("y") = wrap(y), Named("cluster") = cluster);
  mat samples(N_sample, dp, fill::zeros);
  
  for(int i = 0; i < N_sample; ++i){
    samples.row(i) = y.row(cluster[i]);
  }
  //Rcout << samples << std::endl;
  //return List::create(Named("samples") = samples);
    //NumericMatrix samples_ = wrap(samples);
    //rownames(samples_) = cluster;
    // NumericMatrix y_ = wrap(y);
    // rownames(y_) = rowname;
    // 
  // if(parameters.containsElementNamed("sd"))
  //   return List::create(Named("samples") = (samples), Named("cluster") = cluster, Named("pi") = (p), Named("parameters") = parameters, Named("y") = (y),Named("M") = M, Named("lambda") = parameters["lambda"], Named("sigma") = parameters["sigma"]);
  // else
  //   return List::create(Named("samples") = (samples), Named("cluster") = cluster, Named("pi") = (p), Named("parameters") = parameters, Named("y") = (y),Named("M") = M, Named("lambda") = parameters["lambda"], Named("Sigma") = parameters["Sigma"]);
  Rcpp::Shield<Rcpp::List> result(Rcpp::List::create(Named("samples") = (samples), Named("cluster") = cluster, Named("pi") = (p), Named("parameters") = parameters, Named("y") = (y)));
  return result;
  //return List::create(Named("samples") = wrap(samples), Named("cluster") = cluster, Named("index") = index, Named("pi") = wrap(p), Named("parameters") = parameters, Named("y") = wrap(y), Named("lambda") = parameters["lambda"]);
}

// [[Rcpp::export]]
List update_DP_normal(NumericMatrix X, List tau, double L = -1, double U = 2){
  List parameters = tau["parameters"];
  double M = tau["M"];
  long N = X.nrow();
  if(strcmp(parameters["distribution"], "normal") == 0){
    int p = as<int>(parameters["p"]);
    //Rcout << p << std::endl;
    double sigma;
    NumericMatrix Sigma;
    //if(p == 1){
    if(parameters.containsElementNamed("sd")){
      sigma = tau["sigma"];
    }else{
      Sigma = as<NumericMatrix>(tau["Sigma"]);
    }
    
    
    IntegerVector cluster = tau["cluster"];
    NumericVector pi = tau["pi"];
    NumericMatrix y = tau["y"];
    bool CDP = parameters["CDP"];
    double lambda = tau["lambda"];
    //if(p == 1){
    if(parameters.containsElementNamed("sd")){
      //double b_1 = innerProduct(y - as<double>(parameters["mu"]), y - as<double>(parameters["mu"])) / 2 / pow(sigma, 2);
      //Rcout << b_1 / y.length() << std::endl;
      //double a_1 = y.length() / 2 + as<double>(parameters["a"]);
      //Rcout << innerProduct(y - as<double>(parameters["mu"]), y - as<double>(parameters["mu"])) / 2 / pow(sigma, 2) << std::endl;
      lambda = rinvgamma(y.length() / 2 + as<double>(parameters["a"]), as<double>(parameters["b"]) + innerProduct(y - as<double>(parameters["mu"]), y - as<double>(parameters["mu"])) / 2 / pow(sigma, 2));
    }else{
      //double b_2 = quadratic_form(y, as<NumericVector>(parameters["mu"]), Sigma) / 2;
      //Rcout << b_2 / (p * pi.length()) << std::endl;
      //double a_2 = p * pi.length() / 2 + as<double>(parameters["a"]);
      //Rcout << quadratic_form(y, as<NumericVector>(parameters["mu"]), Sigma) / 2 << std::endl;
      lambda = rinvgamma(y.length() / 2 + as<double>(parameters["a"]), as<double>(parameters["b"]) + quadratic_form(y, as<NumericVector>(parameters["mu"]), Sigma) / 2);
    }
    
   //return List::create();
   for(int i = 0 ; i < X.rows(); ++i){
     NumericVector log_density(y.nrow());
     NumericVector X_row = X(i, _);
     for(int j = 0 ; j < y.nrow(); ++j){
       NumericVector log_dens;
       NumericVector mu = y(j,_);
       //Rcout << mu[0] << std::endl;
       double pj = pi[j];
       if(parameters.containsElementNamed("sd")){
       //if(p == 1){
         log_dens = dnorm(X_row, mu[0], sigma, true);
         log_density[j] = log_dens[0] + log(pj);
       }else{
         NumericMatrix X_ = NumericMatrix(1, X_row.length(), X_row.begin());
         log_dens = wrap(dmvnorm(as<mat>(X_), mu, as<mat>(Sigma), true));
         //Rcout << 456 << std::endl;
         
         log_density[j] = log_dens[0] + log(pj);
         //Rcout << log_density[j] << std::endl;
       }
     }
     //density = set_tol(density, tol);
     log_density = log_density - max(log_density);
     NumericVector density = exp(log_density) / sum(exp(log_density));
     int selected = sample(seqC(0, y.nrow() - 1), 1, false, density)[0];
     cluster[i] = selected;
   }
   
    tau["cluster"] = cluster;
    //Rcout << "check cluster" << std::endl;
    
    IntegerVector nk = table(cluster);
    for(int i = 0; i < y.nrow(); ++i){
      const char * n = std::to_string(i).c_str();
      if(!nk.containsElementNamed(n)){
        nk[n] = 0;
      }
    }
    
    
    CharacterVector p_names = rownames(y);
    nk = nk[p_names];
    NumericMatrix beta(nk.length() - 1, 2);
    beta(_,0) = as<NumericVector>(nk[seqC(0, (nk.length() - 2))]);
    IntegerVector rev_nk = rev(nk);
    
    IntegerVector rev_cumsum = cumsum(rev_nk);
    IntegerVector non_rev_cumsum = rev(rev_cumsum);
    beta(_,1) = as<NumericVector>(non_rev_cumsum[seqC(1, (nk.length() - 1))]);
    NumericVector Vh(nk.length() - 1);
    for(int i = 0; i < nk.length() - 1; ++i){
      double new_beta_sample = rbeta(1, 1 + beta(i, 0), M + beta(i, 1))[0];
      while(new_beta_sample == 1.0){
        new_beta_sample = rbeta(1, 1 + beta(i, 0), M + beta(i, 1))[0];
      }
      Vh[i] = new_beta_sample;
    }
    Vh.push_back(1);
    //return List::create(Named("Vh") = Vh, Named("nk") = nk, Named("pi") = pi);
    pi = NumericVector(nk.length());
    pi[0] = Vh[0];
    NumericVector log_negative_1_cumsum = cumsum(log(1 - Vh));
    //Rcout << negative_1_cumprod << std::endl;
    for(int i = 1; i < nk.length() - 1; ++i){
      pi[i] = Vh[i] * exp(log_negative_1_cumsum[i - 1]);
    }
    
    if (nk.length() > 1){
      pi[nk.length() - 1] = exp(log_negative_1_cumsum[nk.length() - 2]);
      M = as<NumericVector>(rtgamma(1, nk.length() - 1, -1 / log_negative_1_cumsum[nk.length() - 2], pow(N, L), pow(N, U)))[0];
    }else{
      pi[0] = 1;
      M = 0;
    }
    pi.names() = p_names;
    //Rcout << 123 << std::endl;
    //return List::create(Named("Vh") = Vh, Named("nk") = nk, Named("pi") = pi);
    //Rcout << 123 << std::endl;
    //Rcout << p << std::endl;
    if(parameters.containsElementNamed("sd")){
      //NumericMatrix a(y.nrow(), 2);
      //double inv_sigma_2 = std::pow(1 / sigma, 2);
      for(int i = 0 ; i < y.nrow(); ++i){
        if(nk[i] == 0){
          double var = std::pow(sigma, 2) * lambda;
          double mean = parameters["mu"];
          
          double sample = rnorm(1, mean, sqrt(var))[0];
          //return List::create(Named("mean") = mean, Named("var") = var, Named("s") = sample);
          y(i,0) = sample;
          //a(i, 0) = mean;
          //a(i, 1) = var;
        }else{
          
          //Rcout << inv_sigma << std::endl;
          double var = std::pow(sigma, 2) * lambda / (lambda * nk[i] + 1);
          //return List::create(Named("var") = var, Named("parameters") = parameters, Named("inv_sigma") = inv_sigma, Named("nk")=nk);
          LogicalVector cluster_index = (cluster == i);
          double mean = (lambda * colSums(row_matrix(X, cluster_index))[0] + as<double>(parameters["mu"])) / (lambda * nk[i] + 1);
          
          double sample = rnorm(1, mean, sqrt(var))[0];
          //return List::create(Named("mean") = mean, Named("var") = var, Named("s") = sample);
          y(i,0) = sample;
          
          //a(i, 0) = mean;
          //a(i, 1) = var;
        }
        
      }
      // dreturn List::create(Named("a") = a,Named("y") = y, Named("nk") = nk, Named("pi") = pi);
    }else{
      //Rcout << y << std::endl;
      //NumericMatrix inv_sigma = solve(Sigma);
      for(int i = 0 ; i < y.nrow(); ++i){
        //Rcout <<i << " "<< nk[i] << std::endl;
        if(nk[i] == 0){
          //Rcout <<i << " "<< nk[i] << std::endl;
          NumericMatrix var = Sigma;
          NumericVector mean = parameters["mu"];
          NumericMatrix r_sample = wrap(rmvnorm(1, mean, lambda * as<mat>(var)));
          //NumericVector sample = (0, _);
          y(i,_) = r_sample(0, _);
        }else{
          //Rcout <<i << " "<< nk[i] << std::endl;
          //Rcout << inv_sigma << std::endl;
          NumericMatrix var = matrix_mul_scalar(Sigma, lambda / (lambda * nk[i] + 1));
          //NumericMatrix var = solve_pos_def(matrix_add(as<NumericMatrix>(parameters["inv_Sigma"]), matrix_mul_scalar(inv_sigma, nk[i])));
          //Rcout << var << std::endl;
          //var = make_symmetric(var);
          //return List::create(Named("var") = var, Named("parameters") = parameters, Named("inv_sigma") = inv_sigma, Named("nk")=nk);
          LogicalVector cluster_index = (cluster == i);
          //NumericMatrix sums_b_h = NumericMatrix(1, var.nrow(), colSums(row_matrix(X, cluster_index)).begin());
          NumericVector sums_b_h = colSums(row_matrix(X, cluster_index));
          
          NumericVector mean = (sums_b_h * lambda + as<NumericVector>(parameters["mu"])) / (lambda * nk[i] + 1);
            
            
          //NumericMatrix mu_par = NumericMatrix(1, var.nrow(), as<NumericVector>(parameters["mu"]).begin());
          //NumericVector mean = matrix_mul_scalar(Sigma, lambda / (lambda * nk[i] + 1));
          //NumericVector mean = matrix_multiply(matrix_add(matrix_multiply(sums_b_h, inv_sigma), matrix_multiply(mu_par, parameters["inv_Sigma"])), var);
          //NumericVector mean = matrix_multiply(matrix_add(matrix_multiply(sums_b_h, inv_sigma), matrix_multiply(mu_par, parameters["inv_Sigma"])), var);
          //var = solve(solve(tau$parameters$sigma) + nk[h] * solve(sigma))
          //mean = (colSums(X[tau$cluster == h,, drop = FALSE]) %*% solve(sigma) + tau$parameters$mu %*% solve(tau$parameters$sigma)) %*% var
          //return List::create(Named("rows") = colSums(row_matrix(X, cluster_index)), Named("inv_sigma") = inv_sigma, Named("mean") = mean, Named("var") = var, Named("parameters") = parameters);
          //Rcout << mean << var << std::endl;
          
          
          NumericMatrix r_sample = wrap(rmvnorm(1, mean, as<mat>(var)));
          //NumericVector sample = (0, _);
          y(i,_) = r_sample(0, _);
          
          //NumericVector neab = matrix_multiply(matrix_multiply(colSums(row_matrix(X, cluster_index)), inv_sigma) + matrix_multiply(parameters["mu"], parameters["inv_sigma"]), var);
        }
        
        //NumericVector sample = rmvnorm(1, mean, var);
        //y(i,_) = sample;
      }
    }
    
    if(CDP == true){
      for(int i = 0; i < p; ++i){
        y (_,i) = y (_,i) - sum(y(_,i) * pi);
      }
    }
    
    NumericMatrix samples(N, p);
    for(int i = 0; i < N; ++i){
      samples(i,_) = y(cluster[i], _);
    }
    
    
    tau["pi"] = pi;
    tau["M"] = M;
    tau["y"] = y;
    tau["cluster"] = cluster;
    tau["samples"] = samples;
    tau["lambda"] = lambda;
    parameters["lambda"] = lambda;
    
    if(parameters.containsElementNamed("sd")){
      double v = parameters["v"];
      double k = parameters["k"];
      NumericVector X_ = X(_, 0);
      NumericVector samples_ = samples(_, 0);
      NumericVector y_ = y(_, 0);
      tau["sigma"] = sqrt(rinvgamma((N + v + y.length()) / 2, (k*v + innerProduct(X_ - samples_, X_ - samples_) + innerProduct(y_ - as<double>(parameters["mu"]), y_ - as<double>(parameters["mu"])) / lambda)/2));
      parameters["sigma"] = tau["sigma"];
      //parameters["inv_sigma_2"] = 1 / pow(as<double>(tau["sigma"]), 2);
    }else{
      NumericMatrix Psi = parameters["Psi"];
      double d = parameters["d"];
      NumericMatrix X_diff_musi = matrix_add(X, matrix_mul_scalar(samples, -1));
      NumericMatrix B_B_A = matrix_multiply(transpose(X_diff_musi), X_diff_musi);
      //Rcout << B_B_A << std::endl;
      NumericVector mu0 = parameters["mu"];
      
      NumericMatrix musi_diff_mu0 = wrap(as<arma::mat>(y).each_row() - as<arma::rowvec>(mu0));
      //Rcout << y << std::endl;
      //Rcout << musi_diff_mu0 << std::endl;
      NumericMatrix mu_mu_A = matrix_multiply(transpose(musi_diff_mu0), musi_diff_mu0);
      //Rcout <<  mu_mu_A << std::endl;
      mu_mu_A = matrix_mul_scalar(mu_mu_A, 1 / lambda);
      //Rcout <<  mu_mu_A << std::endl;
      NumericMatrix Iwish_para = matrix_add(matrix_add(Psi, B_B_A), mu_mu_A);
      //Rcout << Iwish_para << std::endl;
      //Rcout << d + N + y.nrow() << std::endl;
      //Rcout << matrix_mul_scalar(Iwish_para, 1 / (d + N + y.nrow() - 4)) << std::endl;
      NumericMatrix covariance = wrap(riwishArma(d + N + y.nrow(), as<arma::mat>(Iwish_para)));
      tau["Sigma"] = covariance;
      parameters["Sigma"] = covariance;
      //parameters["inv_Sigma"] = solve_pos_def(covariance);
    }

    tau["parameters"] = parameters;
    return tau;
  }
  return tau;
}