#ifndef ARMADILLO_H_
#define ARMADILLO_H_
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#endif

#ifndef DIST_H_
#define DIST_H_
#include <RcppDist.h>
// [[Rcpp::depends(RcppArmadillo, RcppDist)]]
#endif

#include <Rcpp.h>
#ifndef UTILS_H_
#define UTILS_H_
#include "utils.h"
#endif
#include <cmath>
#ifndef DP_SAMPLER_H_
#define DP_SAMPLER_H_
#include "DP_sampler.h"
#endif

using namespace Rcpp;

// [[Rcpp::export]]
List DP(List parameters, double M, long N_truncated, long N_sample, bool CDP = true){
  
  NumericVector p(N_truncated);
  IntegerVector index = seq(0, N_truncated - 1);
  IntegerVector cluster;
  NumericMatrix samples;
  if(N_truncated > 1){
    NumericVector b = rbeta(N_truncated, 1, M);
    NumericVector index(N_truncated);
    p[0] = b[0];
    for(int i = 1; i <= N_truncated-1; i++){
      p[i] = b[i] * prodC(1.0 - as<NumericVector>(b[seqC(0,i-1)]));
    }
  }else{
    p[0] = 1;
  }
  p.attr("names") = index;
  List samples_parameters = DP_sampler(N_truncated, parameters);
  NumericMatrix y = samples_parameters["y"];
  
  if(strcmp(parameters["distribution"], "normal") == 0){
    int dp = as<int>(parameters["p"]);
    parameters["CDP"] = CDP;
    CharacterVector rowname(index.begin(), index.end());
    rownames(y) = rowname;
    cluster = sample(index, N_sample, true, p);
    if(CDP == true){
      for(int i = 0; i < dp; ++i){
        y (_,i) = y (_,i) - sum(y(_,i) * p);
      }
    }
    NumericMatrix samples(N_sample, dp);
    for(int i = 0; i < N_sample; ++i){
      samples(i,_) = y(cluster[i], _);
    }
    rownames(samples) = cluster;
    return List::create(Named("samples") = samples, Named("cluster") = cluster, Named("pi") = p, Named("parameters") = parameters, Named("y") = y, Named("M") = M);
  }
  return List::create(Named("samples") = samples, Named("cluster") = cluster, Named("index") = index, Named("pi") = p, Named("parameters") = parameters, Named("y") = y);
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
    if(p == 1){
      sigma = tau["sigma"];
    }else{
      Sigma = as<NumericMatrix>(tau["Sigma"]);
    }
    
    
    IntegerVector cluster = tau["cluster"];
    NumericVector pi = tau["pi"];
    NumericMatrix y = tau["y"];
    bool CDP = parameters["CDP"];
    NumericVector log_density(y.nrow());
    NumericVector log_dens;
    //Rcout << 123 << std::endl;
    //Rcout << density << std::endl;
    for(int i = 0 ; i < X.rows(); ++i){
      NumericVector X_row = X(i, _);
      
      for(int j = 0 ; j < y.nrow(); ++j){
        NumericVector mu = y(j,_);
        //Rcout << mu[0] << std::endl;
        double pj = pi[j];
        if(p == 1){
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
    if(p == 1){
      //NumericMatrix a(y.nrow(), 2);
      double inv_sigma_2 = std::pow(1 / sigma, 2);
      for(int i = 0 ; i < y.nrow(); ++i){
        if(nk[i] == 0){
          double var = std::pow(as<double>(parameters["sigma"]), 2);
          double mean = parameters["mu"];
          
          double sample = rnorm(1, mean, sqrt(var))[0];
          //return List::create(Named("mean") = mean, Named("var") = var, Named("s") = sample);
          y(i,0) = sample;
          //a(i, 0) = mean;
          //a(i, 1) = var;
        }else{
          
          //Rcout << inv_sigma << std::endl;
          double var = 1 / (as<double>(parameters["inv_sigma_2"]) + inv_sigma_2 * nk[i]);
          //return List::create(Named("var") = var, Named("parameters") = parameters, Named("inv_sigma") = inv_sigma, Named("nk")=nk);
          LogicalVector cluster_index = (cluster == i);
          double mean = ((colSums(row_matrix(X, cluster_index))[0] * inv_sigma_2) + (as<double>(parameters["mu"]) * as<double>(parameters["inv_sigma_2"]))) * var;
          
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
      NumericMatrix inv_sigma = solve(Sigma);
      for(int i = 0 ; i < y.nrow(); ++i){
        //Rcout <<i << " "<< nk[i] << std::endl;
        if(nk[i] == 0){
          //Rcout <<i << " "<< nk[i] << std::endl;
          NumericMatrix var = parameters["Sigma"];
          NumericVector mean = parameters["mu"];
          NumericMatrix r_sample = wrap(rmvnorm(1, mean, as<mat>(var)));
          //NumericVector sample = (0, _);
          y(i,_) = r_sample(0, _);
        }else{
          //Rcout <<i << " "<< nk[i] << std::endl;
          //Rcout << inv_sigma << std::endl;
          NumericMatrix var = solve_pos_def(matrix_add(as<NumericMatrix>(parameters["inv_Sigma"]), matrix_mul_scalar(inv_sigma, nk[i])));
          //Rcout << var << std::endl;
          //var = make_symmetric(var);
          //return List::create(Named("var") = var, Named("parameters") = parameters, Named("inv_sigma") = inv_sigma, Named("nk")=nk);
          LogicalVector cluster_index = (cluster == i);
          NumericMatrix sums_b_h = NumericMatrix(1, var.nrow(), colSums(row_matrix(X, cluster_index)).begin());
          
          NumericMatrix mu_par = NumericMatrix(1, var.nrow(), as<NumericVector>(parameters["mu"]).begin());
          
          //Rcout << sums_b_h <<std::endl;
          NumericVector mean = matrix_multiply(matrix_add(matrix_multiply(sums_b_h, inv_sigma), matrix_multiply(mu_par, parameters["inv_Sigma"])), var);
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
    return tau;
  }
  return tau;
}