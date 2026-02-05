/*
 *  SBMTrees: Sequential imputation with Bayesian Trees Mixed-Effects models
 *  Copyright (C) 2024 Jungang Zou
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy is available at
 *  https://www.R-project.org/Licenses/GPL-2
 */

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
  //Rcout << y << std::endl;
  std::string dist = as<std::string>(parameters["distribution"]);

  if(dist == "normal"){
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
    if(parameters.containsElementNamed("sd"))
      return List::create(Named("samples") = samples, Named("cluster") = cluster, Named("pi") = p, Named("parameters") = parameters, Named("y") = y, Named("M") = M, Named("lambda") = parameters["lambda"], Named("sigma") = parameters["sigma"]);
    else
      return List::create(Named("samples") = samples, Named("cluster") = cluster, Named("pi") = p, Named("parameters") = parameters, Named("y") = y, Named("M") = M, Named("lambda") = parameters["lambda"], Named("Sigma") = parameters["Sigma"]);
  }
  return List::create(Named("samples") = samples, Named("cluster") = cluster, Named("index") = index, Named("pi") = p, Named("parameters") = parameters, Named("y") = y, Named("lambda") = parameters["lambda"]);
}

// [[Rcpp::export]]
List update_DP_normal(NumericMatrix X, List tau, double L = 0, double U = 0.5){
  List parameters = tau["parameters"];
  double M = tau["M"];
  long N = X.nrow();
  if(strcmp(parameters["distribution"], "normal") == 0){
    int p = as<int>(parameters["p"]);
    double sigma = 1;
    NumericMatrix Sigma;
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
    if(parameters.containsElementNamed("sd")){
      lambda = rinvgamma(y.length() / 2 + as<double>(parameters["a"]), as<double>(parameters["b"]) + innerProduct(y - as<double>(parameters["mu"]), y - as<double>(parameters["mu"])) / 2 / pow(sigma, 2));
    }else{
      lambda = rinvgamma(y.length() / 2 + as<double>(parameters["a"]), as<double>(parameters["b"]) + quadratic_form(y, as<NumericVector>(parameters["mu"]), Sigma) / 2);
    }
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
   //Rcout << "complete update cluster probability" << std::endl;
   //Rcout << "update sampling probability pi" << std::endl;

    //Rcout << "update table count" << std::endl;
    IntegerVector nk = table(cluster);
    for(int i = 0; i < y.nrow(); ++i){
      std::string n_str = std::to_string(i);
      const char *n = n_str.c_str();
      if(!nk.containsElementNamed(n)){
        nk[n] = 0;
      }
    }


    CharacterVector p_names = rownames(y);
    nk = nk[p_names];

    //Rcout << "update beta" << std::endl;
    NumericMatrix beta(nk.length() - 1, 2);
    beta(_,0) = as<NumericVector>(nk[seqC(0, (nk.length() - 2))]);
    IntegerVector rev_nk = rev(nk);

    IntegerVector rev_cumsum = cumsum(rev_nk);
    IntegerVector non_rev_cumsum = rev(rev_cumsum);
    beta(_,1) = as<NumericVector>(non_rev_cumsum[seqC(1, (nk.length() - 1))]);
    NumericVector Vh(nk.length() - 1);
    for(int i = 0; i < nk.length() - 1; ++i){
      double new_beta_sample = R::rbeta(1 + beta(i, 0), M + beta(i, 1));
      int beta_count = 0;
      while(new_beta_sample == 1.0 && beta_count < 10){
        //Rcout <<"beta_count:" << beta_count << std::endl;
        beta_count += 1;
        new_beta_sample = R::rbeta(1 + beta(i, 0), M + beta(i, 1));
      }
      if(beta_count >= 10)
        new_beta_sample = (1 + beta(i, 0)) / (1 + beta(i, 0) + M + beta(i, 1));
      Vh[i] = new_beta_sample;
    }
    Vh.push_back(1);
    //Rcout << "finish update beta" << std::endl;

    //Rcout << "update pi" << std::endl;
    pi = NumericVector(nk.length());
    pi[0] = Vh[0];
    NumericVector log_negative_1_cumsum = cumsum(log(1 - Vh));
    //Rcout << negative_1_cumprod << std::endl;
    for(int i = 1; i < nk.length() - 1; ++i){
      pi[i] = Vh[i] * exp(log_negative_1_cumsum[i - 1]);
    }

    if (nk.length() > 1){
      //Rcout << "nk>1" << std::endl;
      pi[nk.length() - 1] = exp(log_negative_1_cumsum[nk.length() - 2]);
      M = as<NumericVector>(rtgamma(1, nk.length() - 1, -1 / log_negative_1_cumsum[nk.length() - 2], pow(N, L), pow(N, U)))[0];
    }else{
      pi[0] = 1;
      M = 0;
    }
    pi.names() = p_names;
    //Rcout << "complete update sampling probability pi" << std::endl;

    //return List::create(Named("Vh") = Vh, Named("nk") = nk, Named("pi") = pi);
    //Rcout << 123 << std::endl;
    //Rcout << p << std::endl;
    if(parameters.containsElementNamed("sd")){
      //Rcout << "update cluster mean" << std::endl;
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
      //Rcout << "complete update cluster mean" << std::endl;
      //Rcout << 123 << std::endl;
      // dreturn List::create(Named("a") = a,Named("y") = y, Named("nk") = nk, Named("pi") = pi);
    }else{
      //Rcout << y << std::endl;
      //NumericMatrix inv_sigma = solve(Sigma);
      //Rcout << "update cluster mean" << std::endl;
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
      //Rcout << "complete update cluster mean" << std::endl;
    }


    //Rcout << "update centralization" << std::endl;
    if(CDP == true){
      for(int i = 0; i < p; ++i){
        y (_,i) = y (_,i) - sum(y(_,i) * pi);
      }
    }
    //return List::create();
    NumericMatrix samples(N, p);
    for(int i = 0; i < N; ++i){
      samples(i,_) = y(cluster[i], _);
    }

    //Rcout << "complete update centralization" << std::endl;

    tau["pi"] = pi;
    tau["M"] = M;
    tau["y"] = y;
    tau["cluster"] = cluster;
    tau["samples"] = samples;
    tau["lambda"] = lambda;
    parameters["lambda"] = lambda;

    //Rcout << "update Sigma/sigma" << std::endl;

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
    //Rcout << "complete update Sigma/sigma" << std::endl;
    tau["parameters"] = parameters;
    return tau;
  }
  return tau;
}

// [[Rcpp::export]]
List resample_tau(List tau){
  IntegerVector cluster = tau["cluster"];
  NumericVector pi = tau["pi"];
  NumericMatrix y = tau["y"];
  NumericMatrix samples = tau["samples"];
  for(int i = 0; i < cluster.length(); ++i){

    //Rcout << seqC(0, pi.length()-1) << std::endl;
    int c = sample(seqC(0, pi.length()-1), 1, true, pi)[0];
    //Rcout << i << " " << cluster.length() << "   " << c << std::endl;
    cluster[i] = c;
    samples(i, _) = y(cluster[i], _);
  }
  //tau["cluster"] = cluster;
  //tau["samples"] = samples;

  return List::create(Named("cluster") = cluster, Named("samples") = samples);
}


// [[Rcpp::export]]
List resample_tau_prob(List tau, double correct_prob){
  IntegerVector cluster = tau["cluster"];
  NumericVector pi = tau["pi"];
  NumericMatrix y = tau["y"];
  NumericMatrix samples = tau["samples"];
  for(int i = 0; i < cluster.length(); ++i){
    if(pi[cluster[i]] < correct_prob){
      int c = sample(seqC(0, pi.length()-1), 1, true, pi)[0];
      //Rcout << i << " " << cluster.length() << "   " << c << std::endl;
      cluster[i] = c;
      samples(i, _) = y(cluster[i], _);
    }
  }
  //tau["cluster"] = cluster;
  //tau["samples"] = samples;

  return List::create(Named("cluster") = cluster, Named("samples") = samples);
}



// [[Rcpp::export]]
NumericMatrix predict_tau(List tau, long N = 1){
  IntegerVector nk = table(as<IntegerVector>(tau["cluster"]));
  CharacterVector cluster_name = nk.names();
  List parameters = tau["parameters"];
  NumericMatrix center = tau["y"];
  double n = as<NumericMatrix>(tau["samples"]).nrow();

  NumericMatrix predicted_value(N, as<long>(parameters["p"]));
  double M = tau["M"];
  for(long i = 0; i < N; ++i){
    //Rcout << i << ": " << N << std::endl;
    NumericVector sample_vec;
    if(runif(1)[0] <= (M / (M + n))){
      //Rcout << 455 << std::endl;
      if(parameters.containsElementNamed("sd")){
        sample_vec = R::rnorm(as<double>(parameters["mu"]), sqrt(as<double>(tau["lambda"])) * as<double>(tau["sigma"]));
      }else{
        sample_vec = rmvnorm(1, as<vec>(parameters["mu"]), as<double>(tau["lambda"]) * as<mat>(tau["Sigma"]));
      }
    }else{
      //Rcout << 1 << std::endl;
      NumericVector prob = as<NumericVector>(nk) / (as<double>(tau["M"]) + n);

      //Rcout << 12 << std::endl;
      std::string s = as<std::string>(sample(cluster_name, 1, true, prob)[0]);

      //Rcout << 123 << std::endl;
      sample_vec = center(std::stoi(s), _);
      //Rcout << sample_vec << std::endl;
    }
    predicted_value(i, _) = sample_vec;
  }

  return predicted_value;
}

// [[Rcpp::export]]
NumericVector probit_probability_tau(List tau, NumericVector y_input){
  NumericMatrix y = tau["y"];     // K x 1
  NumericVector pi = tau["pi"];   // length K

  int n = y_input.size();
  NumericVector prob(n, 0.0);

  for(int k = 0; k < y.nrow(); ++k){
    prob = prob + pi[k] * Rcpp::pnorm(y_input + y(k,0));
  }
  return prob;
}



// [[Rcpp::export]]
NumericVector normal_probability_tau(List tau, NumericVector y, double sigma){
  NumericVector mu = as<NumericVector>(tau["y"]);
  NumericVector pi = as<NumericVector>(tau["pi"]);
  NumericVector prob(y.length());
  for(int i = 0; i < pi.length(); ++i){
    //sRcout << as<NumericVector>(y + mu[i] << std::endl;
    prob = prob + pi[i] * Rcpp::dnorm(y, mu[i], sigma);
  }
  return prob;
}


// [[Rcpp::export]]
NumericVector normal_loglik_tau(List tau, NumericVector y, double sigma) {
  NumericVector mu = as<NumericVector>(tau["y"]);   // means
  NumericVector pi = as<NumericVector>(tau["pi"]);  // mixture weights
  int n = y.size();
  int K = pi.size();

  NumericVector loglik(n); // log density for each y[i]

  for (int i = 0; i < n; ++i) {
    // store log(pi_k) + log φ(y[i] | mu[k], sigma^2)
    NumericVector log_terms(K);
    for (int k = 0; k < K; ++k) {
      log_terms[k] = std::log(pi[k]) + R::dnorm(y[i], mu[k], sigma, true);
    }
    // log-sum-exp trick
    double m = max(log_terms);
    double sum_exp = 0.0;
    for (int k = 0; k < K; ++k) {
      sum_exp += std::exp(log_terms[k] - m);
    }
    loglik[i] = m + std::log(sum_exp);
  }

  return loglik; // vector of log-densities
}




// [[Rcpp::export]]
double normal_loglik_tau_scalar(List tau, double y, double sigma) {
  NumericMatrix ymat = tau["y"];   // K x 1
  NumericVector pi   = tau["pi"];  // length K
  int K = pi.size();

  NumericVector log_terms(K);
  for(int k = 0; k < K; ++k){
    if(pi[k] <= 0.0){
      log_terms[k] = R_NegInf;
    } else {
      log_terms[k] = std::log(pi[k]) + R::dnorm(y, ymat(k,0), sigma, true);
    }
  }

  double m = max(log_terms);
  if(!R_finite(m)) return m; // all -Inf case
  double sum_exp = 0.0;
  for(int k = 0; k < K; ++k){
    if(R_finite(log_terms[k])) sum_exp += std::exp(log_terms[k] - m);
  }
  return m + std::log(sum_exp);
}


