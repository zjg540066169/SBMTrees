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

#ifndef UTILS_H_
#define UTILS_H_
#include "utils.h"
#endif

#ifndef RCPP_H_
#define RCPP_H_
#include <Rcpp.h>
#endif

#ifndef BART_M_H_
#define BART_M_H_
#include "bart_model.h"
#endif


#ifndef BM_H_
#define BM_H_
#include "bm.h"
#endif


#ifndef UPDATE_H_
#define UPDATE_H_
#include "DP_lambda.h"
#include "cal_random_effects.h"
#include "update_B.h"
#include "update_Covariance.h"
#include <cmath>
#endif

using namespace Rcpp;

class bmtrees : public bm{
public:
  ~bmtrees() override { delete tree; }

  bmtrees(NumericVector Y, NumericMatrix X, NumericMatrix Z, CharacterVector subject_id, IntegerVector row_id, bool binary = false, bool CDP_residual = false, bool CDP_re = false, double tol=1e-40, int ntrees = 200, double k = 2.0, double pi_CDP = 0.99, bool train = true) : bm(Y, X, Z, subject_id, row_id, binary, CDP_residual, CDP_re, tol, ntrees, k, pi_CDP, train) {     // Constructor
    tree = nullptr;
    if(train){
      tree = new bart_model(this->X, this->Y - re - tau_samples, 100L, false, false, false,  ntrees, R_NilValue, k);
      if(binary){
        tree -> update(1, 50, 50, 1, false, 10L);
      }else
        tree -> update(50, 50, 1, false, 10L);

      tree_pre = colMeans(tree->predict(this->X));
      sigma = tree->get_sigma();
    }
  }

  void update_sigma() override{
    NumericVector Y_hat = re + tau_samples + tree_pre;
    if(!binary){
      double rss = sum(pow(this->Y - Y_hat, 2));
      if(!CDP_residual)
        sigma = tree -> get_invchi(N, rss);
      else
        sigma = tau["sigma"];
    }else{
      sigma = 1;
    }
  }

  void update_tree(){

    NumericVector Y_ = Y - re - tau_samples;
    tree->set_data(X, Y_);
    List tree_obj = tree -> update(sigma, 1, 1, 1, false, 1L);
    tree_pre = tree_obj["yhat.train.mean"];
    // if(CDP_re || CDP_residual)
    //   tree_pre_mean = mean(tree_pre);
    // else
    //   tree_pre_mean = 0;
    // tree_pre = tree_pre - tree_pre_mean;
  }

  void set_tree(List tree){
    //this->tree = tree;
    tree_pre = tree["yhat.train.mean"];

    sigma = 1;
    if(!binary){
      sigma = tree["sigma.mean"];
    }
  }

  List get_tree(){
    return this->tree -> get_tree_object();;
  }

  void update_all(bool verbose = false) override{
    if(verbose)
      Rcout << "update BART" << std::endl;

    update_tree();
    if(verbose)
      Rcout << "update residual" << std::endl;
    NumericVector residual_tem = Y - re - tree_pre;
    NumericMatrix residual(N, 1, residual_tem.begin());

    if(CDP_residual){
      if(verbose)
        Rcout << "update DP residual" << std::endl;
      tau["sigma"] = sigma;
      tau =  update_DP_normal(residual, tau, 0, 0.5);
      tau_samples = tau["samples"];
      M = tau["M"];
      sigma = tau["sigma"];
      if(binary)
        sigma = 1;
    }else{
      if(verbose)
        Rcout << "update sigma" << std::endl;
      update_sigma();
    }

    if(CDP_re){
      if(verbose)
        Rcout << "update DP random effects" << std::endl;
      B_tau["Sigma"] = Covariance;
      if(d == 1){
        B_tau["sigma"] = Covariance;
      }

      B_tau = update_DP_normal(B, B_tau, 0, 0.5);
      B_tau_samples = as<NumericMatrix>(B_tau["samples"]);
      M_re = B_tau["M"];
      Covariance = as<NumericMatrix>(B_tau["Sigma"]);

    }else{
      Covariance = update_Covariance(clone(B), clone(B_tau_samples), inverse_wishart_matrix, d + 2, n_subject);
    }
    cov = as<arma::mat>(Covariance);
    if(verbose)
      Rcout << "M_re:" << M_re << "  " << "M:" << M << std::endl;

    if(verbose)
      Rcout << "update B" << std::endl;
    B = update_B(Y - tau_samples - tree_pre, z, subject_id, B_tau_samples, subject_to_B, Covariance, sigma);

    if(verbose)
      Rcout << "update random effects" << std::endl;
    re = cal_random_effects(z, subject_id, B, subject_to_B, Covariance, binary);

    if(binary){
      for(int i = 0; i < N; ++i){
        if(Y_original[i] == 0){
          NumericVector mean_y = rtruncnorm(1, tree_pre[i] + re[i] + tau_samples[i], sigma, R_NegInf, 0);
          Y[i] = mean_y[0];
        }else{
          NumericVector mean_y = rtruncnorm(1, tree_pre[i] + re [i] + tau_samples[i], sigma, 0, R_PosInf);
          Y[i] = mean_y[0];
        }
      }
      if(verbose)
        Rcout << "update probit outcome" << std::endl;
    }
    row_id_to_re = create_row_id_to_re(row_id, re);
    new_subject_B.clear();
  }


  List posterior_sampling() override{
    return List::create(
      Named("tree") = tree->get_tree_object(),
      Named("M") = M,
      Named("M_re") = M_re,
      Named("sigma") = sigma,
      Named("Sigma") = Covariance,
      Named("B") = B,
      Named("tau_samples") = tau_samples,
      Named("B_tau_samples") = B_tau_samples,
      Named("re") = re,
      Named("tree_pre") = tree_pre,
      Named("y_predict") = tree_pre + re + tau_samples,
      Named("tau") = tau,
      Named("B_tau") = B_tau
    );
  }

  NumericVector predict_fix(NumericMatrix X_test) override{
    NumericVector X_hat_test = colMeans(tree -> predict(X_test, false));
    return X_hat_test;
  };



private:

  bart_model* tree;

  NumericVector tree_pre;

};

// // [[Rcpp::export]]
// List test(NumericVector Y, NumericMatrix X, NumericMatrix Z, CharacterVector subject_id, IntegerVector row_id, bool binary, bool nCDP_residual = false, bool CDP_re = false){
//   //Rcout << "123" << std::endl;
//   bmtrees a = bmtrees(Y, X, clone(Z), subject_id, row_id, binary, nCDP_residual, CDP_re);
//   a.update_all(true);
//   //a.update_all();
//   //a.update_all();
//   //Rcout << "456" << std::endl;
//   //a.update_tree();
//   //a.update_sigma();
//   //a.update(true);
//   List samples = a.posterior_sampling();
//   NumericVector expect_Y = as<NumericVector>(samples["tree_pre"]) + as<NumericVector>(samples["re"]);
//   return List::create(Named("sample") = a.posterior_sampling(), Named("X_hat_test") = a.predict(X, clone(Z), subject_id, row_id),Named("X_hat_test2") = a.predict(X, clone(Z), subject_id, row_id), Named("Y") = samples["y_predict"], Named("expect_Y") = expect_Y);
//   //return a.posterior_sampling();
//   //List post = a.posterior_sampling();
//   //NumericVector predict_y = a.predict(X, Z, subject_id);
//   //return a.posterior_sampling();
//   //return a.a_test;
//   //return List::create(Named("sample") = a.posterior_sampling(), Named("re") = a.get_re(), Named("samples") = a.get_tau_samples(), Named("Y") = a.get_Y());
// }
