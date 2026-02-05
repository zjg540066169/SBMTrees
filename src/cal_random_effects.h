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

#ifndef UTILS_H_
#define UTILS_H_
#include "utils.h"
#endif

#ifndef UPDATE_H_
#define UPDATE_H_
#include "DP_lambda.h"
#include <cmath>
#endif


#include <Rcpp.h>
#include <cmath>

using namespace Rcpp;
using namespace arma;

NumericVector cal_random_effects(
    NumericMatrix Z,
    CharacterVector subject_id,
    NumericMatrix B,
    std::unordered_map<std::string, int> subject_to_B,
    NumericMatrix Covariance,
    bool binary = false
){
  long N = subject_id.length();
  NumericVector re(N);
  arma::mat cov = as<arma::mat>(Covariance);

  // Cache draws for new subjects (temporary, per call)
  std::unordered_map<std::string, NumericVector> new_subject_B;

  for(long i = 0; i < N; ++i){
    std::string sid = Rcpp::as<std::string>(subject_id[i]);
    NumericVector zi = Z(i, _);

    // in-sample subject: use fitted Bi
    auto it_train = subject_to_B.find(sid);
    if(it_train != subject_to_B.end()){
      int b_pos = it_train->second;
      NumericVector Bi = B(b_pos, _);
      re[i] = innerProduct(zi, Bi);
      continue;
    }

    // out-of-sample subject: sample once per subject, reuse across rows
    auto it_new = new_subject_B.find(sid);
    if(it_new == new_subject_B.end()){
      NumericVector mean0(zi.size());              // zeros
      NumericMatrix Bdraw = wrap(rmvnorm(1, mean0, cov));
      NumericVector Bi = Bdraw(0, _);
      new_subject_B.emplace(sid, Bi);
      re[i] = innerProduct(zi, Bi);
    } else {
      NumericVector Bi = it_new->second;
      re[i] = innerProduct(zi, Bi);
    }
  }

  return re;
}

NumericVector cal_random_effects_tau(
    NumericMatrix Z,
    CharacterVector subject_id,
    NumericMatrix B,
    std::unordered_map<std::string, int> subject_to_B,
    List B_tau,
    NumericMatrix Covariance,
    bool binary = false
){
  long N = subject_id.length();
  NumericVector re(N);
  arma::mat cov = as<arma::mat>(Covariance);

  // Temporary cache for *new* subjects in this prediction call:
  // key: subject id, value: sampled Bi vector
  std::unordered_map<std::string, NumericVector> new_subject_B;

  for(long i = 0; i < N; ++i){
    std::string sid = Rcpp::as<std::string>(subject_id[i]);
    NumericVector zi = Z(i, _);

    // 1) If subject is in training (has stored B), use it
    auto it_train = subject_to_B.find(sid);
    if(it_train != subject_to_B.end()){
      int b_pos = it_train->second;
      NumericVector Bi = B(b_pos, _);
      re[i] = innerProduct(zi, Bi);
      continue;
    }

    // 2) Otherwise, it's a new subject: reuse the same Bi within this call
    auto it_new = new_subject_B.find(sid);
    if(it_new == new_subject_B.end()){
      // First time we see this new subject -> sample Bi once
      NumericVector B_mu = predict_tau(B_tau, 1)(0, _);     // draw cluster mean for B
      NumericMatrix Bdraw = wrap(rmvnorm(1, B_mu, cov));    // draw Bi | mu
      NumericVector Bi = Bdraw(0, _);

      // store it
      new_subject_B.emplace(sid, Bi);

      re[i] = innerProduct(zi, Bi);
    } else {
      // reuse cached Bi
      NumericVector Bi = it_new->second;
      re[i] = innerProduct(zi, Bi);
    }
  }

  return re;
}



// NumericVector cal_random_effects_tau(NumericMatrix Z, CharacterVector subject_id, NumericMatrix B, std::unordered_map<std::string, int> subject_to_B, List B_tau, NumericMatrix Covariance, bool binary = false){
//
//   long N = subject_id.length();
//   NumericVector re(N);
//   mat cov = as<mat>(Covariance);
//   //Rcout << 123 << std::endl;
//   for(int i = 0; i < N; ++i){
//     String subject = subject_id[i];
//     NumericVector Bi;
//     //NumericVector B_mu = predict_tau(B_tau, 1)(0,_);
//     //Bi = rmvnorm(1, as<vec>(B_mu), as<mat>(Covariance));
//
//     auto it = subject_to_B.find(std::string(subject_id[i]));
//     if (it != subject_to_B.end()) {
//       int b_pos = (int)(it->second);
//       Bi = B(b_pos,_);
//       NumericVector zi = Z(i,_);
//       double zi_Bi = innerProduct(zi, Bi);
//       re[i] = zi_Bi;
//     }else{
//       NumericVector zi = Z(i,_);
//       NumericVector B_mu = predict_tau(B_tau, 1)(0,_);
//       //NumericMatrix B = wrap(rmvnorm(1, as<vec>(NumericVector(zi.length())), cov));
//       NumericMatrix B = wrap(rmvnorm(1, B_mu, cov));
//       NumericVector Bi = B(0,_);
//       double zi_Bi = innerProduct(zi, Bi);
//       re[i] = zi_Bi;
//     }
//
//
//     //
//   }
//   //Rcout << 1213 << std::endl;
//   return re;
// }







// // [[Rcpp::export]]
// vec cal_random_effects_arma(mat Z, CharacterVector subject_id, mat B, List subject_to_B){
//
//   long N = subject_id.length();
//   vec re(N);
//   //Rcout << 123 << std::endl;
//   for(int i = 0; i < N; ++i){
//     //Rcout << 456 << std::endl;
//     String subject = subject_id[i];
//     int b_pos = subject_to_B[subject];
//     vec zi = Z.row(i);
//     vec Bi = B.row(b_pos);
//     double zi_Bi = dot(zi, Bi);
//     re[i] = zi_Bi;
//   }
//   //Rcout << 1213 << std::endl;
//   return re;
// }
