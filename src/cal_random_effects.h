#ifndef ARMADILLO_H_
#define ARMADILLO_H_
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#endif

#ifndef UTILS_H_
#define UTILS_H_
#include "utils.h"
#endif



#include <Rcpp.h>
#include <cmath>

using namespace Rcpp;
using namespace arma;

NumericVector cal_random_effects(NumericMatrix Z, CharacterVector subject_id, NumericMatrix B, std::unordered_map<std::string, int> subject_to_B){
  
  long N = subject_id.length();
  NumericVector re(N);
  //Rcout << 123 << std::endl;
  for(int i = 0; i < N; ++i){
    String subject = subject_id[i];
    auto it = subject_to_B.find(std::string(subject_id[i]));
    if (it != subject_to_B.end()) {
    
    //Rcout << 456 << std::endl;
    
    int b_pos = (int)(it->second);
    NumericVector zi = Z(i,_);
    NumericVector Bi = B(b_pos,_);
    double zi_Bi = innerProduct(zi, Bi);
    re[i] = zi_Bi;
    }
  } 
  //Rcout << 1213 << std::endl;
  return re;
}

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