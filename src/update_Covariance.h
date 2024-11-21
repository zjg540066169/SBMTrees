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

using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix update_Covariance(NumericMatrix B, NumericMatrix Mu, NumericMatrix inverse_wishart_matrix, double df, long N_subject){
  NumericMatrix B_diff_mu = matrix_add(B, matrix_mul_scalar(Mu, -1));
  NumericMatrix B_B_A = matrix_multiply(transpose(B_diff_mu), B_diff_mu);
  //Rcout <<  B_B_A  << std::endl;
  NumericMatrix Iwish_para = matrix_add(inverse_wishart_matrix, B_B_A);
  //Rcout << Iwish_para << std::endl;
  //Rcout <<  arma::chol(as<arma::mat>(Iwish_para)) << std::endl;
  //Iwish_para = fix_riwish(Iwish_para);
  //Iwish_para = make_symmetric(Iwish_para);
  NumericMatrix covariance = wrap(riwishArma(df + N_subject, as<arma::mat>(Iwish_para)));
  //Rcout << covariance << std::endl;
  //return fix_riwish(covariance);
  return covariance;
}