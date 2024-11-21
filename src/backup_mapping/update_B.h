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
#include <cmath>

using namespace Rcpp;

NumericMatrix update_B(NumericVector re, NumericMatrix Z, CharacterVector subject_id, NumericMatrix Mu, std::unordered_map<std::string, int> subject_to_B, NumericMatrix Covariance, double sigma){
  NumericMatrix B(subject_to_B.size(), Z.ncol());
                              
  CharacterVector subject_names = unique(subject_id);
  //Rcout << "inv_covariance" << std::endl;
  NumericMatrix inv_covariance = solve_pos_def(Covariance);
  //Rcout << as<NumericVector>(Mu) << std::endl;
  //NumericMatrix inv_covariance = solve_pos_def(make_symmetric(make_nonsingular(Covariance)));
  //NumericMatrix inv_covariance = solve(Covariance);
  //return inv_covariance;
  //Rcout << "B" << std::endl;
  //Rcout << inv_covariance << "  |  ";
  for(int i = 0; i < subject_names.length(); ++i){
    auto it = subject_to_B.find(std::string(subject_names[i]));
    if (it != subject_to_B.end()) {
      String name = subject_names[i];
      long B_position = (int)(it->second);
      //std::cout << name << " ";
      //Rcout << i << " " << B_position << std::endl;
      CharacterVector subject = {name};
      LogicalVector index = character_vector_equals(subject_id, subject);
      //return List::create(Named("index") = index);
      NumericMatrix Zi = row_matrix(Z, index);
      NumericVector Ri = re[index];
      //return List::create(Named("Zi") = Zi, Named("Ri") = Ri);
      
      //double alpha_s = alpha[B_position];
      //Rcout << inv_covariance;
      NumericMatrix a = matrix_add(inv_covariance, matrix_mul_scalar(matrix_multiply(transpose(Zi), Zi), 1 / pow(sigma, 2)));
      NumericMatrix var = solve_pos_def(a);
      
      NumericVector mui_v = Mu(B_position, _);
      NumericMatrix mui = NumericMatrix(1, mui_v.length(), mui_v.begin());
      NumericMatrix ri = NumericMatrix(Ri.length(), 1, Ri.begin());
      NumericVector mu = wrap(matrix_multiply( matrix_add(matrix_multiply(mui, inv_covariance), transpose(matrix_mul_scalar(matrix_multiply(transpose(Zi), ri), 1/ pow(sigma,2)))), var));
      
      //Rcout << mu << " | " << var(0, 0) << " ||| ";
      
      NumericVector B_s = wrap(rmvnorm(1, mu, as<arma::mat>(var)));
      B(B_position, _) = B_s;
    }else{ 
      
    }
    
    
  //for(int i = 0; i < 1; ++i){
    
  }
  //Rcout <<std::endl;
  return B;
}