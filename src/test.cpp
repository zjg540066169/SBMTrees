#include <RcppArmadillo.h>
#include <omp.h>
#include <RcppDist.h>
#include "DP_lambda_arma.h"

// [[Rcpp::depends(RcppArmadillo, RcppDist)]]
// [[Rcpp::export]]
Rcpp::List parallelMCMC(int n_chains, int n_dim) {
  Rcpp::List results(n_chains);
  

  
  // Set the number of threads (adjust based on your machine's capabilities)
  omp_set_num_threads(4);
  
  // Run parallel MCMC sampling
#pragma omp parallel
{
#pragma omp for
  for (int i = 0; i < n_chains; i++) {
        List par = List::create(Named("p") = n_dim);
        results[i] = DP_arma(par, 2, 100, 10000, true);
    } 
}

return results;
}
