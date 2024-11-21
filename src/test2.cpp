#include <Rcpp.h>

#ifndef PARALLELFOR_H_
#define PARALLELFOR_H_
// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>
#endif

#ifndef BMTREES_H_
#define BMTREES_H_
#include "bmtrees.h"
#endif

#ifndef UPDATE_H_
#define UPDATE_H_
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
#endif

#include <vector>
#include <thread>
#include <algorithm> // for std::for_each
#include <mutex>
#include <condition_variable>
#include <queue>

using namespace Rcpp;
using namespace RcppParallel;


/*** R
# predict.DP_LMM_BART = function(chain, X_test, Z_test, subject_id_test){
#   
#   if(is.null(dim(Z_test))){d = 1} else{d = dim(Z_test)[2]}
#   subject_to_B = chain$subject_to_B
#   
#   
#  invisible(capture.output( X_hat <-  colMeans(predict(chain$post_trees, unname(as.matrix(X_test))))))
#   
#   
#   B = matrix(colMeans(chain$post_B), ncol = d)
#   
#   re = sapply(1:length(subject_id_test), function(i){
#     id = subject_id_test[i]
#     pos = unlist(subject_to_B[as.character(id)]) + 1
#     Bi = B[pos, ]
#     Z_test[i,] %*% Bi
#   })
#   e = colMeans(chain$post_tau_samples)
#   
#   return(unname(X_hat + re + e))
# }

library(foreach)
library(doParallel)
library(Rcpp)
#sourceCpp("./DP_LMM_BART2.cpp")
# model_fitting = function(chain_collection, X, Y, Z, subject_id, type, R, ncores, seed = NULL){
#   registerDoParallel(ncores)
#   chain_collection = foreach(i = 1:(dim(X)[2]), .packages = c("Rcpp", "dbarts", "mvtnorm", "stats", "truncnorm", "matrixcalc", "pedmod"), .export = c("chain_collection")) %dopar% {
#     # fit outcome model
#       if(i == p){
#         chain = chain_collection[i]
#         chain.update_X(X)
#         chain.update_Y(Y)
#         chain.update(FALSE, seed)
#         return(chain)
#       }
#       if(sum(R(_, i + 1)) != 0){
#         X_train = X[,1:i]
#         y_train = X[,i+1]
#         chain = chain_collection[i]
#         chain.update_X(X_train)
#         chain.update_Y(y_train)
#         chain.update(FALSE, seed)
#         return(chain)
#       }
#     return(chain)
#   }
#   stopImplicitCluster()
#   #cat("finish fitting")
#   return(chain_collection)
# }







library(stats)
library(BART3)
library(stats)
library(MCMCpack)
library(mvtnorm)
library(cascsim)
library(truncnorm)
library(hbmem)

matrix_multiply = function(a, b){
  return(a%*%b)
}

matrix_mul_scalar = function(a, b){
  return(a * b)
}

matrix_add = function(a, b){
  return(a + b)
}

matrix_minus = function(a, b){
  return(a - b)
}

vector_mul_generate_matrix = function(A){
  A %*% t(A)
}

make_symmetric = function(s){
  s[lower.tri(s)] = t(s)[lower.tri(s)]
  return(s)
}

update_tree = function(X, y){
  #print(y)
  #print(class(X))
  #X = read_csv("X.csv")
  #y = read_csv("y.csv")
  
  if(sd(X[,1]) == 0 & dim(X)[2] > 1){
    X = X[,2:dim(X)[2]]
  }
  
  tree = gbart(X, y, verbose = 0, ndpost=1, nskip = 1)
  return(tree)
  #tree_pre = tree$yhat.train.mean#bart_fit[[1]]
}

predict_tree = function(trees, X_test){
  #print(X_test)
  if(sd(X_test[,1]) == 0 & dim(X_test)[2] > 1){
    X_test = X_test[,2:dim(X_test)[2]]
  }
  #print(X_test)
  invisible(capture.output( X_hat <-  colMeans(predict(trees, unname(as.matrix(X_test))))))
  X_hat
}


matrix_multiply = function(a, b){
  return(a%*%b)
}

matrix_add = function(a, b){
  return(a + b)
}

vector_mul_generate_matrix = function(A){
  if(class(A) == "numeric"){
    return(A %*% t(A))
  }else{
    return(t(A) %*% (A))
  }
  
}

make_symmetric = function(s){
  s[lower.tri(s)] = t(s)[lower.tri(s)]
  return(s)
}

make_nonsingular = function(s){
  d = dim(s)[1]
  tol = 1e-30^(1/d)
  s + diag(d) * tol 
}

inner_prod = function(a, b){
  t(a) %*% b
}

tree_update = function(trees, data, ncores, seed = NULL){
  if(!is.null(seed))
    set.seed(seed)
  registerDoParallel(ncores)
  trees = foreach(i = 1:(length(trees)), .multicombine = T, .combine = list) %dopar% {
    X = data[[i]]$X
    y = data[[i]]$Y
    if(sd(X[,1]) == 0 & dim(X)[2] > 1){
      X = X[,2:dim(X)[2]]
    }
    tree = gbart(X, y, verbose = 0, ndpost=1, nskip = 1)
    #print(y)
    return(tree)
  }
  
  stopImplicitCluster()
  return(trees)
}

*/

// Thread pool class
class ThreadPool {
public:
  ThreadPool(size_t num_threads) : stop(false) {
    for (size_t i = 0; i < num_threads; ++i) {
      workers.emplace_back([this] {
        while (true) {
          std::function<void()> task;
          {
            std::unique_lock<std::mutex> lock(queue_mutex);
            condition.wait(lock, [this] { return stop || !tasks.empty(); });
            if (stop && tasks.empty()) {
              return;
            }
            task = std::move(tasks.front());
            tasks.pop();
          }
          task();
        }
      });
    }
  }
  
  template<class F>
  void enqueue(F&& f) {
    {
      std::unique_lock<std::mutex> lock(queue_mutex);
      tasks.emplace(std::forward<F>(f));
    }
    condition.notify_one();
  }
  
  ~ThreadPool() {
    {
      std::unique_lock<std::mutex> lock(queue_mutex);
      stop = true;
    }
    condition.notify_all();
    for (std::thread& worker : workers) {
      worker.join();
    }
  }
  
private:
  std::vector<std::thread> workers;
  std::queue<std::function<void()>> tasks;
  std::mutex queue_mutex;
  std::condition_variable condition;
  bool stop;
};


// [[Rcpp::export]]
void sequential_imputation(NumericMatrix X, NumericVector Y, LogicalVector type, NumericMatrix Z, CharacterVector subject_id, bool verbose = true, bool nCDP_residual = true, bool CDP_re = true, Nullable<long> seed = R_NilValue, double tol = 1e-20, int  ncores = 0) {
  Rcpp::Environment base("package:base");
  Rcpp::Environment G = Rcpp::Environment::global_env();
  //Rcpp::Function predict_dbarts_residual = G["predict.DP_LMM_BART"];
  //Rcpp::Function model_fitting = G["model_fitting"];
  Rcpp::Environment stats("package:stats");
  Rcpp::Function dnorm_cpp = stats["dnorm"];
  Rcpp::Function pnorm_cpp = stats["pnorm"];
  if(seed.isNotNull()){
    Rcpp::Function set_seed_r = base["set.seed"];
    set_seed_r(seed);
    //Rcout << seed;
  }
  
  int n = X.nrow();
  int p = X.cols();
  List imputation_X_DP = List::create();
  List imputation_Y_DP = List::create();
  int  skip_indicator = -1;
  //List chain_collection = List::create();
  std::vector<bmtrees> chain_collection; 
  
  
  CharacterVector X_names = colnames(X);
  //Rcout << X_names << std::endl;
  for(int i = 0; i < p; ++i){
    if(i == p - 1){
      chain_collection.emplace_back(bmtrees(clone(Y), clone(X), clone(Z), clone(subject_id), 0, nCDP_residual, CDP_re, tol, seed));
      break;
    }
    NumericMatrix X_train = X(_, Range(0,i));
    NumericVector y_train = X(_, i + 1);
    chain_collection.emplace_back(bmtrees(clone(y_train), clone(X_train), clone(Z), subject_id, type[i+1], nCDP_residual, CDP_re, tol, seed));
  }
  
  if (verbose){
    Rcout << "Complete initialization" << std::endl;
    Rcout << std::endl;
  }
  
  Rcpp::Function tree_update = G["tree_update"];
  List tree_list;
  List tree_train_data_list;
  for(int i = 0; i < p; ++i){
    //Rcout << i << p << std::endl;
    if(i == p - 1){
      chain_collection[i].update_X(clone(X));
      tree_list.push_back(chain_collection[i].get_tree());
      tree_train_data_list.push_back(chain_collection[i].get_tree_training_data());
    }else{
      
        NumericMatrix X_train = X(_, Range(0,i));
        NumericVector y_train = X(_, i + 1);
        chain_collection[i].update_X(clone(X_train));
        chain_collection[i].update_Y(clone(y_train));
        tree_list.push_back(chain_collection[i].get_tree());
        tree_train_data_list.push_back(chain_collection[i].get_tree_training_data());
      
    }
  }
  tree_list = tree_update(tree_list, tree_train_data_list, ncores, seed);
  int flag = 0;
  for(int i = 0; i < p; ++i){
    //Rcout << i << p << std::endl;
    if(i == p - 1){
      chain_collection[i].set_tree(tree_list[flag]);
      flag ++;
    }else{
      
        chain_collection[i].set_tree(tree_list[flag]);
        flag ++;
      
    }
  }
  
  
  
  unsigned int num_cores = ncores;
  ThreadPool pool(num_cores);
  
  for (auto& obj : chain_collection) {
    pool.enqueue([&obj] {
      obj.update(true);
    });
  }
 


}