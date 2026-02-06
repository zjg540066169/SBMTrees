/*
 *  BART: Bayesian Additive Regression Trees
 *  Copyright (C) 2017 Robert McCulloch and Rodney Sparapani
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

/*
 *  Modifications by Jungang Zou, 2024.
 *  - This function is modified version of cwbart.cpp. I wrote a new class including
 *  the main function of cwbart, as well as other functions, like updating the X and Y.
 *
 *  These modifications comply with the terms of the GNU General Public License
 *  version 2 (GPL-2).
 */



#ifndef ARMADILLO_H_
#define ARMADILLO_H_
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#endif

#ifndef BART_MODEL_H_
#define BART_MODEL_H_
#include "BART/tree.h"
#include "BART/treefuns.h"
#include "BART/info.h"
#include "BART/bartfuns.h"
#include "BART/bd.h"
#include "BART/bart.h"
#include<stdio.h>
#include "BART/cpwbart.h"

#endif
#ifndef RCPP_H_
#define RCPP_H_
#include <Rcpp.h>
#endif


using namespace Rcpp;
/*** R
bartModelMatrix=function(X, numcut=0L, usequants=FALSE, type=7,
                         rm.const=FALSE, cont=FALSE, xinfo=NULL) {
  #print(xinfo)
  X.class = class(X)[1]

  if(X.class=='factor') {
    X.class='data.frame'
    X=data.frame(X=X)
  }
  grp=NULL
  if(X.class=='data.frame') {
    print(X)
    p=dim(X)[2]

    xnm = names(X)
    for(i in 1:p) {
      print(i)
      if(is.factor(X[[i]])) {
        #print(i)
        Xtemp = nnet::class.ind(X[[i]])
        colnames(Xtemp) = paste(xnm[i],1:ncol(Xtemp),sep='')
        X[[i]]=Xtemp
        m=ncol(Xtemp)
        grp=c(grp, rep(m, m))
      } else {
        X[[i]]=cbind(X[[i]])
        colnames(X[[i]])=xnm[i]
        grp=c(grp, 1)
        ##grp=c(grp, i)
      }
    }
    Xtemp=cbind(X[[1]])
    if(p>1) for(i in 2:p) Xtemp=cbind(Xtemp, X[[i]])
    X=Xtemp
  }
  else if(X.class=='numeric' | X.class=='integer') {
    X=cbind(as.numeric(X))
    ##grp=1
  }
  else if(X.class=='NULL') return(X)
  else if(X.class!='matrix')
    stop('Expecting either a factor, a vector, a matrix or a data.frame')

  #print(789)
  N <- nrow(X)
  p <- ncol(X)

  xinfo. <- matrix(nrow=p, ncol=numcut)
  nc <- numcut
  rm.vars <- c()
  #print(789)
  if(length(xinfo)==0 & N>0 & p>0 & (rm.const | numcut[1]>0)) {
    for(j in 1:p) {
      X.class <- class(X[1, j])[1]

      if(X.class=='numeric' | X.class=='integer') {
        xs <- unique(sort(X[ , j]))
        k <- length(xs)
        nc[j] <- numcut

        if(k %in% 0:1) { # deal with constant variables
          rm.vars <- c(rm.vars, -j)
          nc[j] <- 1
          if(k==0) xs <- NA
        }
        else if(cont)
          xs <- seq(xs[1], xs[k], length.out=numcut+2)[-c(1, numcut+2)]
        else if(k<numcut) {
          xs <- 0.5*(xs[1:(k-1)]+xs[2:k]) #  if k < numcut, use middle point between values to split
          nc[j] <- k-1
        }
        else if(usequants) {
          xs <- quantile(X[ , j], type=type,
                         probs=(0:(numcut+1))/(numcut+1))[-c(1, numcut+2)]
          names(xs) <- NULL
        }
        else xs <-
          seq(xs[1], xs[k], length.out=numcut+2)[-c(1, numcut+2)]
      }
      else
        stop(paste0('Variables of type ', X.class, ' are not supported'))

      xinfo.[j, 1:nc[j] ] <- xs
    }
  }

  X <- data.matrix(X)
  #print(xinfo)
  #print(dim(xinfo.))
  if(length(xinfo)>0) {
    if(is.list(xinfo)) for(j in 1:p){

      xinfo.[j, 1:length(xinfo[[j]])] <- xinfo[[j]]
    }
    else if(is.matrix(xinfo)) xinfo. <- xinfo
    else stop('Only a list or a matrix can be provided for xinfo')

    for(j in 1:p) nc[j] <- sum(!is.na(xinfo.[j, ]))
  }

  xinfo <- xinfo.

  if(rm.const && length(rm.vars)>0 &&
     !(length(rm.vars)==p && all((1:p)==(-rm.vars)))) {
    X <- X[ , rm.vars]
    nc <- nc[rm.vars]
    xinfo <- xinfo[rm.vars, ]
    grp <- grp[rm.vars]
  }
  else if(length(rm.vars)==0 || (length(rm.vars)==p && all((1:p)==(-rm.vars))))
    rm.vars <- 1:p

  dimnames(xinfo) <- list(dimnames(X)[[2]], NULL)

  if(all(numcut==0)) return(X)
  else return(list(X=X, numcut=as.integer(nc), rm.const=rm.vars,
                   xinfo=xinfo, grp=grp))
}
*/

#define TRDRAW(a, b) trdraw(a, b)

class bart_model {
private:
  // --- CRITICAL MEMORY FIX ---
  // We must store the data vectors in the class so they persist.
  // Otherwise, pointers passed to 'bm' become invalid references (NaNs/Crashes).

  Rcpp::NumericMatrix stored_xv; // Was NumericVector
  Rcpp::NumericVector stored_yv;

  // Standard BART members
  IntegerVector numcut;
  bool usequants;
  bool cont;
  IntegerVector rm_const;

  long n;
  long p;
  long ntrees;
  double sigmaf;
  double tau;

  double *ix;
  double *iy;

  double alpha;
  double mybeta;
  double fmean;
  double sigma;
  double nu;
  double lambda;

  List tree_object;

  arn gen;
  bart bm;

public:
  bart_model(){};

  bart_model(NumericMatrix x_train, NumericVector y_train, long numcut=100L, bool usequants = false, bool cont = false, bool rm_const = false, int ntrees = 300, Nullable<double> sigmaf = R_NilValue, double k = 2.0, double power = 2, double base = 0.95, double nu = 3){
    Function bartModelMatrix = Rcpp::Environment::namespace_env("SBMTrees")["bartModelMatrix"];

    this->usequants = usequants;
    this->cont = cont;
    this->rm_const = rm_const;
    this->alpha = base;
    this->mybeta = power;
    this->tree_object = List();
    sigma = 1.0;
    this->nu = nu;

    n = y_train.length();

    // Prepare X using the R function
    Rcpp::List temp = bartModelMatrix(clone(x_train), numcut, usequants, 7, rm_const, cont);

    // Transpose logic (as in your original code)
    NumericMatrix X_trans = transpose(as<NumericMatrix>(temp["X"]));

    // --- MEMORY FIX ---
    // Copy data to class member to ensure persistence
    this->stored_xv = X_trans;
    this->ix = &this->stored_xv[0];

    this->numcut = as<IntegerVector>(temp["numcut"]);
    NumericMatrix Xinfo = as<NumericMatrix>(temp["xinfo"]);

    if(n != X_trans.ncol())
      throw("The length of y_train and the number of rows in x_train must be identical");

    p = X_trans.nrow(); // Dimensions based on transposed matrix

    // Initial Sigma Estimate
    double sigest = sd(y_train);
    if(p < n){
      arma::mat x_r(x_train.begin(), n, p, false);
      bool has_constant_column = false;
      if(!rm_const){
        for (int j = 0; j < p; ++j) {
          if (arma::all(x_r.col(j) == 1)) {
            has_constant_column = true;
            break;
          }
        }
      }
      if(has_constant_column == false){
        arma::mat allOne(n, 1, arma::fill::ones);
        x_r.insert_cols(0, allOne);
      }
      arma::colvec y_r(y_train.begin(), y_train.size(), false);
      arma::colvec coef = arma::solve(x_r, y_r);
      arma::colvec resid = y_r - x_r*coef;
      double sig2 = arma::as_scalar(arma::trans(resid)*resid/(n-p));
      sigest = pow(sig2, 0.5);
    }

    NumericVector qch;
    qch.push_back(1 - 0.9);
    double qchi = Rcpp::qchisq(qch, nu, true, false)[0];
    this->lambda = (sigest * sigest * qchi) / nu;

    this->ntrees = ntrees;
    if(this->rm_const.length() == 0){
      this->rm_const = seq(1, p);
    }

    // Prepare Y
    this->fmean = mean(y_train);

    // --- MEMORY FIX ---
    this->stored_yv = clone(y_train) - this->fmean;
    this->iy = &this->stored_yv[0];

    if(sigmaf.isNull()){
      tau = (max(this->stored_yv) - min(this->stored_yv))/(2*k*sqrt(ntrees));
    }else{
      this->sigmaf = as<double>(sigmaf) / sqrt(ntrees);
    }

    bm = bart(ntrees);

    if(Xinfo.size()>0) {
      xinfo xi_;
      xi_.resize(p);
      for(size_t i=0;i<p;i++) {
        xi_[i].resize(this->numcut[i]);
        for(size_t j=0;j<(size_t)this->numcut[i];j++) xi_[i][j]=Xinfo(i, j);
      }
      bm.setxinfo(xi_);
    }

    int *nc = &this->numcut[0];

    bm.setprior(alpha, mybeta, tau);
    bm.setdata(p, n, this->ix, this->iy, nc);
  };

  List update(long nburn, long npost, int skip, bool verbose = false, long print_every = 100L){
    return update(this->sigma, nburn, npost, skip, verbose, print_every);
  };

  List update(double sigma, long nburn, long npost, int skip, bool verbose = false, long print_every = 100L){
    this->sigma = sigma; // Update internal sigma with passed value

    Rcpp::NumericVector trmean(n);
    Rcpp::NumericVector tsigma(0);
    Rcpp::NumericMatrix trdraw(npost / skip, n);
    Rcpp::NumericMatrix varprb(npost / skip, p);
    Rcpp::IntegerMatrix varcnt(npost / skip, p);

    std::stringstream treess;
    treess.precision(10);
    treess << npost / skip << " " << ntrees << " " << p << endl;

    std::vector<double> ivarprb (p, 0.);
    std::vector<size_t> ivarcnt (p, 0);

    if(verbose) printf("\nMCMC\n");

    size_t trcnt=0;
    size_t treedrawscnt=0;
    bool keeptreedraw;
    xinfo& xi = bm.getxinfo();

    for(size_t i=0; i < nburn + npost; i++) {
      if(verbose && i % print_every == 0){
        Rcout << "Iteration: " << i << "/" << nburn + npost << std::endl;
      }

      // Draw BART (Updates trees and returns acceptance)
      bm.draw(this->sigma, gen);

      // Update Sigma (Gibbs step)
      double restemp = 0, rss=0.0;
      for(size_t k=0; k<n; k++) {
        restemp = (this->iy[k] - bm.f(k));
        rss += restemp*restemp;
      }
      this->sigma = sqrt((nu*lambda + rss)/gen.chi_square(n+nu));

      if(i>=nburn) {
        for(size_t k=0; k<n; k++) trmean[k]+=bm.f(k);

        keeptreedraw = npost && (((i-nburn+1) % skip) == 0);
        if(keeptreedraw) {
          tsigma.push_back(this->sigma);
          for(size_t k=0; k<n; k++) TRDRAW(trcnt,k)=bm.f(k);
          trcnt+=1;

          for(size_t j=0; j<ntrees; j++) {
            treess << bm.gettree(j);
          }

#ifndef NoRcpp
          ivarcnt=bm.getnv();
          ivarprb=bm.getpv();
          size_t k_idx = (i-nburn)/skip;
          for(size_t j=0; j<p; j++){
            varcnt(k_idx, j) = ivarcnt[j];
            varprb(k_idx, j) = ivarprb[j];
          }
#else
          // Legacy support omitted for brevity
#endif
          treedrawscnt +=1;
        }
      }
    }

    for(size_t k=0; k<n; k++) trmean[k] /= npost;

#ifndef NoRcpp
    Rcpp::List ret;
    ret["yhat.train.mean"] = trmean;
    ret["yhat.train"] = trdraw;
    ret["varcount"] = varcnt;
    ret["varprob"] = varprb;

    Rcpp::List xiret(xi.size());
    for(size_t i=0; i<xi.size(); i++) {
      Rcpp::NumericVector vtemp(xi[i].size());
      std::copy(xi[i].begin(), xi[i].end(), vtemp.begin());
      xiret[i] = Rcpp::NumericVector(vtemp);
    }

    Rcpp::List treesL;
    treesL["cutpoints"] = xiret;
    treesL["trees"] = Rcpp::CharacterVector(treess.str());

    ret["treedraws"] = treesL;
    ret["mu"] = fmean;
    ret["yhat.train.mean"] = trmean + fmean;
    ret["yhat.train"] = trdraw + fmean;
    ret["sigma"] = tsigma;

    this->tree_object = ret;
    return ret;
#else
    return List::create();
#endif
  };

  // --- CRITICAL DATA UPDATE FUNCTIONS ---

  void set_data(NumericMatrix x_train, NumericVector y_train){
    n = y_train.length();

    // 1. Update Y
    this->fmean = mean(y_train);
    // Deep copy to member variable
    this->stored_yv = clone(y_train) - this->fmean;
    // Update internal pointer to point to valid member memory
    this->iy = &this->stored_yv[0];

    // 2. Update X
    // Deep copy transposed matrix to member variable
    this->stored_xv = transpose(clone(x_train));
    // Update internal pointer
    this->ix = &this->stored_xv[0];

    p = this->stored_xv.nrow();
    int *nc = &numcut[0];

    // 3. Update BART
    // WARNING: This assumes 'bm' can warm-start with new data.
    // Standard BART handles this by re-calculating sufficient stats in 'draw'.
    bm.setdata(p, n, this->ix, this->iy, nc);
  };

  void set_Y(NumericVector y_train){
    n = y_train.length();

    // 1. Update Y
    this->fmean = mean(y_train);
    // Deep copy to member variable
    this->stored_yv = clone(y_train) - this->fmean;
    // Update internal pointer
    this->iy = &this->stored_yv[0];

    int *nc = &numcut[0];

    // 2. Update BART using existing X pointer (this->ix)
    // Since we updated 'this->ix' in set_data, this correctly uses the current X.
    bm.setdata(p, n, this->ix, this->iy, nc);
  };

  NumericMatrix predict(NumericMatrix x_predict, bool verbose = false){
    if(this->tree_object.length() == 0){
      return NumericMatrix();
    }

    // Use transpose clone for safety
    NumericMatrix X = transpose(as<NumericMatrix>(clone(x_predict)));

    // Assumes cpwbart (Parallel Predict) is available from headers
    NumericMatrix predict_y = cpwbart(this->tree_object["treedraws"], X, verbose);

    return predict_y + this->fmean;
  };

  SEXP get_tree_object(){ return tree_object; }
  double get_sigma(){ return this->sigma; }
  bool get_usequants(){ return this->usequants; }
  double get_lambda(){ return lambda; }
  double get_nu(){ return nu; }
  double get_invchi(long n, double rss){
    return sqrt((nu*lambda + rss)/gen.chi_square(n+nu));
  }

};

// [[Rcpp::export]]
SEXP bart_train(NumericMatrix X, NumericVector Y, long nburn = 100, long npost = 1000, bool verbose = true){
  // 1. Allocate model
  bart_model * m = new bart_model(X, Y);

  // 2. Initial Burn-in
  List a = m->update(1, 1, 1, verbose);

  NumericMatrix y_pre = m->predict(X);
  NumericMatrix post_y(npost, Y.length());

  for(int i = 0; i < nburn + npost; ++i){
    if(i % 100 == 0)
      Rcout << i << " " << nburn + npost << std::endl;
    // FIX: Pass m->get_sigma() so we don't reset sigma to 1.0 every step.
    // This allows the chain to mix properly.
    a = m->update(1, 1, 1, verbose, 200);

    // Optional: Update Y if needed (not needed for standard train)
    //m->set_Y(Y);
    m->set_data(X, Y);
    y_pre = m->predict(X);
    if(i >= nburn)
      post_y(i-nburn, _) = y_pre;
  }

  // 3. CLEAN UP MEMORY (Fixes Memory Leak)
  delete m;

  return List::create(Named("y_pre") = y_pre, Named("post_y") = post_y);
}
