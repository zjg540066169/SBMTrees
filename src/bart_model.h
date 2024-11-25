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
    print("xinfo")
    print(xinfo)
    if(is.list(xinfo)) for(j in 1:p){
      
      xinfo.[j, 1:length(xinfo[[j]])] <- xinfo[[j]]
      
    } 
    else if(is.matrix(xinfo)) xinfo. <- xinfo
    else stop('Only a list or a matrix can be provided for xinfo')
    print("xinfo.")
    print(xinfo.)
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

class bart_model{
public:
  bart_model(){};
  
  bart_model(NumericMatrix x_train, NumericVector y_train, long numcut=100L, bool usequants = false, bool cont = false, bool rm_const = false, int ntrees = 300, Nullable<double> sigmaf = R_NilValue, double k = 2.0, double power = 2, double base = 0.95, double nu = 3){
    //Rcout << 123 << std::endl;
    G = Rcpp::Environment::global_env();
    //Rcout << 123 << std::endl;
    Function bartModelMatrix = G["bartModelMatrix"];
    //Rcout << 123 << std::endl;
    this->usequants = usequants;
    this->cont = cont;
    this->rm_const = rm_const;
    this->alpha = base;
    this->mybeta = power;
    this->tree_object = List();
    sigma = 1;
    this->nu = nu;
    
    n = y_train.length();
    Rcpp::List temp = bartModelMatrix(clone(x_train), numcut, usequants, 7, rm_const, cont);
    NumericMatrix X = transpose(as<NumericMatrix>(temp["X"]));
    
    this->numcut = as<IntegerVector>(temp["numcut"]);
    NumericMatrix Xinfo = as<NumericMatrix>(temp["xinfo"]);
    
    //Rcout << 123 << std::endl;
    //Rcout << "bartModelMatrix" << std::endl;
    if(n!=X.ncol())
      throw("The length of y_train and the number of rows in x_train must be identical");
    
    p = X.nrow();
    
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
      arma::colvec coef = arma::solve(x_r, y_r);      // fit model y ~ X
      arma::colvec resid = y_r - x_r*coef;            // residuals
      
      double sig2 = arma::as_scalar( arma::trans(resid)*resid/(n-p));
      sigest = pow(sig2, 0.5);
      //Rcout << sigest << std::endl;
    }
    NumericVector qch;
    qch.push_back(1 - 0.9);
    double qchi = Rcpp::qchisq(qch, nu, true, false)[0];
    this->lambda = (sigest * sigest * qchi) / nu;
    
    
    this->ntrees = ntrees;
    if(this->rm_const.length() == 0){
      this->rm_const = seq(1, p);
    }
    
    this->fmean = mean(y_train);
    NumericVector y = clone(y_train)-this->fmean;
    //Rcout << y.length() << std::endl;
    if(sigmaf.isNull()){
      tau = (max(y)-min(y))/(2*k*sqrt(ntrees));
    }else{
      this->sigmaf = as<double>(sigmaf) / sqrt(ntrees);
    }
    
    bm = bart(ntrees);
    
    if(Xinfo.size()>0) {
      xinfo xi_;
      xi_.resize(p);
      for(size_t i=0;i<p;i++) {
        xi_[i].resize(this->numcut[i]);
        //Rcpp::IntegerVector cutpts(Xinfo[i]);
        for(size_t j=0;j<(size_t)this->numcut[i];j++) xi_[i][j]=Xinfo(i, j);
      }
      bm.setxinfo(xi_);
    }
    
    Rcpp::NumericVector  xv(X);
    ix = &xv[0];
    Rcpp::NumericVector  yv(y); 
    iy = &yv[0];
    int *nc = &this->numcut[0];
    
    
    //heterbart bm(m);
    bm.setprior(alpha,mybeta,tau);
    bm.setdata(p,n,ix,iy, nc);
    //Rcout << "finish initialization" << std::endl;
  };
 
  List update(long nburn, long npost, int skip, bool verbose = false, long print_every = 100L){
    Rcpp::NumericVector trmean(n); //train
    Rcpp::NumericVector tsigma(0); //train
    Rcpp::NumericMatrix trdraw(npost / skip,n);
    Rcpp::NumericMatrix varprb(npost / skip,p);
    Rcpp::IntegerMatrix varcnt(npost / skip,p);
    
    std::stringstream treess;  //string stream to write trees to  
    treess.precision(10);
    treess << npost / skip << " " << ntrees << " " << p << endl;
    std::vector<double> ivarprb (p,0.);
    std::vector<size_t> ivarcnt (p,0);
    
    
    if(verbose)
      printf("\nMCMC\n");
    size_t trcnt=0;
    size_t treedrawscnt=0; //count kept bart draws
    bool keeptreedraw;
    xinfo& xi = bm.getxinfo();
    //Rcout << "sigma:" << this->sigma<<std::endl;
    for(int i=0; i < nburn + npost;i++) {
      if(verbose){
        if(i % print_every == 0){
          printf("iteration %d",i);
          Rcout << "/"<<nburn + npost <<std::endl;
        }
        
        
      }
      //draw bart
      //Rcout << "before draw" << std::endl;
      bm.draw(sigma,gen);
      //Rcout << "after draw" << std::endl;
      double restemp = 0, rss=0.0;
      for(size_t k=0;k<n;k++) {restemp=(iy[k]-bm.f(k)); rss += restemp*restemp;}
      //Rcout << "rss:" << rss << std::endl;
      sigma = sqrt((nu*lambda + rss)/gen.chi_square(n+nu));
      
      if(i>=nburn) {
        for(size_t k=0;k<n;k++) trmean[k]+=bm.f(k);
        keeptreedraw = npost && (((i-nburn+1) % skip) ==0);
        if(keeptreedraw) {
          //Rcout << sigma << std::endl;
          tsigma.push_back(sigma);
          for(long k=0;k<n;k++) TRDRAW(trcnt,k)=bm.f(k);
          trcnt+=1;
          for(size_t j=0;j<ntrees;j++) {
            treess << bm.gettree(j);
          }
#ifndef NoRcpp
          ivarcnt=bm.getnv();
          ivarprb=bm.getpv();
          size_t k=(i-nburn)/skip;
          for(size_t j=0;j<p;j++){
            varcnt(k,j)=ivarcnt[j];
            //varcnt(i-burn,j)=ivarcnt[j];
            varprb(k,j)=ivarprb[j];
            //varprb(i-burn,j)=ivarprb[j];
          }
#else
          varcnt.push_back(bm.getnv());
          varprb.push_back(bm.getpv());
#endif
          
          treedrawscnt +=1;
        }
      }
    }
    if(verbose){
      printf("iteration %zu complete\n", nburn + npost);
    }
    for(size_t k=0;k<n;k++) trmean[k]/=npost;
    
#ifndef NoRcpp
    Rcpp::List ret;
    //ret["sigma"]=sdraw;
    //Rcout << trmean << std::endl;
    ret["yhat.train.mean"]=trmean;
    ret["yhat.train"]=trdraw;
    //ret["varcount"]=varcount;
    ret["varcount"]=varcnt;
    ret["varprob"]=varprb;
    Rcpp::List xiret(xi.size());
    for(size_t i=0;i<xi.size();i++) {
      Rcpp::NumericVector vtemp(xi[i].size());
      std::copy(xi[i].begin(),xi[i].end(),vtemp.begin());
      xiret[i] = Rcpp::NumericVector(vtemp);
    }
    
    Rcpp::List treesL;
    treesL["cutpoints"] = xiret;
    treesL["trees"]=Rcpp::CharacterVector(treess.str());
    //   if(treesaslists) treesL["lists"]=list_of_lists;
    ret["treedraws"] = treesL;
    ret["mu"] = fmean;
    ret["yhat.train.mean"] = trmean+fmean;
    ret["yhat.train"] = trdraw+fmean;
    ret["sigma"] = tsigma;
    this->tree_object = ret;
    return ret;
#else
    
#endif
    
    
    //return List::create();
  };
  
  List update(double sigma, long nburn, long npost, int skip, bool verbose = false, long print_every = 100L){
    this->sigma = sigma;
    Rcpp::NumericVector trmean(n); //train
    Rcpp::NumericMatrix trdraw(npost / skip,n);
    Rcpp::NumericMatrix varprb(npost / skip,p);
    Rcpp::IntegerMatrix varcnt(npost / skip,p);
    
    std::stringstream treess;  //string stream to write trees to  
    treess.precision(10);
    treess << npost / skip << " " << ntrees << " " << p << endl;
    std::vector<double> ivarprb (p,0.);
    std::vector<size_t> ivarcnt (p,0);
    
    
    if(verbose)
      printf("\nMCMC\n");
    size_t trcnt=0;
    size_t treedrawscnt=0; //count kept bart draws
    bool keeptreedraw;
    xinfo& xi = bm.getxinfo();
    //Rcout << "sigma:" << this->sigma<<std::endl;
    for(size_t i=0;i< nburn + npost;i++) {
      if(verbose){
        if(i % print_every == 0){
          printf("iteration %zu",i);
          Rcout << "/"<<nburn + npost <<std::endl;
        }
        
        
      }
      //draw bart
      //Rcout << "before draw" << std::endl;
      bm.draw(sigma,gen);
      //Rcout << "after draw" << std::endl;
      
      if(i>=nburn) {
        for(size_t k=0;k<n;k++) trmean[k]+=bm.f(k);
        keeptreedraw = npost && (((i-nburn+1) % skip) ==0);
        if(keeptreedraw) {
          //Rcout << sigma << std::endl;
          for(size_t k=0;k<n;k++) TRDRAW(trcnt,k)=bm.f(k);
          trcnt+=1;
          for(size_t j=0;j<ntrees;j++) {
            treess << bm.gettree(j);
          }
#ifndef NoRcpp
          ivarcnt=bm.getnv();
          ivarprb=bm.getpv();
          size_t k=(i-nburn)/skip;
          for(size_t j=0;j<p;j++){
            varcnt(k,j)=ivarcnt[j];
            //varcnt(i-burn,j)=ivarcnt[j];
            varprb(k,j)=ivarprb[j];
            //varprb(i-burn,j)=ivarprb[j];
          }
#else
          varcnt.push_back(bm.getnv());
          varprb.push_back(bm.getpv());
#endif
          
          treedrawscnt +=1;
        }
      }
    }
    if(verbose){
      printf("iteration %zu complete\n", nburn + npost);
    }
    for(size_t k=0;k<n;k++) trmean[k]/=npost;
    
#ifndef NoRcpp
    Rcpp::List ret;
    //ret["sigma"]=sdraw;
    //Rcout << trmean << std::endl;
    ret["yhat.train.mean"]=trmean;
    ret["yhat.train"]=trdraw;
    //ret["varcount"]=varcount;
    ret["varcount"]=varcnt;
    ret["varprob"]=varprb;
    Rcpp::List xiret(xi.size());
    for(size_t i=0;i<xi.size();i++) {
      Rcpp::NumericVector vtemp(xi[i].size());
      std::copy(xi[i].begin(),xi[i].end(),vtemp.begin());
      xiret[i] = Rcpp::NumericVector(vtemp);
    }
    
    Rcpp::List treesL;
    treesL["cutpoints"] = xiret;
    treesL["trees"]=Rcpp::CharacterVector(treess.str());
    //   if(treesaslists) treesL["lists"]=list_of_lists;
    ret["treedraws"] = treesL;
    ret["mu"] = fmean;
    ret["yhat.train.mean"] = trmean+fmean;
    ret["yhat.train"] = trdraw+fmean;
    ret["sigma"] = sigma;
    this->tree_object = ret;
    return ret;
#else
    
#endif
    
    
    //return List::create();
  };
  
  
  void set_data(NumericMatrix x_train, NumericVector y_train){
    n = y_train.length();
    //Function bartModelMatrix = G["bartModelMatrix"];
    this->fmean = mean(y_train);
    NumericVector y = clone(y_train)-this->fmean;
    //xinfo & xi = bm.getxinfo();
    //long ncu = max(numcut);
    //Rcpp::List temp = bartModelMatrix(clone(x_train), ncu, usequants, R_NilValue, rm_const, cont, xi);
    
    //setxinfo 
    NumericMatrix X = transpose(clone(x_train));
    //numcut = as<IntegerVector>(temp["numcut"]);
    // //NumericMatrix X = transpose(x_train);
    p = X.nrow();
    Rcpp::NumericVector  xv(X);
    // free(ix);
    // free(iy);
    double * ix_new = &xv[0];
    Rcpp::NumericVector  yv(y);
    double * iy_new = &yv[0];
    int *nc = &numcut[0];
    bm.setdata(p,n,ix_new, iy_new, nc);
    
    
    
    //update(sigma, 1, 0, 1, false);
    //delete [] ix;
    //delete [] iy;
    //ix = ix_new;
    //iy = iy_new;
    //bm.setdata(p,n,ix,iy, nc);
  };
  
  NumericMatrix predict(NumericMatrix x_predict, bool verbose = false){
    //Function bartModelMatrix = G["bartModelMatrix"];
    if(this->tree_object.length() == 0){
      return NumericMatrix();
    }
    //xinfo & xi = bm.getxinfo();
    //long ncu = max(numcut);
    //Rcpp::List temp = bartModelMatrix(clone(x_predict), ncu, usequants, R_NilValue, rm_const, cont, xi);
    NumericMatrix X = transpose(as<NumericMatrix>(clone(x_predict)));
    //Rcout << xi << std::endl;
    //return xi;
    //return this->tree_object["treedraws"];
    NumericMatrix predict_y = cpwbart(this->tree_object["treedraws"], X, verbose);
    //Rcout << predict_y << std::endl;
    //Rcout << "predict_Y" << std::endl;
    return predict_y + this->fmean;
  };
  
  SEXP get_tree_object(){
    return tree_object;
  }
  
  double get_sigma(){
    return this->sigma;
  }
  
  bool get_usequants(){
    return this->usequants;
  }
  
  double get_invchi(long n, double rss){
    return sqrt((nu*lambda + rss)/gen.chi_square(n+nu));
  }
  
  double get_lambda(){
    return lambda;
  }
  
  double get_nu(){
    return nu;
  }
  
  
private:
  Environment G;
  Environment base;
  
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
};

// [[Rcpp::export]]
SEXP bart_train(NumericMatrix X, NumericVector Y, long nburn = 100, long npost = 1000, bool verbose = true){
  bart_model * m;
  m = new bart_model(X, Y);
  List a;
  for(int i = 0 ; i < nburn + npost; ++i){
    a = m -> update(1, 1, 1, verbose);
    m->set_data(X, Y);
    NumericMatrix y_pre = m->predict(X);
  }
  
  
  //Rcout<< "after"<<std::endl;
  //system("pause");
  
  //a = m->update(nburn, npost, 3, verbose);
  //List b = m.update(nburn, npost, 3, verbose);
  //return m -> get_tree_object();
  return List::create(Named("p") = m ->predict(X, true), Named("y") = a);
  //return List();
  //return m.update(sigma, 100, 1, 1, true);
}


/*
RCPP_MODULE(bart_model_module){
  using namespace Rcpp;
  class_<bart_model>("bart_model")
    .constructor<NumericMatrix, NumericVector, long, bool, int, Nullable<double>>()
    .method("update", &bart_model::update)
    .method("set_data", &bart_model::set_data)
    .method("predict", &bart_model::predict);
}*/
