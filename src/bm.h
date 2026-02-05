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




#ifndef UPDATE_H_
#define UPDATE_H_
#include "DP_lambda.h"
#include "cal_random_effects.h"
#include "update_B.h"
#include "update_Covariance.h"
#include <cmath>
#endif

using namespace Rcpp;

class bm{
public:
  virtual ~bm() = default;

  bm(NumericVector Y, NumericMatrix X, NumericMatrix Z, CharacterVector subject_id, IntegerVector row_id, bool binary = false, bool CDP_residual = false, bool CDP_re = false, double tol=1e-40, int ntrees = 200, double k = 2.0, double pi_CDP = 0.99, bool train = true) {     // Constructor
    Rcout << 2123 << std::endl;
    if(train){
      this->tol = tol;
      this->CDP_residual = CDP_residual;
      this->CDP_re = CDP_re;

      d = as<NumericMatrix>(Z).ncol();

      this->binary = binary;
      this->Y_original = clone(Y);
      this->Y = clone(Y);
      this->X = clone(X);
      this->z = clone(Z);
      this->subject_id = subject_id;
      this->row_id = row_id;

      N = Y.length();
      p = X.ncol();
      n_subject = unique(subject_id).length();

      subject_to_B = create_subject_to_B(subject_id);
      row_id_to_id = create_row_id_to_row(row_id);

      if(binary){
        NumericVector y_ = clone(this->Y);
        y_ = y_ * 2 - 1;
        this->Y = clone(y_);
      }

      tau_samples = NumericVector(N);
      B_tau_samples = NumericMatrix(n_subject, d);
      Function get_inverse_wishart_matrix2 = Environment::namespace_env("SBMTrees")["get_inverse_wishart_matrix2"];
      //Environment pkg_env = Environment::global_env();
      //Function get_inverse_wishart_matrix2 = pkg_env["get_inverse_wishart_matrix2"];
      List lmm = get_inverse_wishart_matrix2(this->X, this->Y, z, subject_id, subject_to_B);

      sigma = as<double>(lmm["sigma"]);
      NumericMatrix coe = as<NumericMatrix>(lmm["coe"]);
      inverse_wishart_matrix = as<NumericMatrix>(lmm["covariance"]);
      if(CDP_re){
        M_re = pow(n_subject, (double)(runif(1, 0, 0.5)[0]));
        B_tau = DP(List::create(Named("p") = d, Named("cov") = as<NumericMatrix>(inverse_wishart_matrix)), M_re, sqrt(n_subject), n_subject, true);
        B_tau_samples = as<NumericMatrix>(B_tau["samples"]);
        Covariance = as<NumericMatrix>(B_tau["Sigma"]);
      }
      if(CDP_residual){
        M = pow(N, (double)(runif(1, 0, 0.5)[0]));
        //Rcout << M << std::endl;
        //Rcout <<  (double)lmm["sigma"] << std::endl;
        if(binary){
          if(CDP_re)
            tau = DP(List::create(Named("p") = 1, Named("sd") = 1, Named("pi") = pi_CDP), M, sqrt(N), N, true);
          else
            tau = DP(List::create(Named("p") = 1, Named("sd") = 1, Named("pi") = pi_CDP), M, sqrt(N), N, true);
        }else{
          if(CDP_re)
            tau = DP(List::create(Named("p") = 1, Named("sd") = (double)lmm["sigma"], Named("pi") = pi_CDP), M, sqrt(N), N, true);
          else
            tau = DP(List::create(Named("p") = 1, Named("sd") = (double)lmm["sigma"], Named("pi") = pi_CDP), M, sqrt(N), N, true);
        }
        tau_samples = NumericVector(as<NumericMatrix>(tau["samples"])(_,0));
        //Rcout <<tau_samples << std::endl;
      }
      if(!CDP_re)
        Covariance = inverse_wishart_matrix;
      B = NumericMatrix(n_subject, d);
      //B = coe;
      cov = as<arma::mat>(Covariance);
      re = cal_random_effects(z, subject_id, B, subject_to_B, Covariance, binary);
      row_id_to_re = create_row_id_to_re(row_id, re);
      new_subject_B.clear();
    }
  }


  NumericVector get_Y(){
    return this->Y;
  }

  NumericVector get_re(){
    return this->re;
  }

  double get_tau_sigma(){
    return sigma;
  }

  long get_R_cluster_number(){
    if(CDP_residual){
      IntegerVector CDP_residual_cluster = as<IntegerVector>(tau["cluster"]);
      //Rcout << CDP_residual_cluster << std::endl;
      IntegerVector nk = table(CDP_residual_cluster);
      long cnt = std::count_if(nk.begin(), nk.end(),
                              [](int v){ return v > 0; });
      //Rcout << nk << std::endl;
      //Rcout << cnt << std::endl;
      return cnt;
    }else{
      return 0;
    }
  }

  long get_RE_cluster_number(){
    if(CDP_re){
      IntegerVector CDP_re_cluster = as<IntegerVector>(B_tau["cluster"]);
      IntegerVector nk = table(CDP_re_cluster);
      long cnt = std::count_if(nk.begin(), nk.end(),
                               [](int v){ return v > 0; });
      //Rcout << cnt << std::endl;
      return cnt;
    }else{
      return 0;
    }
  }

  void update_X_Y(NumericMatrix X, NumericVector Y){
    this->X = clone(X);
    this->Y_original = clone(Y);
    this->Y = clone(Y);
    if(binary){
      NumericVector y_ = clone(this->Y);
      y_ = y_ * 2 - 1;
      this->Y = clone(y_);
    }
  }


  virtual void update_sigma() = 0;
  virtual void update_all(bool verbose = false) = 0;
  virtual List posterior_sampling() = 0;
  virtual NumericVector predict_fix(NumericMatrix X_test) = 0;


  List get_tree_training_data(){
    return List::create(Named("X") = X, Named("Y") = Y - re - tau_samples);
  }


  List get_CDP_re_data(bool verbose = false){
    if(CDP_re){
      if(verbose)
        Rcout << "update DP" << std::endl;
      B_tau["Sigma"] = Covariance;
      if(d == 1){
        B_tau["sigma"] = Covariance;
      }
      return(List::create(Named("B") = B, Named("B_tau") = B_tau));
    }else{
      return(List::create());
    }
  }

  void set_CDP_re_data(List B_tau){
    this->B_tau = B_tau;
    B_tau_samples = as<NumericMatrix>(B_tau["samples"]);
    M_re = B_tau["M"];
  }


  NumericVector get_tau_samples(){
    return tau_samples;
  }

  NumericVector get_tau_mu(){
    return tau["y"];
  }

  NumericVector get_tau_pi(){
    return tau["pi"];
  }

  NumericVector get_B_tau_samples(){
    return B_tau_samples;
  }

  NumericVector get_B_tau_mu(){
    return B_tau["y"];
  }

  NumericVector get_B_tau_pi(){
    return B_tau["pi"];
  }

  NumericVector get_B_tau_lambda(){
    return B_tau["lambda"];
  }


  NumericVector predict_expectation(NumericMatrix X_test, NumericMatrix Z_test, CharacterVector subject_id_test, IntegerVector row_id_test, bool keep_re = true){
    int n = X_test.nrow();

    NumericVector y_pre = predict_fix(X_test);
    for(int i = 0 ; i < n ; ++i){
      y_pre[i] = y_pre[i] + predict_re(Z_test(i,_), std::string(subject_id_test[i]), row_id_test[i]);
    }
    if(this->binary == true){
      if(CDP_residual){
        return probit_probability_tau(tau, y_pre);
      }else{
        return Rcpp::pnorm(y_pre);
      }
    }else{
      return y_pre;
    }
  }

  double predict_re(const NumericVector& Z_test_i,
                    const std::string& subject_id_i,
                    int row_id_i) {

    const std::string rid = std::to_string(row_id_i);

    // 1) row-level cache (training or test)
    if (auto it = row_id_to_re.find(rid); it != row_id_to_re.end())
      return it->second;

    NumericVector Bi;

    // 2) in-sample subject => use fitted B
    if (auto it = subject_to_B.find(subject_id_i); it != subject_to_B.end()) {
      Bi = B(it->second, _);

      // 3) out-of-sample subject but already drawn this iteration
    } else if (auto it = new_subject_B.find(subject_id_i); it != new_subject_B.end()) {
      Bi = it->second;

      // 4) out-of-sample subject => draw once, cache
    } else {
      if (CDP_re) {
        NumericVector mu = predict_tau(B_tau, 1)(0, _);
        arma::mat draw = rmvnorm(1, as<arma::vec>(mu), cov);
        Bi = wrap(draw.row(0));
      } else {
        arma::vec mu0(d, arma::fill::zeros);
        arma::mat draw = rmvnorm(1, mu0, cov);
        Bi = wrap(draw.row(0));
      }
      new_subject_B.emplace(subject_id_i, Bi);
    }

    const double re_i = innerProduct(Z_test_i, Bi);
    row_id_to_re.emplace(rid, re_i);
    return re_i;
  }



  double predict_eta(int row_id_i) {
    if (CDP_residual){
      const std::string rid = std::to_string(row_id_i);
      // 1) row-level cache (training or test)
      if (auto it = row_id_to_id.find(rid); it !=row_id_to_id.end())
        return tau_samples[it->second];
      return predict_tau(tau, 1)(0, 0);
    }else{
      return 0;
    }
  }




  NumericVector predict_sample(NumericMatrix X_test, NumericMatrix Z_test, CharacterVector subject_id_test, IntegerVector row_id_test, bool keep_re = true){
    int n = X_test.nrow();

    NumericVector y_pre = predict_fix(X_test);
    for(int i = 0 ; i < n ; ++i){
      y_pre[i] = y_pre[i] + predict_re(Z_test(i,_), std::string(subject_id_test[i]), row_id_test[i]);
      y_pre[i] = y_pre[i] + predict_eta(row_id_test[i]);
      y_pre[i] = R::rnorm(y_pre[i], sigma);
      if(this->binary == true){
        y_pre[i] = y_pre[i] >= 0;
      }
    }
    return y_pre;
  }

  double predict_probability_log(double Y_test,
                                 double fX_test,
                                 const NumericVector& Zi,
                                 const std::string& subject_id_test,
                                 int row_id_test) {

    auto clamp01 = [](double p){
      const double eps = 1e-15;
      if(p < eps) return eps;
      if(p > 1.0 - eps) return 1.0 - eps;
      return p;
    };

    auto logsumexp = [](const std::vector<double>& lt){
      double m = -INFINITY;
      for(double v : lt) m = std::max(m, v);
      if(!R_finite(m)) return m;
      double s = 0.0;
      for(double v : lt) s += std::exp(v - m);
      return m + std::log(s);
    };

    const std::string rid = std::to_string(row_id_test);

    // -------- Case 1: in-sample row y_ij (condition on stored b_i and eta_ij) --------
    if (auto it_id = row_id_to_id.find(rid); it_id != row_id_to_id.end()) {
      const int idx = it_id->second;

      double re_i = 0.0;
      if (auto it_re = row_id_to_re.find(rid); it_re != row_id_to_re.end()) re_i = it_re->second;

      const double eta_i = (CDP_residual ? tau_samples[idx] : 0.0);
      const double mu = fX_test + re_i + eta_i;

      if (!binary) {
        return R::dnorm(Y_test, mu, sigma, true);
      } else {
        if (Y_test != 0.0 && Y_test != 1.0) Rcpp::stop("Y_test must be 0 or 1 for binary probit.");
        const double p = clamp01(R::pnorm(mu, 0.0, 1.0, true, false));
        return (Y_test == 1.0) ? std::log(p) : std::log1p(-p);
      }
    }

    // Helpers for mixture atoms/weights:
    NumericMatrix eta_atom(1,1); eta_atom(0,0) = 0.0;
    NumericVector eta_w(1); eta_w[0] = 1.0;
    if (CDP_residual) { eta_atom = as<NumericMatrix>(tau["y"]); eta_w = as<NumericVector>(tau["pi"]); }

    // -------- Case 2: out-of-sample row but in-sample subject y_ij' (integrate out eta only) --------
    if (auto itB = subject_to_B.find(subject_id_test); itB != subject_to_B.end()) {
      const double re_i = innerProduct(Zi, B(itB->second, _));
      const double base = fX_test + re_i;

      if (!binary) {
        std::vector<double> lt(eta_atom.nrow());
        for (int h = 0; h < eta_atom.nrow(); ++h) {
          const double w = eta_w[h];
          lt[h] = (w <= 0.0) ? -INFINITY
          : std::log(w) + R::dnorm(Y_test, base + eta_atom(h,0), sigma, true);
        }
        return logsumexp(lt);
      } else {
        if (Y_test != 0.0 && Y_test != 1.0) Rcpp::stop("Y_test must be 0 or 1 for binary probit.");
        double p1 = 0.0;
        for (int h = 0; h < eta_atom.nrow(); ++h) {
          const double w = eta_w[h];
          if (w > 0.0) p1 += w * R::pnorm(base + eta_atom(h,0), 0.0, 1.0, true, false);
        }
        p1 = clamp01(p1);
        return (Y_test == 1.0) ? std::log(p1) : std::log1p(-p1);
      }
    }

    // -------- Case 3: out-of-sample subject y_i'j' (integrate out eta and b) --------
    arma::vec z = as<arma::vec>(Zi);
    const double zSz = arma::as_scalar(z.t() * cov * z);

    NumericMatrix mu_atom(1, d); std::fill(mu_atom.begin(), mu_atom.end(), 0.0);
    NumericVector mu_w(1); mu_w[0] = 1.0;
    if (CDP_re) { mu_atom = as<NumericMatrix>(B_tau["y"]); mu_w = as<NumericVector>(B_tau["pi"]); }

    if (!binary) {
      const double sd = std::sqrt(zSz + sigma * sigma);

      std::vector<double> lt(eta_atom.nrow() * mu_atom.nrow());
      int t = 0;
      for (int h = 0; h < eta_atom.nrow(); ++h) {
        const double we = eta_w[h];
        if (we <= 0.0) { t += mu_atom.nrow(); continue; }
        const double eta_h = eta_atom(h,0);

        for (int hb = 0; hb < mu_atom.nrow(); ++hb) {
          const double wb = mu_w[hb];
          if (wb <= 0.0) { lt[t++] = -INFINITY; continue; }

          const double zmu = innerProduct(Zi, mu_atom(hb, _));
          const double mean = fX_test + zmu + eta_h;

          lt[t++] = std::log(we) + std::log(wb) + R::dnorm(Y_test, mean, sd, true);
        }
      }
      return logsumexp(lt);

    } else {
      if (Y_test != 0.0 && Y_test != 1.0) Rcpp::stop("Y_test must be 0 or 1 for binary probit.");

      const double denom = std::sqrt(zSz + 1.0); // probit: latent variance fixed to 1
      double p1 = 0.0;

      for (int h = 0; h < eta_atom.nrow(); ++h) {
        const double we = eta_w[h];
        if (we <= 0.0) continue;
        const double eta_h = eta_atom(h,0);

        for (int hb = 0; hb < mu_atom.nrow(); ++hb) {
          const double wb = mu_w[hb];
          if (wb <= 0.0) continue;

          const double zmu = innerProduct(Zi, mu_atom(hb, _));
          const double arg = (fX_test + zmu + eta_h) / denom;

          p1 += we * wb * R::pnorm(arg, 0.0, 1.0, true, false);
        }
      }

      p1 = clamp01(p1);
      return (Y_test == 1.0) ? std::log(p1) : std::log1p(-p1);
    }
  }

  LogicalVector replace_indicator(double correct_prob){
    LogicalVector replace(X.nrow());
    if(CDP_residual){
      NumericVector CDP_residual_pi = tau["pi"];
      IntegerVector CDP_residual_cluster = tau["cluster"];
      for(int i = 0 ; i < X.nrow(); ++i){
        int k = CDP_residual_cluster[i];
        std::string n_str = std::to_string(k);
        const char *n = n_str.c_str();
        if(CDP_residual_pi.containsElementNamed(n)){
          if(CDP_residual_pi[n] < correct_prob){
            replace[i] = true;
          }
        }
      }
    }
    if(CDP_re){
      NumericVector CDP_re_pi = B_tau["pi"];
      IntegerVector CDP_re_cluster = B_tau["cluster"];
      //Rcout << CDP_re_cluster.length() << std::endl;
      for(int i = 0 ; i < X.nrow(); ++i){
        String name = subject_id[i];
        auto it = subject_to_B.find(std::string(name));
        if (it != subject_to_B.end()){
          long B_position = (int)(it->second);
          int k = CDP_re_cluster[B_position];
          std::string n_str = std::to_string(k);
          const char *n = n_str.c_str();
          if(CDP_re_pi.containsElementNamed(n)){
            if(CDP_re_pi[n] < correct_prob){
              replace[i] = true;
            }
          }
        }
      }
    }
    return replace;
  }


protected:
  double tol;

  Environment G;

  int d;
  bool binary;
  NumericVector Y_original;
  NumericVector Y;
  NumericMatrix X;
  NumericMatrix z;
  CharacterVector subject_id;
  IntegerVector row_id;

  long N;
  int p;
  int n_subject;

  std::unordered_map<std::string, int> subject_to_B;
  std::unordered_map<std::string, int> row_id_to_id;
  std::unordered_map<std::string, double> row_id_to_re;
  std::unordered_map<std::string, NumericVector> new_subject_B;

  NumericMatrix inverse_wishart_matrix;
  NumericMatrix Covariance;
  NumericMatrix B;

  double M_re = 0;
  double M = 0;
  double sigma = 1;
  List tau;
  List tau_prop;

  List B_tau;
  NumericVector tau_samples;
  NumericMatrix B_tau_samples;

  bool CDP_residual;
  bool CDP_re;
  NumericVector re;
  arma::vec re_arma;

  NumericVector random_test;
  NumericVector re_test;

  arma::mat cov;

};
