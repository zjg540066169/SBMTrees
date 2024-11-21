#include <Rcpp.h>
#ifndef UTILS_H_
#define UTILS_H_
#include "utils.h"
#endif

#ifndef DP_SAMPLER_H_
#define DP_SAMPLER_H_
#include "DP_sampler.h"
#endif

using namespace Rcpp;

/*** R
library(mvtnorm)
library(MCMCpack)
library(cascsim)
library(hbmem)


matrix_multiply = function(a, b){
  return(a%*%b)
}

matrix_add = function(a, b){
  return(a + b)
}

make_symmetric = function(s){
  s[lower.tri(s)] = t(s)[lower.tri(s)]
  return(s)
}
*/

// K: number of DP
// L: number of sample in each DP
// [[Rcpp::export]]
List nDP(List parameters, CharacterVector cluster_id, double M_alpha = 2, double M_beta = 2, long L = 55, long K = 35, bool CDP = false, Nullable<long> seed = R_NilValue){
  
  Environment base = Environment("package:base");
  
  if(seed.isNotNull()){
    Function set_seed_r = base["set.seed"];
    set_seed_r(seed);
  }
  
  IntegerVector index(seqC(0, K - 1));
  IntegerVector index_L(seqC(0, L - 1));
  long N = cluster_id.length();
  // observation-wise weights
  // each row is weights for each DP
  NumericMatrix w(K, L);
  //return List::create(Named("w") = w);
  for(int j = 0 ; j < K; ++j){
    if(L > 1){
      NumericVector b = rbeta(L, 1, M_beta);
      b[L - 1] = 1;
      w(j, 0) = b[0];
      
      for(int i = 1; i < L; i++){
        w(j, i) = b[i] * prodC(1.0 - as<NumericVector>(b[seqC(0,i-1)]));
      }
    }else{
      w(j, 0) = 1;
    }
  }
  
  //rownames(w) = index;
  //colnames(w) = index_L;
  
  
  // weights to choose which DP
  NumericVector pi(K);
  if(K > 1){
    NumericVector b = rbeta(K, 1, M_alpha);
    b[K - 1] = 1;
    pi[0] = b[0];
    for(int i = 1; i < K; i++){
      pi[i] = b[i] * prodC(1.0 - as<NumericVector>(b[seqC(0,i-1)]));
    }
  }else{
    pi[0] = 1;
  }
  pi.attr("names") = index;
  
  
  parameters["CDP"] = CDP;
  parameters["K"] = K;
  parameters["L"] = L;
  
  
  
  //if(strcmp(parameters["distribution"], "normal") == 0){
  
  
  List samples_parameters = DP_sampler(K * L, parameters);
  NumericMatrix sample_DP = samples_parameters["y"];
  
  
  long dp = parameters["p"];
  //group data
  List theta(K);
  for(int j = 0 ; j < K; ++j){
    
    NumericVector v = w(j,_);
    v.attr("names") = index_L;
    
    NumericMatrix data = sample_DP(Range(j * L, (j + 1) * L - 1), _);
    rownames(data) = seqC(0, L - 1);
    if(CDP){
      for(int i = 0; i < dp; ++i){
        data (_,i) = data (_,i) - sum(data(_,i) * v);
      }
    }
    theta[j] = List::create(Named("y") = data, Named("pi") = v);
  }
  //return List::create(Named("theta") = theta, Named("w") = w, Named("pi") = pi, Named("sample_DP") = sample_DP);
  
  CharacterVector unique_cluster_id = cluster_id[duplicated(cluster_id) == 0];
  //Rcout << unique_cluster_id << std::endl;
  /*List cluster_map;
   
   for(int i = 0 ; i < unique_cluster_id.length(); ++i){
   long sampled = sample(index, 1, true, pi)[0];
   String unique_cluster_id_ = unique_cluster_id[i];
   cluster_map.push_back(sampled, unique_cluster_id_);
   }*/
  
  IntegerVector cluster_map;
  for(int i = 0 ; i < unique_cluster_id.length(); ++i){
    long sampled = sample(index, 1, true, pi)[0];
    cluster_map.push_back(sampled);
    //String unique_cluster_id_ = unique_cluster_id[i];
    //cluster_map.push_back(sampled, unique_cluster_id_);
  }
  
  cluster_map.attr("names") = unique_cluster_id;
  
  IntegerVector cluster = cluster_map[cluster_id];
  
  
  IntegerVector with_cluster(cluster_id.length());
  //return List::create(Named("with_cluster") = with_cluster, Named("cluster") = cluster, Named("cluster_map") = cluster_map, Named("theta") = theta, Named("w") = w, Named("pi") = pi, Named("sample_DP") = sample_DP);
  
  
  for(int i = 0 ; i < cluster_id.length(); ++i){
    List d = theta[cluster[i]];
    NumericVector pi_individual = d["pi"];
    long sampled = sample(index_L, 1, true, pi_individual)[0];
    with_cluster[i] = sampled;
  }
  with_cluster.attr("names") = cluster_id;
  
  NumericMatrix samples(N, dp);
  for(int i = 0 ; i < N; ++i){
    List d = theta[cluster[i]];
    NumericMatrix data_individual = d["y"];
    samples(i, _) = data_individual(with_cluster[i], _);
  }
  rownames(samples) = with_cluster;
  //return List::create(Named("samples") = samples, Named("with_cluster") = with_cluster, Named("cluster") = cluster, Named("cluster_map") = cluster_map, Named("theta") = theta, Named("w") = w, Named("pi") = pi, Named("sample_DP") = sample_DP);
  
  
  return List::create(Named("w") = w, Named("cluster_id") = cluster_id, Named("parameters") = parameters, Named("M_alpha") = M_alpha, Named("M_beta") = M_beta, Named("samples") = samples, Named("within_cluster") = with_cluster, Named("cluster") = cluster, Named("cluster_map") = cluster_map, Named("p") = pi, Named("y") = theta);
  
}

// [[Rcpp::export]]
List update_nDP_normal(NumericMatrix X, List tau, double L_alpha = -1, double U_alpha = 2, double L_beta = -1, double U_beta = 2, double tol = 1e-40){
  //Rcout << 123;
  
  Environment base = Environment("package:base");
  
  //Function set_seed = base["set.seed"];
  //set_seed(123);
  
  
  List parameters = tau["parameters"];
  //Rcout << 123;
  double M_alpha = tau["M_alpha"];
  //Rcout << 123;
  double M_beta = tau["M_beta"];
  bool CDP = parameters["CDP"];
  
  //Rcout << 123;
  IntegerVector cluster_map = tau["cluster_map"];
  NumericMatrix samples = tau["samples"];
  //Rcout << 123 <<std::endl;
  IntegerVector within_cluster = tau["within_cluster"];
  IntegerVector cluster = tau["cluster"];
  NumericVector pi = tau["p"];
  CharacterVector cluster_id = tau["cluster_id"];
  List y = tau["y"];
  CharacterVector unique_cluster_id = cluster_map.names();
  
  
  long N = X.nrow();
  long N_subject = cluster_map.length();
  long p = as<int>(parameters["p"]);
  long n_obs_per_subject = max(table(cluster_id));
  long K = parameters["K"];
  long L = parameters["L"];
  NumericMatrix w = tau["w"];
  
  IntegerVector index(seqC(0, K - 1));
  IntegerVector index_L(seqC(0, L - 1));
  CharacterVector index_L_ch(seqC(0, L - 1));
  
  //Rcout << 123;
  
  Environment mvtnorm("package:mvtnorm");
  Function rmvnorm = mvtnorm["rmvnorm"];
  Function dmvnorm = mvtnorm["dmvnorm"];
  Environment cascsim("package:hbmem");
  Function rtgamma = cascsim["rtgamma"];
  //Environment matlib("package:matlib");
  
  Environment G = Rcpp::Environment::global_env();
  
  Function matrix_multiply = G["matrix_multiply"];
  Function make_symmetric = G["make_symmetric"];
  Function matrix_add = G["matrix_add"];
  Function solve = base["solve"];
  Environment stats("package:stats");
  Function dnorm = stats["dnorm"];
  
  
  
  
  
  if(strcmp(parameters["distribution"], "normal") == 0){
    
    //Function set_seed_r = base["set.seed"];
    //set_seed_r(std::floor(std::fabs(seed)));
    
    double sigma;
    NumericMatrix Sigma;
    if(p == 1){
      sigma = tau["sigma"];
    }else{
      Sigma = as<NumericMatrix>(tau["Sigma"]);
    }
    
    
    
    for (int j=0; j<cluster_map.length(); ++j)
    {
      //for(int j = 0; j < cluster_map.length(); ++j){
      //Rcout << j << std::endl;
      CharacterVector subject_id {unique_cluster_id[j]};
      LogicalVector subject_X_index = character_vector_equals(cluster_id, subject_id);
      NumericMatrix x_sub = row_matrix(X, subject_X_index);
      NumericVector log_prob(K);
      for(int k = 0; k < K; ++k){
        log_prob[k] = 0;
        List y_k = y[k];
        NumericMatrix theta_k = y_k["y"];
        NumericVector p_k = y_k["pi"];
        for(int i = 0; i < x_sub.rows(); ++i){
          double w_p = 0;
          NumericVector xn(L);
          for(int l = 0 ; l < L; ++l){
            if(p == 1){
              xn[l] = log(p_k[l]) + as<NumericVector>(dnorm(x_sub(i, 0),  theta_k(l, 0), sigma, true))[0];
            }else{
              xn[l] = log(p_k[l]) + as<NumericVector>(dmvnorm(x_sub(i, _),  theta_k(l, _), Sigma, true))[0];
            }
          }
          log_prob[k] += max(xn) + log(sum(exp(xn - max(xn))));
        }
        log_prob[k] += log(pi[k]);
      }
      
      NumericVector pr = exp(log_prob - max(log_prob));
      
      //pr = set_tol(pr, tol);
      
      pr = pr / sum(pr);
      //Rcout << pr << std::endl;
      cluster_map[j] = sample(index, 1, true, pr)[0];
    }
    cluster = cluster_map[cluster_id];
    
    for(int i = 0; i < N; ++i){
      long cl = cluster[i];
      NumericVector pr (L);
      List y_cl = y[cl];
      NumericMatrix theta_cl = y_cl["y"];
      if(p == 1){
        pr = dnorm(X(i, 0), theta_cl(_, 0), sigma, true);
      }else{
        for(int l = 0 ; l < L; ++l){
          pr[l] = as<NumericVector>(dmvnorm(X(i, _), theta_cl(l, _), Sigma, true))[0];
        }
      }
      //pr = set_tol(pr, tol);
      NumericVector y_cl_pi = y_cl["pi"];
      pr = pr + log(y_cl_pi);
      pr = pr - max(pr);
      pr = exp(pr) / sum(exp(pr));
      //pr = pr / sum(pr);
      within_cluster[i] = sample(index_L, 1, true, pr)[0];
    }
    
    
    IntegerVector nk = table(cluster_map);
    for(int i = 0; i < K; ++i){
      const char * n = std::to_string(i).c_str();
      if(!nk.containsElementNamed(n)){
        nk[n] = 0;
      }
    }
    CharacterVector p_names = pi.names();
    nk = nk[p_names];
    
    if (nk.length() > 1){
      NumericMatrix beta(nk.length() - 1, 2);
      beta(_,0) = as<NumericVector>(nk[seqC(0, (nk.length() - 2))]) + 1;
      IntegerVector rev_nk = rev(nk);
      IntegerVector rev_cumsum = cumsum(rev_nk);
      IntegerVector non_rev_cumsum = rev(rev_cumsum);
      beta(_,1) = as<NumericVector>(non_rev_cumsum[seqC(1, (nk.length() - 1))]) + M_alpha;
      
      NumericVector Vh(nk.length() - 1);
      for(int i = 0; i < nk.length() - 1; ++i){
        double new_beta_sample = rbeta(1, beta(i, 0), beta(i, 1))[0];
        while(new_beta_sample == 1.0){
          new_beta_sample = rbeta(1, beta(i, 0), beta(i, 1))[0];
        }
        Vh[i] = new_beta_sample;
      }
      Vh.push_back(1);
      pi = NumericVector(nk.length());
      pi[0] = Vh[0];
      NumericVector log_negative_1_cumsum = cumsum(log(1 - Vh));
      //Rcout << negative_1_cumprod << std::endl;
      for(int i = 1; i < nk.length(); ++i){
        pi[i] = Vh[i] * exp(log_negative_1_cumsum[i - 1]);
      }
      //return List::create(Named("a") = log_negative_1_cumsum[nk.length() - 2], Named("Vh") = Vh, Named("pi") = pi, Named("M_alpha") = M_alpha, Named("log_negative_1_cumsum") = log_negative_1_cumsum);
      //set_seed(123);
      M_alpha = as<NumericVector>(rtgamma(1, K - 1, -1 / (log_negative_1_cumsum[nk.length() - 2]), pow(N_subject, L_alpha), pow(N_subject, U_alpha)))[0];
    }else{
      pi[0] = 1;
      M_alpha = 0;
    }
    pi.names() = p_names;
    
    //return List::create(Named("M_alpha") = M_alpha, Named("nk") = nk, Named("pi") = pi, Named("cluster_map") = cluster_map, Named("cluster") = cluster, Named("within_cluster") = within_cluster);
    
    
    //set_seed(123);
    double M_beta_second = 0;
    for(int k = 0 ; k < K; ++k){
      ////////////////
      //set_seed(123);
      NumericVector Vh(L - 1);
      NumericVector p_l(L);
      if(L > 1){
        if(nk[k] > 0){
          LogicalVector cluster_index = (cluster == k); 
          IntegerVector nk_within = table(as<IntegerVector>(within_cluster[cluster_index]));
          
          for(int l = 0; l < L; ++l){
            const char * n = std::to_string(l).c_str();
            if(!nk_within.containsElementNamed(n)){
              nk_within[n] = 0;
            }
          }
          
          nk_within = nk_within[index_L_ch];
          //if (k == 4)
          //  return List::create(Named("index_L") = index_L, Named("cluster") = cluster, Named("nk_within") = nk_within, Named("K") = k);
          
          NumericMatrix beta(L - 1, 2);
          beta(_,0) = as<NumericVector>(nk_within[seqC(0, (nk_within.length() - 2))]) + 1;
          IntegerVector rev_nk_within = rev(nk_within);
          
          IntegerVector rev_cumsum_within = cumsum(rev_nk_within);
          IntegerVector non_rev_cumsum_within = rev(rev_cumsum_within);
          
          
          beta(_,1) = as<NumericVector>(non_rev_cumsum_within[seqC(1, (nk_within.length() - 1))]) + M_beta;
          //if (k == 4)
          //  return List::create(Named("beta") = beta, Named("nk_within") = nk_within, Named("K") = k);
          
          
          for(int i = 0; i < L - 1; ++i){
            double new_beta_sample = rbeta(1, beta(i, 0), beta(i, 1))[0];
            while(new_beta_sample == 1.0){
              new_beta_sample = rbeta(1, beta(i, 0), beta(i, 1))[0];
            }
            Vh[i] = new_beta_sample;
            M_beta_second += max_d(log(1 - Vh[i]), log(tol));
          }
          Vh.push_back(1);
          //return List::create(Named("Vh") = Vh);
          p_l[0] = Vh[0];
          NumericVector log_negative_1_cumsum = cumsum(log(1 - Vh));
          for(int i = 1; i < L; ++i){
            p_l[i] = Vh[i] * exp(log_negative_1_cumsum[i - 1]);
          }
          
        }else{
          Vh = rbeta(L - 1, 1, M_beta);
          Vh.push_back(1);
          p_l[0] = Vh[0];
          for(int i = 1; i < L; i++){
            M_beta_second += max_d(log(1 - Vh[i - 1]), log(tol));
            p_l[i] = Vh[i] * prodC(1.0 - as<NumericVector>(Vh[seqC(0,i-1)]));
          }
          //Rcout << k << std::endl;
          //return List::create(Named("Vh") = Vh, Named("p_l") = p_l);
        }
      }else{
        Vh.push_back(1);
        p_l[0] = 1;
      }
      //return List::create(Named("beta") = beta, Named("nk_within") = nk_within, Named("K") = k);
      
      Vh.names() = index_L;
      p_l.names() = index_L;
      List y_k = y[k];
      y_k["pi"] = p_l;
      y[k] = y_k;
      w(k, _) = p_l;
      //within_weights[k] = List::create(Named("Vh") = Vh, Named("p_l") = p_l);
    }
    
    //Rcout << M_beta_second << std::endl;
    // update M_beta
    ////////////
    //set_seed(123);
    M_beta = as<NumericVector>(rtgamma(1, K * (L - 1), -1/M_beta_second, pow(n_obs_per_subject, L_beta), pow(n_obs_per_subject, U_beta)))[0];
    //return List::create(Named("y") = y, Named("w") = w, Named("M_beta") = M_beta);
    if(p == 1){
      double inv_sigma_2 = std::pow(1 / sigma, 2);
      for(int k = 0 ; k < K ; ++k){
        //set_seed(123);
        List y_k = y[k];
        NumericMatrix y_k_l = y_k["y"];
        NumericVector y_k_pi = y_k["pi"];
        LogicalVector cluster_index = (cluster == k); 
        for(int l = 0 ; l < L; ++l){
          LogicalVector within_cluster_index = (within_cluster == l);
          LogicalVector select_index = logic_and(cluster_index, within_cluster_index);
          long n_h = count_if(select_index);
          if(n_h > 0){
            double var = 1 / (as<double>(parameters["inv_sigma_2"]) + inv_sigma_2 * n_h);
            
            double mean = ((colSums(row_matrix(X, select_index))[0] * inv_sigma_2) + (as<double>(parameters["mu"]) * as<double>(parameters["inv_sigma_2"]))) * var;
            //Rcout << mean << " " << sqrt(var) << std::endl;
            double sample = rnorm(1, mean, sqrt(var))[0];
            //set_seed(123);
            //Rcout << sample << std::endl;
            y_k_l(l, 0) = sample;
          }else{
            double var = parameters["sigma"];
            double mean = parameters["mu"];
            //Rcout << mean << " " << (var) << std::endl;
            double sample = rnorm(1, mean, var)[0];
            //set_seed(123);
            //Rcout << sample << std::endl;
            y_k_l(l, 0) = sample;
          }
        }
        
        if(CDP == true){
          for(int i = 0; i < p; ++i){
            y_k_l (_,i) = y_k_l (_,i) - sum(y_k_l(_,i) * y_k_pi);
          }
        }
        
        y_k["y"] = y_k_l;
        y[k] = y_k;
      }
    }else{
      NumericMatrix inv_sigma = solve(Sigma);
      for(int k = 0 ; k < K ; ++k){
        List y_k = y[k];
        NumericMatrix y_k_l = y_k["y"];
        NumericVector y_k_pi = y_k["pi"];
        LogicalVector cluster_index = (cluster == k);
        for(int l = 0 ; l < L; ++l){
          LogicalVector within_cluster_index = (within_cluster == l);
          LogicalVector select_index = logic_and(cluster_index, within_cluster_index);
          if(any(select_index)){
            NumericMatrix var = solve(matrix_add(parameters["inv_sigma"], matrix_mul_scalar(inv_sigma, count_if(select_index))));
            var = make_symmetric(var);
            NumericVector mean = matrix_multiply(matrix_add(matrix_multiply(colSums(row_matrix(X, select_index)), inv_sigma), matrix_multiply(parameters["mu"], parameters["inv_sigma"])), var);
            NumericVector sample = rmvnorm(1, mean, var);
            y_k_l(l,_) = sample;
          }else{
            NumericMatrix var = parameters["sigma"];
            NumericVector mean = parameters["mu"];
            NumericVector sample = rmvnorm(1, mean, var);
            y_k_l(l, _) = sample;
          }
        }
        if(CDP == true){
          for(int i = 0; i < p; ++i){
            y_k_l (_,i) = y_k_l (_,i) - sum(y_k_l(_,i) * y_k_pi);
          }
        }
        y_k["y"] = y_k_l;
        //Rcout << y_k_l << std::endl;
        y[k] = y_k;
      }
    }
    //return List::create(Named("y") = y,Named("M_alpha") = M_alpha, Named("nk") = nk, Named("pi") = pi, Named("cluster_map") = cluster_map, Named("cluster") = cluster, Named("within_cluster") = within_cluster);
    //return List::create();
    for(int i = 0 ; i < N; ++i){
      List d = y[cluster[i]];
      NumericMatrix data_individual = d["y"];
      samples(i, _) = data_individual(within_cluster[i], _);
    }
    
    tau["cluster"] = cluster;
    tau["cluster_map"] = cluster_map;
    tau["within_cluster"] = within_cluster;
    tau["M_alpha"] = M_alpha;
    tau["M_beta"] = M_beta;
    tau["y"] = y;
    tau["w"] = w;
    tau["samples"] = samples;
    
    return tau;
    
    
    //return List::create(Named("M_alpha") = M_alpha, Named("k") = K);
    //return List::create(Named("M_alpha") = M_alpha, Named("pi") = pi, Named("cluster_map") = cluster_map, Named("cluster") = cluster, Named("within_cluster") = within_cluster);
    
    
    
    //Rcout << "complete update cluster map" << std::endl;
    /*
     //update location parameters
     if(p == 1){
     double inv_sigma_2 = std::pow(1 / sigma, 2);
     for(int k = 0 ; k < K ; ++k){
     List y_k = y[k];
     NumericMatrix y_k_l = y_k["y"];
     NumericVector y_k_pi = y_k["pi"];
     for(int l = 0 ; l < L; ++l){
     LogicalVector cluster_index = (cluster == k); 
     LogicalVector within_cluster_index = (within_cluster == l);
     LogicalVector select_index = logic_and(cluster_index, within_cluster_index);
     if(any(select_index)){
     double var = 1 / (as<double>(parameters["inv_sigma_2"]) + inv_sigma_2 * count_if(select_index));
     double mean = ((colSums(row_matrix(X, select_index))[0] * inv_sigma_2) + (as<double>(parameters["mu"]) * as<double>(parameters["inv_sigma_2"]))) * var;
     double sample = rnorm(1, mean, sqrt(var))[0];
     y_k_l(l, 0) = sample;
     }else{
     double var = parameters["sigma"];
     double mean = parameters["mu"];
     double sample = rnorm(1, mean, sqrt(var))[0];
     y_k_l(l, 0) = sample;
     }
     }
     
     if(CDP == true){
     for(int i = 0; i < p; ++i){
     y_k_l (_,i) = y_k_l (_,i) - sum(y_k_l(_,i) * y_k_pi);
     }
     }
     
     y_k["y"] = y_k_l;
     y[k] = y_k;
     }
     }
     
     else{
     NumericMatrix inv_sigma = solve(Sigma);
     for(int k = 0 ; k < K ; ++k){
     List y_k = y[k];
     NumericMatrix y_k_l = y_k["y"];
     NumericVector y_k_pi = y_k["pi"];
     for(int l = 0 ; l < L; ++l){
     LogicalVector cluster_index = (cluster == k);
     LogicalVector within_cluster_index = (within_cluster == l);
     LogicalVector select_index = logic_and(cluster_index, within_cluster_index);
     if(any(select_index)){
     NumericMatrix var = solve(matrix_add(parameters["inv_sigma"], matrix_mul_scalar(inv_sigma, count_if(select_index))));
     var = make_symmetric(var);
     NumericVector mean = matrix_multiply(matrix_add(matrix_multiply(colSums(row_matrix(X, select_index)), inv_sigma), matrix_multiply(parameters["mu"], parameters["inv_sigma"])), var);
     NumericVector sample = rmvnorm(1, mean, var);
     y_k_l(l,_) = sample;
     }else{
     NumericMatrix var = parameters["sigma"];
     NumericVector mean = parameters["mu"];
     NumericVector sample = rmvnorm(1, mean, var);
     y_k_l(l, _) = sample;
     }
     }
     if(CDP == true){
     for(int i = 0; i < p; ++i){
     y_k_l (_,i) = y_k_l (_,i) - sum(y_k_l(_,i) * y_k_pi);
     }
     }
     y_k["y"] = y_k_l;
     y[k] = y_k;
     }
     }
     
     
     for(int i = 0 ; i < N; ++i){
     List d = y[cluster[i]];
     NumericMatrix data_individual = d["y"];
     samples(i, _) = data_individual(within_cluster[i], _);
     }
     
     tau["cluster"] = cluster;
     tau["cluster_map"] = cluster_map;
     tau["within_cluster"] = within_cluster;
     tau["M_alpha"] = M_alpha;
     tau["M_beta"] = M_beta;
     tau["y"] = y;
     tau["samples"] = samples;
     
     return tau;
     }
     return tau;
     */
  }
  return tau;
}