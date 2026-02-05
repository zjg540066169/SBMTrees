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
 *  - To make it easier to compile, I move the function definitions in the separate 
 *  .cpp file to this file, and merge them with declaration.
 *
 *  These modifications comply with the terms of the GNU General Public License 
 *  version 2 (GPL-2).
 */


#include "bartfuns.h"




//--------------------------------------------------
//make xinfo = cutpoints
void makexinfo(size_t p, size_t n, double *x, xinfo& xi, size_t numcut)
{
  int* nc = new int[p];
  for(size_t i=0; i<p; ++i) nc[i]=numcut;
  makexinfo(p, n, x, xi, nc);
  delete [] nc;
}

void makexinfo(size_t p, size_t n, double *x, xinfo& xi, int *nc)
{
   double xinc;

   //compute min and max for each x
   std::vector<double> minx(p,INFINITY);
   std::vector<double> maxx(p,-INFINITY);
   double xx;
   for(size_t i=0;i<p;i++) {
      for(size_t j=0;j<n;j++) {
         xx = *(x+p*j+i);
         if(xx < minx[i]) minx[i]=xx;
         if(xx > maxx[i]) maxx[i]=xx;
      }
   }
   //make grid of nc cutpoints between min and max for each x.
   xi.resize(p);
   for(size_t i=0;i<p;i++) {
      xinc = (maxx[i]-minx[i])/(nc[i]+1.0);
      xi[i].resize(nc[i]);
      for(size_t j=0;j<nc[i];j++) xi[i][j] = minx[i] + (j+1)*xinc;
   }
}
//--------------------------------------------------
//--------------------------------------------------
//compute n and \sum y_i for left and right give bot and v,c
void getsuff(tree& x, tree::tree_p nx, size_t v, size_t c, xinfo& xi, dinfo& di, size_t& nl, double& syl, size_t& nr, double& syr)
{
   double *xx;//current x
   nl=0; syl=0.0;
   nr=0; syr=0.0;

   for(size_t i=0;i<di.n;i++) {
      xx = di.x + i*di.p;
      if(nx==x.bn(xx,xi)) { //does the bottom node = xx's bottom node
         if(xx[v] < xi[v][c]) {
               nl++;
               syl += di.y[i];
          } else {
               nr++;
               syr += di.y[i];
          }
      }
   }

}
//lh, replacement for lil that only depends on sum y.
double lh(size_t n, double sy, double sigma, double tau)
{
   double s2 = sigma*sigma;
   double t2 = tau*tau;
   double k = n*t2+s2;
   return -.5*log(k) + ((t2*sy*sy)/(2.0*s2*k));
}
//--------------------------------------------------

//--------------------------------------------------
//compute n and \sum y_i for left and right bots
void getsuff(tree& x, tree::tree_p l, tree::tree_p r, xinfo& xi, dinfo& di, size_t& nl, double& syl, size_t& nr, double& syr)
{
   double *xx;//current x
   nl=0; syl=0.0;
   nr=0; syr=0.0;

   for(size_t i=0;i<di.n;i++) {
      xx = di.x + i*di.p;
      tree::tree_cp bn = x.bn(xx,xi);
      if(bn==l) {
         nl++;
         syl += di.y[i];
      }
      if(bn==r) {
         nr++;
         syr += di.y[i];
      }
   }
}
//--------------------------------------------------

//--------------------------------------------------

//--------------------------------------------------





//--------------------------------------------------

//--------------------------------------------------


double log_sum_exp2(std::vector<double>& v){
    double mx=v[0],sm=0.;
    for(size_t i=0;i<v.size();i++) if(v[i]>mx) mx=v[i];
    for(size_t i=0;i<v.size();i++){
      sm += exp(v[i]-mx);
    }
    return mx+log(sm);
}

//--------------------------------------------------
//draw variable splitting probabilities from Dirichlet (Linero, 2018)
void draw_s(std::vector<size_t>& nv, std::vector<double>& lpv, double& theta, rn& gen){
  size_t p=nv.size();
// Now draw s, the vector of splitting probabilities
  std::vector<double> _theta(p);
  for(size_t j=0;j<p;j++) _theta[j]=theta/(double)p+(double)nv[j];
  //gen.set_alpha(_theta);
  lpv=gen.log_dirichlet(_theta);
}

//--------------------------------------------------
//draw Dirichlet sparsity parameter from posterior using grid
void draw_theta0_bartfuns(bool const_theta, double& theta, std::vector<double>& lpv,
		 double a, double b, double rho, rn& gen){
  // Draw sparsity parameter theta_0 (Linero calls it alpha); see Linero, 2018
  // theta / (theta + rho ) ~ Beta(a,b)
  // Set (a=0.5, b=1) for sparsity
  // Set (a=1, b=1) for non-sparsity
  // rho = p usually, but making rho < p increases sparsity
  if(!const_theta){
    size_t p=lpv.size();
    double sumlpv=0.,lse;
    
    std::vector<double> lambda_g (1000,0.);
    std::vector<double> theta_g (1000,0.);
    std::vector<double> lwt_g (1000,0.);
    for(size_t j=0;j<p;j++) sumlpv+=lpv[j];
    for(size_t k=0;k<1000;k++){
      lambda_g[k]=(double)(k+1)/1001.;
      theta_g[k]=(lambda_g[k]*rho)/(1.-lambda_g[k]);
      double theta_log_lik=lgamma(theta_g[k])-(double)p*lgamma(theta_g[k]/(double)p)+(theta_g[k]/(double)p)*sumlpv;
      double beta_log_prior=(a-1.)*log(lambda_g[k])+(b-1.)*log(1.-lambda_g[k]);
//      cout << "SLP: " << sumlogpv << "\nTLL: " << theta_log_lik << "\nBLP: " << beta_log_prior << '\n';
      lwt_g[k]=theta_log_lik+beta_log_prior;      
    }
    
    double mx=lwt_g[0],sm=0.;
    for(size_t i=0;i<lwt_g.size();i++) if(lwt_g[i]>mx) mx=lwt_g[i];
    for(size_t i=0;i<lwt_g.size();i++){
      sm += std::exp(lwt_g[i]-mx);
    }
    lse= mx+log(sm);
    for(size_t k=0;k<1000;k++) {
      lwt_g[k]=exp(lwt_g[k]-lse);
//      cout << "LWT: " << lwt_g[k] << '\n';
    }
    gen.set_wts(lwt_g);    
    theta=theta_g[gen.discrete()];
  } 
}

