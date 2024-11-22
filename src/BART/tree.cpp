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

#include "tree.h"


//--------------------------------------------------
//functions
std::ostream& operator<<(std::ostream& os, const tree& t)
{
  tree::cnpv nds;
  t.getnodes(nds);
  os << nds.size() << std::endl;
  for(size_t i=0;i<nds.size();i++) {
    os << nds[i]->nid() << " ";
    os << nds[i]->getv() << " ";
    os << nds[i]->getc() << " ";
    os << nds[i]->gettheta() << std::endl;
  }
  return os;
}
std::istream& operator>>(std::istream& is, tree& t)
{
  size_t tid,pid; //tid: id of current node, pid: parent's id
  std::map<size_t,tree::tree_p> pts;  //pointers to nodes indexed by node id
  size_t nn; //number of nodes
  
  t.tonull(); // obliterate old tree (if there)
  
  //read number of nodes----------
  is >> nn;
  if(!is) {
    //cout << ">> error: unable to read number of nodes" << endl;
    return is;
  }
  
  //read in vector of node information----------
  std::vector<node_info> nv(nn);
  for(size_t i=0;i!=nn;i++) {
    is >> nv[i].id >> nv[i].v >> nv[i].c >> nv[i].theta;
    if(!is) {
      //cout << ">> error: unable to read node info, on node  " << i+1 << endl;
      return is;
    }
  }
  //first node has to be the top one
  pts[1] = &t; //careful! this is not the first pts, it is pointer of id 1.
  t.setv(nv[0].v); t.setc(nv[0].c); t.settheta(nv[0].theta);
  t.p=0;
  
  //now loop through the rest of the nodes knowing parent is already there.
  for(size_t i=1;i!=nv.size();i++) {
    tree::tree_p np = new tree;
    np->v = nv[i].v; np->c=nv[i].c; np->theta=nv[i].theta;
    tid = nv[i].id;
    pts[tid] = np;
    pid = tid/2;
    // set pointers
    if(tid % 2 == 0) { //left child has even id
      pts[pid]->l = np;
    } else {
      pts[pid]->r = np;
    }
    np->p = pts[pid];
  }
  return is;
}