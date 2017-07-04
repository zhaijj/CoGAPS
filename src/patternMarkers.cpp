#include <iostream>       // for use with standard I/O
#include <string>         // for string processing
#include <fstream>        // for output to files
#include <limits>         // for extracting numerical limits of C++


// ------ incorporated to use Cogaps_options ------------
#include <vector>
#include <iomanip>
#include <boost/algorithm/string.hpp>
// ------------------------------------------------------
#include "randgen.h"   // for incorporating a random number generator.
#include "Matrix.h"  // for incorporating a Matrix class
#include "AtomicSupport.h"  // for incorporating an Atomic class
#include "GAPSNorm.h"  // for incorporating calculation of statistics in cogaps.
#include "GibbsSampler.h" // for incorporating the GibbsSampler which
// does all the atomic space to matrix conversion
// and sampling actions.
#include <RcppArmadillo.h>
// ------------------------------------------------------

//namespace bpo = boost::program_options;
using namespace std;
using namespace gaps;
using std::vector;
using namespace Rcpp;

boost::mt19937 rng(43);

// https://github.com/Bioconductor-mirror/CoGAPS/blob/b1b986984c001bd7562c9e32e219a5095fef3a2a/src/GibbsSampler-atomic.cpp#L50
// double ** Amatrix - convert to/from NumericMatrix?

// Amatrix: genes x pattern_weights
// Pmatrix: patterns x samples
// threshold: "unique" - we'll just use this
// lp: don't want
// full: don't want
// return: partition of genes into list

// for debugging from R purposes
vector < vector < double > > NumericMat2Vecs(NumericMatrix X){
  vector < vector < double > > out;
  int nrow = X.nrow();
  int ncol = X.ncol();
  out.resize(nrow);

  for (int ii=0; ii < nrow; ii++){
    out[ii].resize(ncol);
  }

  for (int ii=0; ii < nrow; ii++){
    Rcpp::NumericVector tmp = X(ii, _);
    for (int jj=0; jj < ncol; jj++){
      double tmp_element = tmp[jj];
      out[ii][jj] = tmp_element;
    }
  }

  return out;
}

vector < std::string > listOfStrings2Vec(List X){
  vector < std::string >  out;
  int len = X.size();
  out.resize(len);
  for (int ii=0; ii < len; ii++){
    out[ii] = as<std::string>(X[ii]);
  }
  return out;
}

vector < vector <std::string> > listOfStringLists2Vecs(List X){
  vector < vector <std::string > > out;
  int len_outer = X.size();
  out.resize(len_outer);
  for(int ii=0; ii < len_outer; ii++){
    vector <std::string> tmp = X[ii];
    int len_inner = tmp.size();
    out[ii].resize(len_inner);
    for (int jj=0; jj < len_inner; jj++){
      out[ii][jj] = tmp[jj];
    }
  }
  return out;
}

double get_max(vector<double> X){
  double max_val = -std::numeric_limits<double>::infinity();
  for (int ii=0; ii < X.size(); ii++){
    if (X[ii] > max_val){ max_val = X[ii]; }
  }
 return max_val;
}

bool findStringInList(string gene, vector<std::string>  gene_list){
  bool found = false;
  for (int ii=0; ii < gene_list.size(); ii++){
    //string gene_now = gene_list[ii];
    if (gene == gene_list[ii]){
      found = true;
    }
  }
  return found;
}

vector<int> which_in(string gene, vector<vector<std::string> >  listP){
  vector < int >  inn;
  for (int ii=0; ii < listP.size(); ii++){
    vector <std::string> tmp_listp = listP[ii];
    bool found = findStringInList(gene, tmp_listp);
    if (found){
      inn.push_back(ii);
    }
  }
  return inn;
}

vector<vector<int> > which_in_lists(vector<std::string>  genes, vector<vector<std::string> > listP){
  vector < vector < int > > inn(genes.size());
  for(int jj=0; jj < genes.size(); jj++){
    string gene = genes[jj];
    for (int ii=0; ii < listP.size(); ii++){
      vector<std::string> tmp = listP[ii];
      bool found = findStringInList(gene, tmp);
      if (found){
        inn[jj].push_back(ii);
      }
    }
  }
  return inn;
}

vector<double> vectorDiff(vector<double>  X, vector<double>  Y){
  vector < double >  out(X.size());
  for (int ii=0; ii < X.size(); ii++){
    out[ii] = X[ii] - Y[ii];
  }
  return out;
}

double get_dot(vector<double>  X, vector<double>  Y){
  double sum = 0;
  for(int ii=0; ii < X.size(); ii++){
    sum += X[ii] * Y[ii];
  }
  return sum;
}

int which_min(vector<double>  x){
  int idx = 0;
  double val = std::numeric_limits<double>::infinity();
  for (int ii=0; ii < x.size(); ii++){
    if (x[ii] < val){
      val = x[ii];
      idx = ii;
    }
  }
  return idx;
}

vector <int> idx_sort(vector <double> x){
  vector <int> y(x.size());
  for (int ii=0; ii < x.size(); ii++){
    y[ii] = ii;
  }
  for (int ii=0; ii < x.size(); ii++){
    for(int jj=ii+1; jj < x.size(); jj++){
      if (x[ y[ii] ] > x[ y[jj] ]){
        int tmp_y;
        tmp_y = y[ii];
        y[ii] = y[jj];
        y[jj] = tmp_y;
      }
    }
  }
  return y;
}

// [[Rcpp::export()]]
List patternMarkersC(NumericMatrix A, NumericMatrix P){
  // get scaling factors for A matrix
  vector < vector < double > > A_mat = NumericMat2Vecs(A);
  vector < vector < double > > P_mat = NumericMat2Vecs(P);

  vector<double> p_scales(P_mat.size());
  for (int ii=0; ii < P_mat.size(); ii++){
    p_scales[ii] = get_max(P_mat[ii]);
  }
  // scale A matrix, assuming P matrix is not scaled
  for (int ii=0; ii < A_mat[0].size(); ii++){
    for (int jj=0; jj < A_mat.size(); jj++){
      A_mat[ii][jj] = A_mat[ii][jj] * p_scales[ii];
    }
  }
  // get scaled A matrix
  vector < vector < double > >  Arowmax(A_mat.size(), vector <double>(A_mat[0].size()));
  for (int ii=0; ii < A_mat.size(); ii++){
    for (int jj=0; jj < A_mat[0].size(); jj++){
      Arowmax[ii][jj] = A_mat[ii][jj]/get_max(A_mat[ii]);
    }
  }
  // get maxes
  vector<double>  pmax(A_mat.size());
  for(int ii=0; ii < A_mat.size(); ii++){
    pmax[ii] = get_max(A_mat[ii]);
  }
  printf ("set pmaxes\n");
  // this is equivalent to setting lp = NA
  vector <vector <double> > sstat (A_mat.size(), vector <double> (A_mat[0].size()));
  // sstat.attr("dimnames") = A.attr("dimnames");
  vector <vector <double > > ssranks (A_mat.size(), vector <double> (A_mat[0].size()));
  // ssranks.attr("dimnames") = A.attr("dimnames");
  vector <vector <std::string > >  ssgenes (A_mat.size(), vector < std::string >(A_mat[0].size()));

  for (int ii=0; ii < A_mat[0].size(); ii++){
    vector <double> lp(A_mat[0].size(), 0);
    lp[ii] = 1;
    for (int jj=0; jj < A_mat.size(); jj++){
      vector <double> tmp_A = A_mat[jj];
      vector <double> vv = vectorDiff(tmp_A, lp);
      sstat[jj][ii] = sqrt(get_dot(vv, vv));
      printf ("SSTAT %d, %d: %f\n", jj, ii, sstat[jj][ii]);
    }
    vector <double> tmp_sstat = sstat[ii];
    vector <int> tmp_sstat_idx(tmp_sstat.size());
    tmp_sstat_idx = idx_sort(tmp_sstat);

    int counter = 0;
    for (int jj=0; jj < tmp_sstat_idx.size(); jj++){
      ssranks[tmp_sstat_idx[jj]][ii] = counter;
      printf("Ranks %d, %d: %d\n", tmp_sstat_idx[jj], ii, counter);
      counter += 1;
    }
    // how do we propogate rownames in cogaps? doesn't impact writing this yet ...
    //ssgenes(_, ii) = tmp_sstat.attr("names");
  } // end for loop over number of patterns

  // // below is equivalent to patternMarkers' unique option
  // NumericVector pIndx(ssranks.nrow());
  // for (int ii=0; ii < ssranks.nrow(); ii++){
  //   pIndx[ii] = which_min(ssranks(ii, _));
  // }

  // NumericVector pIndx_unique = clone(pIndx);
  // std::unique(pIndx_unique.begin(), pIndx_unique.end());
  // std::sort(pIndx_unique.begin(), pIndx_unique.end());
  // List gBYp;
  // for (int ii=0; ii < pIndx_unique.size(); ii++){
  //   gBYp.push_back(pIndx[ pIndx == pIndx_unique[ii] ]);
  // }

  // List ssgenes_th;
  // for (int ii=0; ii < A.ncol(); ii++){
  //   CharacterVector ssgenes_tmp = ssgenes(_, ii);
  //   ssgenes_th.push_back(ssgenes_tmp[which_in_lists(ssgenes_tmp, gBYp[ii])]);
  // }
  // return ssgenes_th;
  List out(4);
  return out;
}

