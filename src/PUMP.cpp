#include <iostream>       // for use with standard I/O
#include <string>         // for string processing
#include <fstream>        // for output to files
#include <limits>         // for extracting numerical limits of C++
#include <iterator>

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
#include "PUMP.h"
// does all the atomic space to matrix conversion
// and sampling actions.
#include <RcppArmadillo.h>
// ------------------------------------------------------

//namespace bpo = boost::program_options;
using namespace std;
using namespace gaps;
using std::vector;
using namespace Rcpp;

// https://github.com/Bioconductor-mirror/CoGAPS/blob/b1b986984c001bd7562c9e32e219a5095fef3a2a/src/GibbsSampler-atomic.cpp#L50
// double ** Amatrix - convert to/from NumericMatrix?

// Amatrix: genes x pattern_weights
// Pmatrix: patterns x samples
// threshold: "unique" - we'll just use this
// lp: don't want
// full: don't want
// return: partition of genes into list

// for debugging from R purposes
// vector < vector < double > > NumericMat2Vecs(NumericMatrix X){
//   vector < vector < double > > out;
//   int nrow = X.nrow();
//   int ncol = X.ncol();
//   out.resize(nrow);

//   for (int ii=0; ii < nrow; ii++){
//     out[ii].resize(ncol);
//   }

//   for (int ii=0; ii < nrow; ii++){
//     Rcpp::NumericVector tmp = X(ii, _);
//     for (int jj=0; jj < ncol; jj++){
//       double tmp_element = tmp[jj];
//       out[ii][jj] = tmp_element;
//     }
//   }

//   return out;
// }

double get_max(vector<double> X){
  double max_val = -std::numeric_limits<double>::infinity();
  for (int ii=0; ii < X.size(); ii++){
    if (X[ii] > max_val){ max_val = X[ii]; }
  }
 return max_val;
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


vector<int> subsetAndDuplicate(vector<int> p, int tmp_p){
  vector<int> out;
  printf ("P size: %zd\n", p.size());
  for (int ii=0; ii < p.size(); ii++){
    if (p[ii] == tmp_p){
      out.push_back(ii);
    } 
  }
  return out;
}

vector<double> get_column(vector<vector<double> > x, int column){
  vector<double> out(x.size());
  for (int ii=0; ii < x.size(); ii++){
    out[ii] = x[ii][column];
  }
  return out;
}


vector<int> get_unique(vector<int> x){
  vector<int> a = x;
  vector<int> out;
  std::sort(a.begin(), a.end());
  std::unique_copy(a.begin(), a.end(), std::back_inserter(out));
  return out;
}

int get_match_counts(vector<int> x, int a){
  int out = 0;
  for (int ii=0; ii < x.size(); ii++){
    if (x[ii] == a){ out += 1; }
  }
  return out;
}

vector<int> patternMarkers(vector<vector<double> > A_mat,
                           vector<vector<double> > P_mat,
                           vector<int> flat_pats){
  // get scaling factors for A matrix
  // vector < vector < double > > A_mat = NumericMat2Vecs(A);
  // vector < vector < double > > P_mat = NumericMat2Vecs(P);

  vector<double> p_scales(P_mat.size());
  for (int ii=0; ii < P_mat.size(); ii++){
    p_scales[ii] = get_max(P_mat[ii]);
  }
  // scale A matrix, assuming P matrix is not scaled
  for (int ii=0; ii < A_mat[0].size(); ii++){
    for (int jj=0; jj < A_mat.size(); jj++){
      A_mat[jj][ii] = A_mat[jj][ii] * p_scales[ii];
    }
  }

  // get scaled A matrix
  vector < vector < double > >  Arowmax(A_mat.size(), vector <double>(A_mat[0].size()));
  for (int ii=0; ii < A_mat.size(); ii++){
    for (int jj=0; jj < A_mat[0].size(); jj++){
      // t(Arowmax) is needed bc apply returns pattern x genes, but genes x pattern is what is desired
      Arowmax[ii][jj] = A_mat[ii][jj]/get_max(A_mat[ii]); 
    }
  }
  // get maxes

  vector<double>  pmax(A_mat.size());
  for(int ii=0; ii < A_mat.size(); ii++){
    pmax[ii] = get_max(A_mat[ii]);
  }
  // this is equivalent to setting lp = NA
  vector <vector <double> > sstat (A_mat.size(), vector <double> (A_mat[0].size()));
  int counter = 0;

  vector<int> mins(Arowmax.size());

  for (int ii=0; ii < Arowmax[0].size(); ii++){
    // if a pattern is flat across samples we skip it
    if (flat_pats[ii] == 1){
      continue;
    }
    vector <double> lp(A_mat[0].size(), 0);
    lp[ii] = 1;
    for (int jj=0; jj < Arowmax.size(); jj++){
      vector <double> tmp_A = Arowmax[jj];
      vector <double> vv = vectorDiff(tmp_A, lp);
      sstat[jj][ii] = sqrt(get_dot(vv, vv));
    }
  } // end for loop over number of patterns

  for (int ii=0; ii < mins.size(); ii++){
    vector<double> tmp = sstat[ii];
    mins[ii] = which_min(tmp);
  }
  return mins;
}


/*
  Note: the trick works by asking how many times does a gene get assigned to the pattern
  assigned by the mean matrix (Amean).  This means you have to do this at the end.

  Want to iterate over the set of snapshots, compute patternMarker for each sample A',
  ask if it is the same as Amean.


  CoGAPS computes means/vars online ...
  - for each iteration, get normed matrices and compute pattern markers (see cogapsR.cpp # 440)
  - return vector of length n.genes where each element is the assigned pattern
  - store pattern markers (class variable?) until end of sampling (see GibbsSampMapler.h)
  - want matrix of integers: genes x pattern, entries are the number of times a gene is assigned
  to a pattern
  - compare each A' patternMarker to Amean (cogapsR.cpp # 504ish)
  - for a given gene (row) divide the row element corresponding to the pattern from
  Amean by the # of iterations
  - return vector of gene PUMPs

*/
