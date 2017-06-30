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

// Amatrix: genes x pattern_weights
// Pmatrix: patterns x samples
// threshold: "unique" - we'll just use this
// lp: don't want
// full: don't want
// return: partition of genes into list

// [[Rcpp::export()]]
double max(NumericVector X){
  double max_val = -std::numeric_limits<double>::infinity();
  for (int ii=0; ii < X.size(); ii++){
    if (max_val > X[ii]){ max_val = X[ii]; }
  }
  max_val;
}

// [[Rcpp::export()]]
List patternMarkers(NumericMatrix A, NumericMatrix P){
  // get scaling factors for A matrix
  NumericVector p_scales(P.nrow());
  for (ii=0; ii < P.nrow(); ii++){
    p_scales[ii] = max(P(ii, _));
  }
  // scale A matrix, assuming P matrix is not scaled
  for (int ii=0; ii < A.ncol(); ii++){
    for (int jj=0; jj < A.nrow(); jj++){
      A(ii, jj) = A(ii, jj) * p_scales[ii];
    }
  }

}

