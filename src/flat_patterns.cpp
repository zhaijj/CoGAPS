#include <iostream>       // for use with standard I/O
#include <string>         // for string processing
#include <fstream>        // for output to files
#include <limits>         // for extracting numerical limits of C++
#include <iterator>
#include <numeric>

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
#include "flat_patterns.h"
#include "PUMP.h"

// does all the atomic space to matrix conversion
// and sampling actions.
#include <RcppArmadillo.h>
// ------------------------------------------------------

/* NEEDED: validate this code works */


//namespace bpo = boost::program_options;
using namespace std;
using namespace gaps;
using std::vector;
using namespace Rcpp;

double cmean(vector<double> x){
  double ave = 0;
  for (int ii=0; ii < x.size(); ii++){
    ave += x[ii];
  }
  ave /= x.size();
  return(ave);
}

vector <double> i_minus(vector<double> x, int idx){
  vector <double> x2 = x;
  x2.erase(x2.begin() + idx);
  return(x2);
}

double pless_than(vector<double> x, int idx, double eps){
  vector <double> x_iless = i_minus(x, idx);
  double x_i = x[idx];
  vector <double> lessthan(x_iless.size());
  for (int ii=0; ii < x_iless.size(); ii++){
    x_iless[ii] -= x_i;
    if (abs(x_iless[ii]) < eps){
      lessthan[ii] = 1;
    } else {
      lessthan[ii] = 0;
    }
  }
  return(cmean(lessthan));
}

vector <double> standardize_row(vector <double> x){
  double x_mean = cmean(x);
  for (int ii=0; ii < x.size(); ii++){
    x[ii] /= x_mean;
  }
  return(x);
}

vector <int> find_flat_patterns(vector<vector< double> > mat, double flat_eps, double p_eps){
  vector <int> flats(mat.size());
  vector<double> ps(mat.size());
  for (int ii=0; ii < mat.size(); ii++){
    vector<double> row = mat[ii];
    vector<double> ps_row(row.size());
    row = standardize_row(row);
    // standardize pattern to have mean = 1
    for (int jj=0; jj < row.size(); jj++){
      ps_row[jj] = pless_than(row, jj, flat_eps);
    }
    ps[ii] = cmean(ps_row);
  }
  for (int ii=0; ii < ps.size(); ii++){
    if (ps[ii] > p_eps){
      flats[ii] = 1;
    } else {
      flats[ii] = 0;
    }
  }
  return(flats);
}

