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




