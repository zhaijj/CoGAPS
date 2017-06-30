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

// [[Rcpp::export()]]
bool findStringInList(string gene, List gene_list){
  bool found = false;
  for (int ii=0; ii < gene_list.size(); ii++){
    string gene_now = gene_list[ii];
    if (gene == gene_now){
      found = true;
    }
  }
  return found;
}

// listP is list of lists

// [[Rcpp::export()]]
List which_in(string gene, List listP){
  List inn;
  for (int ii=0; ii < listP.size(); ii++){
    bool found = findStringInList(gene, listP[ii]);
    if (found){
      inn.push_back(ii);
    }
  }
  return inn;
}

// [[Rcpp::export()]]
List assignGenes(List genes, List listP){
  List assigns;
  for (int ii=0; ii < genes.size(); ii++){
    assigns.push_back(which_in(genes[ii], listP));
  }
  assigns.attr("names") = genes;
  return assigns;
}

