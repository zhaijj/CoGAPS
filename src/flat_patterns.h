#ifndef FLAT_PATTERNS_H_
#define FLAT_PATTERNS_H_

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <limits>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <numeric>

double mean(vector<double> x);

vector <double> i_minus(vector<double> x, int idx);

double pless_than(vector<double> x, int idx, double eps);

vector <int> find_flat_patterns(vector<vector<double> > mat, double flat_eps, double p_eps);

#endif
