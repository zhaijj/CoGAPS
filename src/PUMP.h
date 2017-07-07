#ifndef PUMP_H_
#define PUMP_H_

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <limits>
#include <cmath>
#include <limits>
#include <stdexcept>

double get_max(vector<double> x);

vector<double>  vectorDiff(vector<double>  X, vector<double>  Y);

double get_dot(vector<double> X, vector<double> Y);

int which_min(vector<double> x);

vector<int> subsetAndDuplicate(vector<int> p, int tmp_p);

vector<double> get_column(vector<vector<double> > x, int column);

vector<int> get_unique(vector<int> x);

int get_match_counts(vector<int> x, int a);

vector<int> patternMarkers(vector<vector<double> > A_mat, vector<vector<double> > P_mat);

#endif
