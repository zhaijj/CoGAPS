#ifndef _GIBBSSAMPLERSEQ_H_
#define _GIBBSSAMPLERSEQ_H_

#include<iostream>
#include<fstream>
#include <string>
#include <vector>
#include <algorithm>
#include "GibbsSampler.h"
#include<limits>

using namespace gaps;

class GibbsSamplerSeq : public GibbsSampler {
  public:

    // ******************** CONSTRUCTOR ********************************************
    GibbsSamplerSeq(unsigned long nEquil, unsigned long nSample, unsigned int nFactor,
                    double alphaA, double alphaP, double nMaxA, double nMaxP,
                    unsigned long nIterA, unsigned long nIterP,
                    double max_gibbsmass_paraA, double max_gibbsmass_paraP,
                    unsigned long long atomicSize,
                    char label_A, char label_P, char label_D, char label_S,
                    vector<vector<double> > &DVector, vector<vector<double> > &SVector,
                    const string &simulation_id);

    ~GibbsSamplerMap() {};
};
#endif
