#ifndef _GIBBSSAMPLERTRANSFORMATION_H_
#define _GIBBSSAMPLERPATTERN_H_

#include "GibbsSamplerMap.h"
#include <vector>
#include <cmath>

class GibbsSamplerPattern : public GibbsSamplerMap {
  protected:
    int _whichPattern;
    // vector of case status (assume 0 is non-case) i.e. {0, 0, 0, 1, 1, 1}
    std::vector<int> _caseStatus;

  public:
    GibbsSamplerPattern(unsigned long nEquil, unsigned long nSample, unsigned int nFactor,
                        double alphaA, double alphaP, double nMaxA, double nMaxP,
                        unsigned long nIterA, unsigned long nIterP,
                        double max_gibbsmass_paraA, double max_gibbsmass_paraP,
                        unsigned long long atomicSize,
                        char label_A, char label_P, char label_D, char label_S,
                        vector<vector<double> > &DVector, vector<vector<double> > &SVector,
                        const string &simulation_id,
                        vector <vector <double> >  &parameters, char the_fixed_matrix, int whichPattern);

    ~GibbsSamplerPattern() {};
    
    // expect that some transformation of the pattern will be linear so
    // allow for passing an a priori function that transforms the pattern
    // to a linear pattern
    // for example, if we expect the pattern to be logistic, we can pass a logistic
    // function to transform the data i.e.
    // GibbSampPatt.update_pattern(GibbsSampPatt::&logit);
    void update_pattern(std::vector<double>(*transformation)(std::vector<double>));

    // transformations
    // logistic transformation (log(p / 1 - p))
    std::vector<double> logit(std::vector<double> data);
    // identity transformation (pattern is already linear)
    std::vector<double> identity(std::vector<double> data);
};

#endif
