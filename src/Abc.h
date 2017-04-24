#ifndef _ABC_H_
#define _ABC_H_

#include <vector>
#include <cmath>
#include <string>
#include <algorithm>
#include <RcppArmadillo.h>
#include <stdexcept>

class Abc {
    private:
        // MH parameters
        double _delta;

        // ABC parameters
        double _epsilon;

        // value 
        Rcpp::NumericVector _theta;
        std::vector<double> _pattern;

        // priors
        double _prior_mean;
        double _prior_sd;

        // data
        Rcpp::NumericVector _T; // time of recording
        Rcpp::NumericMatrix _D; // data

        // prior density
        std::string _prior_choice;
        Rcpp::NumericVector _prior(Rcpp::NumericVector param);

        // proposal distribution
        std::string _proposal_choice;
        Rcpp::NumericVector _proposal(); // draw from proposal
        Rcpp::NumericVector _proposal(Rcpp::NumericVector param1, // density
                                      Rcpp::NumericVector param2);

        // updated matrices
        Rcpp::NumericMatrix A_update;
        Rcpp::NumericMatrix P_update;

        // epsilon mcmc per Bortot et al 2007
        bool _epsilon_mcmc;
        double _epsilon_rate;
        double _epsilon_prior();
        Rcpp::NumericVector _epsilon_prior(double param);
        double _epsilon_propose(); // draw new epsilon
        Rcpp::NumericVector _epsilon_propose(double param1,  // density of epsilon
                                             double param2);

    public:
        // constructor
        Abc(std::vector<std::vector<double> >& data, 
            std::vector<double> timeRecorded,
            std::string prior="normal",
            std::string proposal="normal",
            bool epsilon_mcmc=false,
            double delta=10.0,
            double epsilon=100.0,
            double epsilon_rate=3.0,
            double prior_mean=0.0,
            double prior_sd=10.0);

        // propose new theta
        void propose(Rcpp::NumericMatrix A, Rcpp::NumericMatrix P);

        // get new pattern and amplitude
        std::vector<double> pattern();
        std::vector<double> amplitude();

        // get parameters
        Rcpp::NumericVector theta();
        Rcpp::NumericVector epsilon();

};

#endif