
#ifndef _ANALYSIS_H_INCLUDED
#define _ANALYSIS_H_INCLUDED

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector> 

#include "sample.hpp"
#include "binning.hpp"

class Analysis{
  private:
    string name;
    int Num_observables;
    ofstream output;
    int Num_replicas;
    int mode; // 0: simple_bin, 1: energy_bin, 2: magnetization_bin
  public:
    Analysis(string, int, int, INFO);
    void set_parameters(string, int, int, INFO);
    void simple_analysis(INFO);
    void energy_analysis(INFO);
    void magnetization_analysis(INFO);
    void correlation_analysis(INFO);

    ~Analysis();
    vector<double> temperature;
    vector<double> field; 
    vector< vector<double> > summ_mean;
    vector< vector<double> > summ_error;
    vector< vector<double> > summ_response;
    vector< vector<double> > summ_Binder;

    vector< vector<double> > summ_bootstrap;
    vector< vector<double> > summ_bootstrap_error;
    vector< vector<double> > summ_bootstrap_response;
    vector< vector<double> > summ_bootstrap_response_error;
    vector< vector<double> > summ_bootstrap_Binder;
    vector< vector<double> > summ_bootstrap_Binder_error;

    void output_open();
    void output_close();
    void output_write();
    void output_write_no_bootstrap();

};


#endif
