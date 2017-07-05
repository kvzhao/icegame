// Try to construct the simplest MC simulation code with OOP
// initialization -> update -> measure -> output -> binning/bootstrap
// measure energy only
// without magnetic field

#ifndef _BINNING_H_INCLUDED
#define _BINNING_H_INCLUDED

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <math.h> // for pow()
using namespace std;
double Avg_vector(int, vector<double>);
double Var_vector(int, vector<double>);
double Std_vector(int, vector<double>);


#include "sample.hpp"
#include "rng.hpp"

class Bin_analysis{
  protected:
    int Num_bins;
    int size_bin;
    int counter_replica;
    int Num_replicas;
    int Num_sites;
    int Num_tests;
    string name;
    double temperature;
    double field;
    int size_data;
  public:
    //Bin_analysis(string, int, int, INFO); 
    void set_parameters(string, int, int, INFO);
    ~Bin_analysis();
    vector< vector<double> > data;
    vector< vector<double> > data2;
    vector< vector<double> > data4;
    vector< vector<double> > data_response;
    vector< vector<double> > data_abs;
    vector< vector<double> > data_Binder;
    vector< vector<double> > data_bootstrap;
    vector< vector<double> > data_response_bootstrap;
    vector< vector<double> > data_Binder_bootstrap;
    //vector< vector<double> > data_list;

    ifstream input;
  
    void input_open();  
    void input_close();
    double get_temperature();
    double get_field();

    void binning();  //do average only, for dataline(list of data)

    double get_mean_value(int);
    double get_squared_value(int);
    double get_abs_mean(int);
    virtual double get_response_function(int) = 0;
    
    double get_Binder_cumulant(int);

    void bootstrap(int);
    double get_mean_bootstrap(int);
    double get_mean_error(int);
    double get_reponse_bootstrap(int);
    double get_response_error(int);
    double get_Binder_bootstrap(int);
    double get_Binder_error(int);

};

class Simple_bin : public Bin_analysis{

  public:
    Simple_bin(string, int, int, INFO);
    double get_response_function(int);
};

class Energy_bin : public Bin_analysis{

  public:
    Energy_bin(string, int, int, INFO);
    double get_response_function(int);

};

class Magnetization_bin : public Bin_analysis{

  public:
    Magnetization_bin(string, int, int, INFO);
    double get_response_function(int);
 
};



#endif
