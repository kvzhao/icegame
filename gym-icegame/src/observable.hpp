// Try to construct the simplest MC simulation code with OOP
// initialization -> update -> measure -> output -> binning/bootstrap
// measure energy only
// without magnetic field

#ifndef _OBSERVABLE_H_INCLUDED
#define _OBSERVABLE_H_INCLUDED

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cstring>
#include <sstream>
#include <math.h> // for pow()
using namespace std;

#include "sample.hpp"
#include "lattice.hpp"

class Observable{
  
  friend class Measurement;
  
  protected:
    Observable(){}
    ~Observable(){}
    int counter_replica;
    void output_open();
    void output_temperature(Sample* );
    void output_field(Sample* , int);
    void output_data(double);
    void output_close();
    void output_dataline(vector<double> &); 

  public:
    void set_obs_ptr(double); 
    virtual void measure(Sample*, Lattice*) = 0;
    string name;
    ofstream output;
    void set_counter(int);
    double (*obs_ptr)(Sample*, Lattice*);
    double (*obs_aptr)(Sample*, Lattice*, vector<double>&);
};

class single_value_obs : public Observable{
  public:
    single_value_obs(string);
    void measure(Sample*, Lattice*);
};

class single_array_obs : public Observable{
  public:
    vector<double> data;
    single_array_obs(string);
    void measure(Sample*, Lattice*);
};

#endif
