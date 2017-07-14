// Try to construct the simplest MC simulation code with OOP
// initialization -> update -> measure -> output -> binning/bootstrap
// measure energy only
// without magnetic field

#ifndef _MEASUREMENT_H_INCLUDED
#define _MEASUREMENT_H_INCLUDED

#include <iostream>
#include <vector>
using namespace std;


#include "sample.hpp"
#include "observable.hpp"
#include "lattice.hpp"

class Measurement{
  
  //friend class SquareIsing;
  
  public: // no need to be modified
    Measurement();
    ~Measurement();
    void go(Sample*, Lattice*);
    void set_counter(int);
    void set_files();
    void set_files_closed();
    void set_files_info(Sample*, int);
    void insert_single_obs(string);
    void insert_array_obs(string); 
  public: // need to be modified when adding new obervable
    //Hamiltonian* model;
    vector<Observable*> machine;
    vector<Observable*>::iterator iter;
    int machine_size;
    int latest;

};

#endif
