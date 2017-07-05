// Try to construct the simplest MC simulation code with OOP
// initialization -> update -> measure -> output -> binning/bootstrap
// measure energy only
// without magnetic field

#ifndef _HEATTREATMENT_H_INCLUDED
#define _HEATTREATMENT_H_INCLUDED

#include <iostream>
//#include <fstream>
#include <vector>
//#include <string>
//#include <sstream>
//#include <time.h>
#include <math.h> // for pow()
using namespace std;


#include "sample.hpp"
#include "hamiltonian.hpp"
#include "lattice.hpp"

class HeatTreatment{
  protected:
    int Num_replicas;
    int Num_thermalization;
    int temperature_period;
  public:
    void import_info(const INFO&);
    virtual void thermalization(Hamiltonian*, Lattice*, vector<Sample *> ) = 0; 
};

class parallel_tempering : public HeatTreatment{
  
  public:
    vector<double> acceptance;
    
    parallel_tempering();  
    parallel_tempering(const INFO&);  

    void init(INFO);
    
    void reset_accptance();

    void thermalization(Hamiltonian*, Lattice*, vector<Sample *>);
    
    void thermalization2(Hamiltonian*, Lattice*, vector<Sample *>);

    void tempering(int, vector<Sample *> );
};

class field_tempering : public HeatTreatment{

  public:
    int index;
    
    vector<double> acceptance;

    field_tempering(INFO, int);

    void reset_accptance();

    void thermalization(Hamiltonian*, Lattice*, vector<Sample *>);

    void thermalization2(Hamiltonian*, Lattice*, vector<Sample *>);

    void tempering(int, vector<Sample *>);
};

class simulated_annealing : public HeatTreatment{
  private:
    int current_replica;
  public:
    simulated_annealing(INFO);
    void thermalization(Hamiltonian*, Lattice*, vector<Sample *> );

    void annealing(int n, int current, vector<Sample*> vs);
   
};

#endif
