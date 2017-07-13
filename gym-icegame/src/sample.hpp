// Try to construct the simplest MC simulation code with OOP
// initialization -> update -> measure -> output -> binning/bootstrap
// measure energy only
// without magnetic field

#ifndef _SAMPLE_H_INCLUDED
#define _SAMPLE_H_INCLUDED

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>

#include "rng.hpp"

using namespace std;



struct Point2d{
  int x,y;
};

struct Point3d{
  int x,y,z;
};


class INFO{
  public:
      // init constructor
    // (L, N, num_neighbors. num_replica, num_bins, steps_per_bin, tempering_periods, num_termalization)
    INFO(int, int, int, int, int, int, int, int); // 8 ints
    INFO(const INFO&);

    int lattice_size;
    int Num_sites;
    int Num_neighbors;
    int Num_replicas;
    int Num_MCsteps;
    int Num_bins;
    int Num_steps_per_bin;
    int tempering_period;
    int Num_thermalization;

    void print_INFO();
};
//INFO MC_info;

class Sample{
  public:
    Sample() {};
    Sample(INFO);
    ~Sample();
  private:
    double temperature;
    vector<double> field;
    //int initial_state;
    //int Num_sites;
    //int L;
	  
  public:
    //int Num_sites;
    int N;
    int L;
    int ctr_output;
    void init(INFO);
    vector<int> Ising;
    vector<double> correlation;

    std::vector<int> get_Ising() {return Ising;};

    void get_configuration(int);
    
    void output_configuration(string);
    void output_configuration(string, int);
    void output_configuration(string, int, int);
    
    void set_energy(double);
    void resize_order(int);
    void resize_magnetization(int);
    void set_magnetization(int, double);
    void set_temperature(double);
    void resize_field(int);
    void set_field(double, int);
    double get_temperature();
    double get_field(int);
    double energy_update;
    vector<double> magnetization_update;
    vector<double> order_update;
    double acceptance_SSF;
    vector<double> product_hs;
    //void initialization();
};


#endif
