// Try to construct the simplest MC simulation code with OOP
// initialization -> update -> measure -> output -> binning/bootstrap
// measure energy only
// without magnetic field

#ifndef _HAMILTONIAN_H_INCLUDED
#define _HAMILTONIAN_H_INCLUDED

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <math.h> // for pow()
#include <numeric>
using namespace std;


#include "sample.hpp"
#include "lattice.hpp"
#include "rng.hpp"


class Hamiltonian{
  protected:
    Hamiltonian(){}
    //virtual ~Hamiltonian() = 0;
  public:
    virtual void MCstep(Sample*, Lattice*) = 0;
    //virtual double spin_correlation(Sample*, Lattice*, vector<double> &) = 0;
    //virtual void set_field(Sample*, double, double, double) = 0;
    double compare_state(Sample*, Sample*);
    int L;
    int N;
    int R;
    int Num_neighbors;

};

class Ising : public Hamiltonian{
  
  public:
    double J1;
    double J2;
    double J3a;
    double J3b;
     
    Ising(INFO, int);
    ~Ising();
    
    void set_J1(double);
    void set_J2(double);
    void set_J3a(double);
    void set_J3b(double);
    double SE(int, Sample*, Lattice*);
    double SE2(int, Sample*, Lattice*);
    double SE3(int, Sample*, Lattice*);

    double ME(int, Sample*);

    void MCstep(Sample*, Lattice*);
   
    bool SSF(Sample*, Lattice*);
    bool SSF2(Sample*, Lattice*);
    bool SSF3(Sample*, Lattice*);

    double total_energy(Sample* , Lattice*);
    double magnetization(Sample*, Lattice*);
    static double spin_correlation(Sample*, Lattice*,vector<double> &);
    
    static double stagger_magnetization(Sample*, Lattice*);

    static double obs_energy(Sample*, Lattice*, vector<double> &);
    static double obs_magnetization(Sample*, Lattice*, vector<double> &);
    static double obs_acceptanceSSF(Sample*, Lattice*, vector<double> &);

    void initialization(Sample*, Lattice*, int);
    void initialize_observable(Sample*, Lattice*, double, vector<double>);
};

class Kagome_ice : public Hamiltonian{

 public:
    double J1;
    double J2;
    double J3a;
    double J3b;
    //static double* Bloch_phase;
    //vector<double> product_hs;
    //static const Coordinates3d easyaxis[3];
    //static void set_easyaxis();
    /*
    class spin_easyaxis{
      public:
        spin_easyaxis();
        double sq2, sq3;
        vector<vector<double> > sublattice;
        //double x0, y0, z0, x1, y1, z1, x2, y2, z2;
    };
    static spin_easyaxis* easyaxis;
    */
    static void get_easyaxis(int, vector<double> &);


    Kagome_ice(INFO, int);
    ~Kagome_ice();
    
    void set_J1(double);
    void set_J2(double);
    void set_J3a(double);
    void set_J3b(double);
    double SE(int, Sample*, Lattice*);
    double SE2(int, Sample*, Lattice*);
    double SE3(int, Sample*, Lattice*);

    double M_111(int, Sample*, Lattice*);
    double M_112(int, Sample*, Lattice*);
    double M_110(int, Sample*, Lattice*);

    double magnetization_111(Sample*, Lattice*);
    double magnetization_112(Sample*, Lattice*);
    double magnetization_110(Sample*, Lattice*);
    
    static double defect_density(Sample*, Lattice*, vector<double> &);
    static double FP_order(Sample*, Lattice*, vector<double> &);
    static double H_order(Sample*, Lattice*, vector<double> &);
    static double domain_order(Sample*, Lattice*, vector<double> &);
    static double magnetization(Sample*, Lattice*, vector<double> &);
    
    static double neutron_map(Sample*, Lattice*, vector<double> &);
    
    double ME(int, Sample*, Lattice*);

    void MCstep(Sample*, Lattice*);
   
    bool SSF(Sample*, Lattice*);
    bool SSF2(Sample*, Lattice*);
    bool SSF3(Sample*, Lattice*);

    double total_energy(Sample* , Lattice*);
    double spin_correlation(Sample*, Lattice*,vector<double> &);
    

    static double obs_energy(Sample*, Lattice*, vector<double> &);
    static double obs_magnetization111(Sample*, Lattice*);
    static double obs_magnetization112(Sample*, Lattice*);
    static double obs_magnetization110(Sample*, Lattice*);
    static double obs_acceptanceSSF(Sample*, Lattice*, vector<double> &);

    void initialization(Sample*, Lattice*, int);
    void initialize_observable(Sample*, Lattice*, double, vector<double>);
    void set_product_hs(Sample*, double, double, double);
    void defect_configuration(Sample*, Lattice*, string, int);
    void degenerate_vortex(Sample*, Lattice*, int);
    void degenerate_qx(Sample*, Lattice*, int);
    void monopole_vortex(Sample*, Lattice*, int);
    void anti_vortex(Sample*, Lattice*, int);
    //void Hstate_Bloch_phase(Lattice*);
};

class Square_ice : public Hamiltonian{
  
  public:
    double J1;
     
    // empty constructor
    //
    Square_ice(INFO, int);
    Square_ice(){};
    ~Square_ice();
    
    void init(INFO info);
    void set_J1(double);
    double SE(int, Sample*, Lattice*);

    double ME(int, Sample*);

    void MCstep(Sample*, Lattice*);
   
    bool SSF(Sample*, Lattice*);

    double total_energy(Sample* , Lattice*);
    double magnetization(Sample*, Lattice*);
    
    static double obs_energy(Sample*, Lattice*, vector<double> &);
    static double obs_magnetization(Sample*, Lattice*, vector<double> &);
    static double obs_acceptanceSSF(Sample*, Lattice*, vector<double> &);
    static double defect_density(Sample*, Lattice*, vector<double> &);
    void initialization(Sample*, Lattice*, int);
    void initialize_observable(Sample*, Lattice*, double, vector<double>);
};


#endif
