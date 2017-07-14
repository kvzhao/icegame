#ifndef _LATTICE_H_INCLUDED
#define _LATTICE_H_INCLUDED

#include <iostream>
#include <vector>
#include <string>
#include <fstream>
//#include <math>
using namespace std;

#include "sample.hpp"


typedef vector<int> vtr1d;
typedef vector<vtr1d> vtr2d;
typedef vector<vtr2d> vtr3d;

struct Coordinates2d{
  double x,y;
};

struct Coordinates3d{
  double x,y,z;
};


class Lattice{
  protected:
    Lattice(){}
    //virtual ~Hamiltonian() = 0;
  public:
    string name;
    virtual void output_lattice() = 0;
    virtual void construct_lattice() = 0;
    virtual void find_neighbors() = 0;
    virtual void show_lattice() = 0;
    virtual int periodic_bd(int, int) = 0;
    int L; // lattice_size
    int N; // Num_sites
    int R; // Num_replicas
    int Num_neighbors;
    vector< vector<int> > NN;
    vector< vector<int> > NNN;
    vector< vector<int> > NNNN;
    vector< vector<int> > NNNN2;
    vector< vector<int> > ANNN;
    
    vtr3d site_index3d;
    vtr2d site_index2d;
    vector<int> sub;
    Point3d* site3d;
    Point2d* site2d;
    Coordinates2d* coor;

    vector<double> Bloch_phase_cos;
    vector<double> Bloch_phase_sin;

    vector< vector<double> > neutron_cos;
    vector< vector<double> > neutron_sin;
    vector< vector<double> > neutron_Q;
};

class Square : public Lattice{
  
  public:
    //Point2d* site;
    Square();
    Square(int, int&, int);    
    ~Square();
    // new function
    void init_from_size(int lattice_size, int num_neighbors); 
    void init(INFO info);
    void construct_lattice();
    void find_neighbors();
    void find_neighbors2();
    void find_neighbors_ANNN();
    void show_lattice();
    void output_lattice();
    int periodic_bd(int, int);
};

class Kagome : public Lattice{
  
  public:
    //vtr3d site_index;
    //vector<Point3d> site;
    //vector<Coordinates> coor;
    //Coordinates coor;
    //Point3d* site;
    //Coordinates* coor;
    //vtr2d* Bravais;
    //Point* site;
    
    Kagome(int, int&, int);    
    ~Kagome();
    //vector<double> Bloch_phase;
    void neutron_initialization(vector<vector<double> > );
    void Hstate_Bloch_phase();
    void construct_lattice();
    void find_neighbors();
    void find_neighbors2();
    void find_neighbors3();
    void show_lattice();
    void output_lattice();
    int periodic_bd(int, int);

};




#endif
