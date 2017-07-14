// Try to construct the simplest MC simulation code with OOP
// initialization -> update -> measure -> output -> binning/bootstrap
// measure energy only
// without magnetic field

#ifndef _CAL_H_INCLUDED
#define _CAL_H_INCLUDED


#include <iostream>
#include <fstream>
#include <vector>
#include <string>
//#include <sstream>
//#include <boost/random.hpp> // for mt19937
//#include <boost/random/mersenne_twister.hpp>
//#include <boost/random/uniform_01.hpp>
//#include <time.h>
#include <math.h> // for pow()
using namespace std;
//double Avg_vector(int, vector<double>);
//double Var_vector(int, vector<double>);
//double Std_vector(int, vector<double>);

//boost::mt19937 rng(time(NULL));
//boost::uniform_01<boost::mt19937> uni01_sampler(rng);
//string make_filename( const string&, int, const string&);

#include "sample.hpp"
#include "observable.hpp"
#include "hamiltonian.hpp"
#include "measurement.hpp"
#include "heattreatment.hpp"
#include "binning.hpp"
#include "timer.hpp"
#include "analysis.hpp"
//#include "rng.hpp"

INFO MC_info;

int main();

#endif



