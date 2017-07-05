#ifndef __SQICE_HPP__
#define __SQUCE_HPP__

// cpp standard 
#include <iostream>
#include <vector>
#include <math.h>
#include <algorithm>
#include <string>

// monte carlo libraries
#include "sample.hpp"
#include "observable.hpp"
#include "hamiltonian.hpp"
#include "measurement.hpp"
#include "heattreatment.hpp"
#include "binning.hpp"
#include "timer.hpp"
#include "analysis.hpp"

// boost.python intefaces
#include <boost/python.hpp>
#include <boost/python/stl_iterator.hpp>
using namespace boost::python; 

/* Function Naming Convention:
    Upper Case for interface
    lower case used as member functions
*/
// 

enum DIR {UP, DOWN, LEFT, RIGHT, UPPER_LEFT, UPPER_RIGHT, LOWER_LEFT, LOWER_RIGHT};


class SQIceGame {
    public:
        // Init Constructor
        SQIceGame (INFO info);

        //

        void InitModel();
        void SetTemperature(double T);
        void MCRun(int mcSteps);

        // 

        void update_ice_config();
        /// this function used fliiping ice config according to states?

        void flip_site(int site);
        void flip_agent_site();

        // return new site
        int go(DIR dir);

        void reset_maps();

        /* get funcs */
        int get_neighbor_site_by_direction(DIR dir);
        vector<int> neighbor_sites();
        vector<int> neighbor_spins();

        void set_agent_site(int site);
        int  get_agent_site();

        // Test function
        void TEST();


    private:
        // private functions
        double _cal_energy_of_state(vector<int> &s);
        double _cal_defect_density_of_state(vector<int> &s);
        int inline _pdb(int site, int d) {return ((site + d) % L + L) % L;};
        void _print_vector(vector<int> v);

        // Physical System
        INFO sim_info;
        Square latt;
        Square_ice model;
        Sample ice_config;

        unsigned int L, N;
        double kT;
        double h1_t, h2_t, h3_t;
        double J1;
        vector<double> mag_fields;

        /* Current states */ 
        vector<int> state_0;
        vector<int> state_t;
        vector<int> state_tp1;

        vector<int> traj_sites;
        vector<int> traj_spins;

        unsigned int updated_counter;

        // RL intefaces
        int agent_site;

        vector<int> agent_map;

        // utilities
        Timer tt;

};

BOOST_PYTHON_MODULE(libicegame)
{
    class_<INFO>("INFO", init<int, int, int, int, int, int, int, int>())
    ;
    class_<SQIceGame>("SQIceGame", init<INFO>())
        .def("init_model", &SQIceGame::InitModel)
        .def("set_temperature", &SQIceGame::SetTemperature)
        .def("mc_run", &SQIceGame::MCRun)

        .def("TEST", &SQIceGame::TEST)
    ;

}

#endif