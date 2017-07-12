#ifndef __SQICE_HPP__
#define __SQUCE_HPP__

// cpp standard 
#include <iostream>
#include <vector>
#include <math.h>
#include <numeric>
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

//#define DEBUG ON

enum DIR {RIGHT, DOWN, LEFT, UP, LOWER_NEXT, UPPER_NEXT, LOWER_RIGHT, LOWER_LEFT, UPPER_LEFT, UPPER_RIGHT, NOWAY};


class SQIceGame {
    public:
        // Init Constructor
        SQIceGame (INFO info);

        void InitModel();
        void SetTemperature(double T);
        void MCRun(int mcSteps);
        vector<double> Draw(int dir_idx);
        void IceMove(int act_idx);

        vector<double> Metropolis();
        inline void FlipTrajectory() {flip_along_traj(traj_norepeat);};
        int Start(int init_site);
        void Reset();

        void UpdateConfig();

        inline int GetAgentSite() {return get_agent_site();};
        vector<double> GetStateTMap();
        vector<double> GetCanvasMap(); 
        vector<double> GetEnergyMap();
        vector<double> GetDefectMap();

        inline int GetUpdatedCounter() {return updated_counter;};
        inline vector<int> GetAcceptedLen() {return accepted_looplength;};
        inline vector<int> GetTrajectory() {return traj_norepeat;};

        bool TimeOut();

        void update_ice_config();
        /// this function used fliiping ice config according to states?

        void flip_site(int site);
        void flip_agent_site();
        void flip_along_traj(const vector<int> &traj);
        void flip_along_short_traj(const vector<int> &traj);

        // return new agent site
        int go(DIR dir);
        DIR how_to_go(int site);

        // return site
        int icemove(bool);

        void update_state_to_config();
        void restore_config_to_state();

        void reset_maps();

        int get_neighbor_site_by_direction(DIR dir);
        DIR get_direction_by_sites(int site, int next_site);
        vector<int> get_neighbor_sites();
        vector<int> get_neighbor_spins();
        vector<int> get_neighbor_candidates(bool same_spin);

        void set_agent_site(int site);
        /* get funcs */
        int  get_agent_site();
        int  get_agent_spin();
        int  get_spin(int site);

        // Test function
        void TEST();

    private:
        // private functions
        double _cal_energy_of_state(const vector<int> &s);
        double _cal_energy_of_site(const vector<int> &s, int site);
        double _cal_defect_density_of_state(const vector<int> &s);
        int _cal_config_difference();
        int inline _pdb(int site, int d, int l) {return ((site + d) % l + l) % l;};

        void long_loop_algorithm();

        void _print_vector(const vector<int> &v);
        double _cal_mean(const vector<int> &s);
        double _cal_stdev(const vector<int> &s);
        bool _is_visited(int site);
        bool _is_traj_continuous();
        bool _is_traj_intersect();
        bool _is_start_end_meets(int site);
        bool _is_long_loop();
        bool _is_short_loop();

        // Physical System
        INFO sim_info;
        Square latt;
        Square_ice model;
        Sample ice_config;

        double config_mean; 
        double config_stdev;

        unsigned int L, N;
        double kT;
        double h1_t, h2_t, h3_t;
        double J1;
        vector<double> mag_fields;

        unsigned int updated_counter;
        unsigned int steps_counter;
        vector<int> accepted_looplength;

        // RL intefaces
        int agent_site;
        int init_agent_site;
        vector<int> traj_sites;
        vector<int> traj_norepeat;
        vector<int> traj_spins;
        vector<DIR> action_list;
        vector<int> sites_counter;

        // MAPS
        /* Current states */ 
        vector<int> state_0;
        vector<int> state_t;
        vector<int> state_tp1;

        vector<int> agent_map;
        vector<int> canvas_traj_map;
        vector<int> canvas_spin_map;
        vector<double> energy_map;
        vector<double> defect_map;
        vector<int> diff_map;

        // utilities
        Timer tt;

};

// Converter for std::vector to python list
template <class T>
struct Vec2List {
    static PyObject* convert (const std::vector<T> &vec) {
        boost::python::list *l = new boost::python::list();
        for (size_t i = 0; i < vec.size(); ++i)
            (*l).append(vec[i]);
        return l->ptr();
    }
};

BOOST_PYTHON_MODULE(libicegame)
{
    class_<INFO>("INFO", init<int, int, int, int, int, int, int, int>())
    ;

    to_python_converter<std::vector<int, class std::allocator<int> >, Vec2List<int> >();
    to_python_converter<std::vector<double, class std::allocator<double> >, Vec2List<double> >();

    class_<SQIceGame>("SQIceGame", init<INFO>())
        .def("init_model", &SQIceGame::InitModel)
        .def("set_temperature", &SQIceGame::SetTemperature)
        .def("start", &SQIceGame::Start)
        .def("reset", &SQIceGame::Reset)
        .def("mc_run", &SQIceGame::MCRun)
        .def("draw", &SQIceGame::Draw)
        .def("icemove", &SQIceGame::IceMove)
        .def("get_agent_site", &SQIceGame::GetAgentSite)
        .def("get_canvas_map", &SQIceGame::GetCanvasMap)
        .def("get_state_t_map", &SQIceGame::GetStateTMap)
        .def("get_energy_map", &SQIceGame::GetEnergyMap)
        .def("get_defect_map", &SQIceGame::GetDefectMap)
        .def("metropolis", &SQIceGame::Metropolis)
        .def("flip_trajectory", &SQIceGame::FlipTrajectory)
        .def("update_config", &SQIceGame::UpdateConfig)
        .def("get_updated_counter", &SQIceGame::GetUpdatedCounter)
        .def("get_accepted_length", &SQIceGame::GetAcceptedLen)
        .def("get_trajectory", &SQIceGame::GetTrajectory)
        .def("timeout", &SQIceGame::TimeOut)

        .def("TEST", &SQIceGame::TEST)
    ;

}

#endif