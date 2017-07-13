#ifndef __SISIM_HPP__
#define __SISIM_HPP__

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <math.h> // for pow()
#include <algorithm>

// monte carlo libraries
#include "sample.hpp"
#include "observable.hpp"
#include "hamiltonian.hpp"
#include "measurement.hpp"
#include "heattreatment.hpp"
#include "binning.hpp"
#include "timer.hpp"
#include "analysis.hpp"

// boost.python interfaces
#include <boost/python.hpp>
#include <boost/python/stl_iterator.hpp>
using namespace boost::python;

/// Naming Convention
// * Camel for external function call, like python
// * CallForPython --> call_from_python
// * 

class SquareIceGame {
    public:
        // Constructor 
        SquareIceGame(INFO in);

        // Initialization Constructor
        // -- Physical Parameters set functions
        void setMagneticField(double h1, double h2, double h3);
        void setCoupling(double J1);
        void setTemperature(double T);

        /* States Transition Function*/ 
        //  Initialize every samples
        void InitConfigs();
        // reset physical configs, heating again, necessary?
        //void reset_configs();

        void MCRun(int mcSteps);

        int CreateDefect(int site);
        int RecoverDefect();
        void Reset();
        int InitAgent(int site);

        double Start(int site);

        // state transition
        void go_to_S0();            // S0 <- St
        void go_to_St();            // S0 -> St

        void update_to_newstate();  // St -> St+1

        int FlipSite(int site);
        int FlipOnSite();
        int StayOnSite();

        // action: return reward
        double Act(int action);
        int whichActionToFlipSite(int site);
        int whichActionToWalkSite(int site);
        int dual_action_flip(int action);
        int dual_action_walk(int action);
        int popActionQueue();
        int popShortLoopActionQueue();

        // <S0, St>
        double MetropolisJudgement();


        // Used for loop algorithm
        bool isShortClosed();
        bool isLongClosed();
        int popVisitedQueue();
        int popTrajectoryQueue();
        // Depreciated
        bool ShortLoopUpdate();

        // metroplis algorithm, return rewards vector
        // r[0]: acceptance status
        // r[1]: fraction of difference
        // r[2]: delta energy
        // r[3]: weights of boltzman factor
        // ...

        /* Move Functions */

        // compute difference between phys and virt, and save to diff_configs
        //
        // Different between S0 and St
        inline double getDiffRatio() {return difference_counter/N;};
        double getDiffEnergy();
        inline int getDiffCounter() {return difference_counter;};
        inline int getNumUpdates() { return num_S0_updates;};

        inline double getS0DefectDensity() {return calculate_S0_defect_density();};
        inline double getStDefectDensity() {return calculate_St_defect_density();};


        // used for state observation
        std::vector<int> getStConfigs();
        std::vector<int> getStNextConfigs();
        std::vector<int> getS0Configs();
        std::vector<int> getdSConfigs();
        int getS0Spin(int site);
        int getStSpin(int site);
        inline int getAgentSpin();
        inline int getAgentSite() {return agent_site;};
        inline int getDefectSpin() {return defect_spin;};
        inline int getDefectSite() {return defect_site;};
        // get spin on St
        int getSpin(int site);

        inline std::vector<int> getAgentPosMap() {return agent_pos_map;};
        inline std::vector<int> getSublattMap() {return sublatt_map;};
        inline std::vector<int> getVisitedMap() {return visited_map;};
        inline std::vector<int> getDefectMap() {return defect_map;};
        inline std::vector<int> getActionMap() {return action_map;};
        inline std::vector<int> getUpSpin() {return up_spin;};
        inline std::vector<int> getDownSpin() {return down_spin;};
        inline std::vector<int> getTrajectoryMap() {return trajectory_map;};
        inline std::vector<int> getTrajectory() {return trajectory_queue;};
        inline std::vector<int> getActionQueue() {return action_queue;};

        /* */

        /* Observation Functions */
        std::vector<int> getNeighborSites(int site);
        std::vector<int> getNeighborSpins(int site);

        // combined methods
        std::vector<int> getNeighborSitesAndSpins(int site);

        // Utilities functions
        void show_info() {MC_info.print_INFO();};

    private:
        // Private Functions // 
        int inline PBD (int site, int d) {return ((site+d) % L + L) % L;}
        
        // flip the site where agent is on, return spin after flipped
        
        // move the direction
        bool go_right();
        bool go_left();
        bool go_up();
        bool go_down();
        bool go_next_up();
        bool go_next_down();

        bool flip_right();
        bool flip_left();
        bool flip_up();
        bool flip_down();
        bool flip_next_up();
        bool flip_next_down();

        /*** 
         *
	        NN[p][0] = site_index2d[xp][y];
    	    NN[p][1] = site_index2d[x][yp];
    	    NN[p][2] = site_index2d[xm][y];
    	    NN[p][3] = site_index2d[x][ym];
    	    NNN[p][0] = site_index2d[xp][yp];
     	    NNN[p][1] = site_index2d[xm][yp];
        	NNN[p][2] = site_index2d[xm][ym];
    	    NNN[p][3] = site_index2d[xp][ym];
        ***/

        bool flip_next_upper_right();
        bool flip_next_upper_left();
        bool flip_next_lower_left();
        bool flip_next_lower_right();

        void update_agent_site(int new_site);
        void seperate_up_or_down(int site);

        double calculate_S0_defect_density();
        double calculate_St_defect_density();

        void reset_maps();

        // Private Data //
        //
        double default_reward;
        
        // Simulation Parameters
        INFO MC_info;
        int L;
        int N;
        
        // physical information or settings
        double kT;               // System Temperature
        double h1_t, h2_t, h3_t; // unit conversion
        double J1, J2, J3a, J3b;
        // Both set by init_configuration
        vector<double> magnetic_field;

        // == Base class pointer in future == // 
        // Lattice Geometry
        Square latt;
        // Hamiltonian
        Square_ice model;

        // Ising spins, physical varaible
        Sample S0;
        Sample St;
        Sample St_next;

        /* Status variables */ 
        int num_S0_updates;
        int defect_site; // init_site on S0
        int init_site;
        vector<int> other_defects;
        int defect_spin; 
        int agent_site; // current site
        int agent_vertex; // +1: upper block; -1: lower block
        int agent_spin; // current spin
        int csite;

        // Auxillary variables
        vector<int> up_spin;
        vector<int> down_spin;

        /* Variable length queue */
        // Notice: Only FlipOnSite counts
        // Notice: visited map only counts agent flip

        /* LxL map of defect position */
        vector<int> defect_map;
        vector<int> agent_pos_map; 
        vector<int> action_map;
        vector<int> sublatt_map;
        vector<int> difference_map;

        // trajectory: walking trajectory
        vector<int> trajectory_map;
        vector<int> trajectory_counter_map;
        vector<int> trajectory_queue;

        vector<int> action_queue;

        // visited: spin flipping only counts
        vector<int> visited_map;
        vector<int> visited_counter_map;
        vector<int> visited_queue;


        // configuration difference from S0 to St+1
        int difference_counter;
        double defect_density;

        // energy difference counter (?), can be wrong
        double dE;

        // Utilities (Not used now)
	    bool data_analysis;
        Timer TT;
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
        .def("print_INFO", &INFO::print_INFO)
        ;
    to_python_converter<std::vector<int, class std::allocator<int> >, Vec2List<int> >();
    to_python_converter<std::vector<double, class std::allocator<double> >, Vec2List<double> >();
    class_<SquareIceGame>("SquareIceGame", init<INFO>())
        .def("init_configs", &SquareIceGame::InitConfigs)
        .def("init_agent", &SquareIceGame::InitAgent)
        .def("create_defect", &SquareIceGame::CreateDefect)
        .def("start", &SquareIceGame::Start)
        .def("reset", &SquareIceGame::Reset)

        .def("set_magnetic_field", &SquareIceGame::setMagneticField)
        .def("set_coupling", &SquareIceGame::setCoupling)
        .def("set_temperature", &SquareIceGame::setTemperature)

        .def("get_spin", &SquareIceGame::getSpin)
        .def("get_agent_site", &SquareIceGame::getAgentSite)
        .def("get_agent_spin", &SquareIceGame::getAgentSpin)
        .def("get_defect_site", &SquareIceGame::getDefectSite)
        .def("get_defect_spin", &SquareIceGame::getDefectSpin)

        .def("get_agent_map", &SquareIceGame::getAgentPosMap)
        .def("get_action_map", &SquareIceGame::getActionMap)
        .def("get_sublatt_map", &SquareIceGame::getSublattMap)
        .def("get_visited_map", &SquareIceGame::getVisitedMap)
        .def("get_defect_map", &SquareIceGame::getDefectMap)
        .def("get_traject_map", &SquareIceGame::getTrajectoryMap)
        .def("get_traject", &SquareIceGame::getTrajectory)
        .def("get_action_sequence", &SquareIceGame::getActionQueue)

        .def("get_up_spin", &SquareIceGame::getUpSpin)
        .def("get_down_spin", &SquareIceGame::getDownSpin)
        .def("get_S0_defect_density", &SquareIceGame::getS0DefectDensity)
        .def("get_St_defect_density", &SquareIceGame::getStDefectDensity)
        .def("get_S0_states", &SquareIceGame::getS0Configs)
        .def("get_St_states", &SquareIceGame::getStConfigs)
        .def("get_dS_states", &SquareIceGame::getdSConfigs)
        .def("get_St_next_states", &SquareIceGame::getStNextConfigs)
        .def("get_neighbor_sites", &SquareIceGame::getNeighborSites)
        .def("get_neighbor_spins", &SquareIceGame::getNeighborSpins)

        .def("get_diff_ratio", &SquareIceGame::getDiffRatio)
        .def("get_diff_energy", &SquareIceGame::getDiffEnergy)
        .def("get_diff_counter", &SquareIceGame::getDiffCounter)
        .def("get_num_updates", &SquareIceGame::getNumUpdates)

        .def("flip_site", &SquareIceGame::FlipSite)
        .def("flip_on_site", &SquareIceGame::FlipOnSite)
        .def("MCRun", &SquareIceGame::MCRun)

        .def("is_short_closed", &SquareIceGame::isShortClosed)
        .def("is_long_closed", &SquareIceGame::isLongClosed)
        .def("pop_visited_queue", &SquareIceGame::popVisitedQueue)
        .def("pop_trajectory_queue", &SquareIceGame::popTrajectoryQueue)

        .def("metropolis_judgement", &SquareIceGame::MetropolisJudgement)
        .def("act", &SquareIceGame::Act)
        .def("pop_action_queue", &SquareIceGame::popActionQueue)
        .def("pop_shortloop_action_queue", &SquareIceGame::popShortLoopActionQueue)
        ;
}
#endif
