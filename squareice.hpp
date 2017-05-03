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

        // state transition
        void go_to_S0();            // S0 <- St
        void go_to_St();            // S0 -> St
        void update_to_newstate();  // St -> St+1

        int FlipSite(int site);
        int FlipOnSite();

        // actions mapping
        // a[0]: UP 
        // a[1]: DOWN
        // a[2]: LEFT
        // a[3]: RIGHT
        // a[4]: NEXTUP
        // a[5]: NEXTDOWN
        // ... retrun isIceRule
        //
        // Many different action element
        bool ActDirection(int action);
        bool ActDirectionIce(int action);
        bool ActSpinSign(int action);
        bool ActOneWay(int action);
        bool ActFlipSite(int site, int action);

        // <S0, St>
        std::vector<double> MetropolisJudgement();

        bool ShortLoopUpdate();
        bool LongLoopUpdate();

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
        inline int getNumShortLoopOccurs() {return num_shortloop_occurs;};

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

        /* */

        /* Observation Functions */
        std::vector<int> getNeighborSites(int site);
        std::vector<int> getNeighborSpins(int site);

        // combined methods
        std::vector<int> getNeighborSitesAndSpins(int site);
        std::vector<int> getNeighborSitesWithSpins(int site);

        std::vector<double> getNeighborNormalizedSites(int site);
        std::vector<double> getNeighborNormalizedSitesAndSpins(int site);
        std::vector<double> getNeighborNormalizedSitesWithSpins(int site);

        // lazy, generate configuration-difference only when called
        std::vector<int> get_diff_configs();

        void set_virt_configs(boost::python::list &ns);

        // Utilities functions
        void show_info() {MC_info.print_INFO();};
        inline int num_of_directional_actions() {return 4+1;};
        inline int num_of_directional_actions_icerule() {return 6+1;};
        inline int num_of_flip_site_actions() {return 2;};

    private:
        // Private Functions // 
        int inline PBD (int site, int d) {return ((site+d) % L + L) % L;}
        
        // flip the site where agent is on, return spin after flipped
        
        /* Random Movement */
        // move and flip the next site
        bool go_spinup();
        bool go_spindown();

        /* Directional Movement */
        // move the direction
        bool go_nn_right();
        bool go_nn_left();
        bool go_nn_up();
        bool go_nn_down();
        // The following two functions depends on-site spin
        bool go_nnn_up();
        bool go_nnn_down();

        // move the direction
        bool go_right();
        bool go_left();
        bool go_up();
        bool go_down();

        /* Same Spin Sign Movement */
        // satisfy ice rule
        bool go_same_spin_up();
        bool go_same_spin_down();
        bool go_same_spin_left();
        bool go_same_spin_right();

        void update_agent_site(int new_site);
        void seperate_up_or_down(int site);

        double calculate_S0_defect_density();
        double calculate_St_defect_density();

        void reset_maps();

        // Private Data //
        
        // Simulation Parameters
        INFO MC_info;
        int L;
        int N;
        
        static const int num_reward_types = 4;

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
        int agent_spin; // current spin

        // Auxillary variables
        vector<int> up_spin;
        vector<int> down_spin;

        /* Variable length queue */
        // Notice: Only FlipOnSite counts
        vector<int> walked_sites;
        // Notice: visited map only counts agent flip
        vector<int> visited_map;

        /* LxL map of defect position */
        vector<int> defect_map;
        vector<int> agent_pos_map; 
        vector<int> action_map;
        vector<int> sublatt_map;
        vector<int> difference_map;
        vector<int> trajectory_map;

        bool is_short_closed();
        bool is_long_closed();
        int num_shortloop_occurs; 
        int num_longloop_occurs;

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
        .def("get_neighbor_sites_and_spins", &SquareIceGame::getNeighborSitesAndSpins)
        .def("get_neighbor_sites_with_spins", &SquareIceGame::getNeighborSitesWithSpins)
        .def("get_neighbor_normalized_sites", &SquareIceGame::getNeighborNormalizedSites)
        .def("get_neighbor_normalized_sites_and_spins", &SquareIceGame::getNeighborNormalizedSitesAndSpins)
        .def("get_neighbor_normalized_sites_with_spins", &SquareIceGame::getNeighborNormalizedSitesWithSpins)

        .def("get_diff_ratio", &SquareIceGame::getDiffRatio)
        .def("get_diff_energy", &SquareIceGame::getDiffEnergy)
        .def("get_diff_counter", &SquareIceGame::getDiffCounter)
        .def("get_num_updates", &SquareIceGame::getNumUpdates)
        .def("get_num_shortloop_occurs", &SquareIceGame::getNumShortLoopOccurs)

        .def("flip_site", &SquareIceGame::FlipSite)
        .def("flip_on_site", &SquareIceGame::FlipOnSite)
        .def("MCRun", &SquareIceGame::MCRun)
        .def("short_loop_update", &SquareIceGame::ShortLoopUpdate)

        .def("metropolis_judgement", &SquareIceGame::MetropolisJudgement)
        .def("act_direction", &SquareIceGame::ActDirection)
        .def("act_direction_ice", &SquareIceGame::ActDirectionIce)
        .def("act_spinsign", &SquareIceGame::ActSpinSign)
        .def("act_flip_site", &SquareIceGame::ActFlipSite)
        ;
}
#endif
