#include "sqice.hpp"
#include <math.h>

SQIceGame::SQIceGame (INFO info) : sim_info(info) {
    model.init(sim_info);
    latt.init(sim_info);
    ice_config.init (sim_info);

    // pre-determined parameters
    N = sim_info.Num_sites;
    L = sim_info.lattice_size;

    h1_t = 0.0;
    h2_t = 0.0;
    h3_t = 0.0;
    J1 = 1.0;
    agent_site = 0;

    mag_fields.push_back(h1_t);
    mag_fields.push_back(h2_t);
    mag_fields.push_back(h3_t);

    // initialze all the member data
    state_0.resize(N, 0);
    state_t.resize(N, 0);
    state_tp1.resize(N, 0);
    

    std::cout << "[GAME] Square Ice Game is created.\n";
}

void SQIceGame::InitModel() {
    model.set_J1(J1);
	model.initialization(&ice_config, &latt, 1);
	model.initialize_observable(&ice_config, &latt, kT, mag_fields);
    std::cout << "[GAME] All physical parameters are initialized\n";
}

void SQIceGame::SetTemperature(double T) {
    if (T >= 0.0)
        kT = T;
    std::cout << "[GAME] Set temperature kT = " << kT << "\n";
}

void SQIceGame::MCRun(int mcSteps) {
    tt.timer_begin();
    for (int i = 0; i < mcSteps; ++i) {
        model.MCstep(&ice_config, &latt);
    }
    tt.timer_end();
    std::cout << "[GAME] Monte Carlo runs " 
              << mcSteps << " steps with "
              << tt.timer_duration() << " seconds.\n"; 
    
    // Check whether it is icestates
    double Etot = model.total_energy(&ice_config, &latt);
    std::cout << "[GAME] Total Energy E = " << Etot << "\n";

    // Get the ising variables
    state_0 = ice_config.Ising;
    state_t = ice_config.Ising;
    std::cout << "[GAME] Average Energy E = " << _cal_energy_of_state(state_0) << "\n";
    std::cout << "[GAME] Defect Density D = " << _cal_defect_density_of_state(state_0) << "\n";
}

// member functions 

void SQIceGame::set_agent_site(int site) {
    if (site >= 0 && site < N) {
        agent_site = site;
    } else {
        std::cout << "[GAME] WORNING, Set Agent on Illegal Site!\n";
    }
}

int SQIceGame::get_agent_site() {
    return agent_site;
}

vector<int> SQIceGame::neighbor_sites() {
    vector<int> neighbors(6);
    neighbors[0] = latt.NN[agent_site][0];
    neighbors[1] = latt.NN[agent_site][1];
    neighbors[2] = latt.NN[agent_site][2];
    neighbors[3] = latt.NN[agent_site][3];
    if (latt.sub[agent_site] == 1) {
        neighbors[4] = latt.NNN[agent_site][0];
        neighbors[5] = latt.NNN[agent_site][2];
    } else {
        neighbors[4] = latt.NNN[agent_site][1];
        neighbors[5] = latt.NNN[agent_site][3];
    }
    return neighbors;
}

vector<int> SQIceGame::neighbor_spins() {
    vector<int> neighbors(6);
    neighbors[0] = state_t[latt.NN[agent_site][0]];
    neighbors[1] = state_t[latt.NN[agent_site][1]];
    neighbors[2] = state_t[latt.NN[agent_site][2]];
    neighbors[3] = state_t[latt.NN[agent_site][3]];
    if (latt.sub[agent_site] == 1) {
        neighbors[4] = state_t[latt.NNN[agent_site][0]];
        neighbors[5] = state_t[latt.NNN[agent_site][2]];
    } else {
        neighbors[4] = state_t[latt.NNN[agent_site][1]];
        neighbors[5] = state_t[latt.NNN[agent_site][3]];
    }
    return neighbors;
}

int SQIceGame::get_neighbor_site_by_direction(DIR dir) {
    int site = agent_site;
    if (dir == RIGHT) {
            site = latt.NN[site][0];
    } else if (dir == DOWN) {
            site = latt.NN[site][1];
    } else if (dir == LEFT) {
            site = latt.NN[site][2];
    } else if (dir == UP) { 
            site = latt.NN[site][3];
    } else if (dir == LOWER_RIGHT) {
            site = latt.NNN[site][0];
    } else if (dir == LOWER_LEFT) {
            site = latt.NN[site][1];
    } else if (dir == UPPER_LEFT) {
            site = latt.NN[site][2];
    } else if (dir == UPPER_RIGHT) {
            site = latt.NN[site][3];
    }
    return site;
}

void SQIceGame::TEST(){
    std::cout << "/=========================================/\n";

    set_agent_site(5000);
    std::cout << "Set Agent site at 5000, agent = " << get_agent_site() << "\n";
    set_agent_site(100);
    std::cout << "Set Agent site at 100, agent = " << get_agent_site() << "\n";
    std::cout << "Its neighbors are ";
    _print_vector(neighbor_sites());
    std::cout << ".... with corresponding spins are ";
    _print_vector(neighbor_spins());

    std::cout << "whose up site is " << get_neighbor_site_by_direction(UP) << "\n";
    std::cout << "whose down site is " << get_neighbor_site_by_direction(DOWN) << "\n";
    std::cout << "whose left site is " << get_neighbor_site_by_direction(LEFT) << "\n";
    std::cout << "whose right site is " << get_neighbor_site_by_direction(RIGHT) << "\n";

    std::cout << "/=========================================/\n";
}

// ###### private function ######
double SQIceGame::_cal_energy_of_state(vector<int> & s) {
    double eng = 0.0;
    for (uint i = 0; i < N; i++) {
        double se;
        if (latt.sub[i] == 1) {
            se = s[i] * (s[latt.NN[i][0]]
                        +s[latt.NN[i][1]]
                        +s[latt.NN[i][2]]
                        +s[latt.NN[i][3]]
                        +s[latt.NNN[i][0]]
                        +s[latt.NNN[i][2]]
                        );
        } else {
            se = s[i] * (s[latt.NN[i][0]]
                        +s[latt.NN[i][1]]
                        +s[latt.NN[i][2]]
                        +s[latt.NN[i][3]]
                        +s[latt.NNN[i][1]]
                        +s[latt.NNN[i][3]]
                        );
        }
        eng += J1 * se;
    }
    eng /= 2.0;
    eng /= N;
    return eng;
}

double SQIceGame::_cal_defect_density_of_state(vector<int> & s) {
    double dd = 0.0;
    for (uint i = 0; i < N; ++i) {
        if (latt.sub[i] == 1) {
            dd += s[i] + s[latt.NN[i][0]] + s[latt.NN[i][1]] + s[latt.NNN[i][0]];
        }
    }
    dd /= 2.0;
    dd /= N;
    return dd;
}

void SQIceGame::_print_vector(vector<int> v) {
    std::cout << "[";
    for (const auto &i: v)
        std::cout << i << ", ";
    std::cout << "]\n";
}