#include "squareice.hpp"
#include <math.h>

SquareIceGame::SquareIceGame(INFO in) : MC_info(in){
    // init model 
    model.init(MC_info);
    // init lattice 
    latt.init(MC_info);

    S0.init(MC_info);
    St.init(MC_info);
    St_next.init(MC_info);

    N = MC_info.Num_sites;
    L = MC_info.lattice_size;

    up_spin.resize(N, 0);
    down_spin.resize(N, 0);

    agent_pos_map.resize(N, 0);
    action_map.resize(N, 0);
    sublatt_map.resize(N, 0);
    defect_map.resize(N, 0);
    visited_map.resize(N, 0);
    visited_counter_map.resize(N, 0);
    difference_map.resize(N, 0);

    trajectory_map.resize(N, 0);
    trajectory_counter_map.resize(N, 0);

    default_reward = 0.0;
    // default setting
    kT = 0.0001;
    h1_t = 0.0;
    h2_t = 0.0;
    h3_t = 0.0;
    J1 = 1.0;
    J2 = 0.0;
    J3a = 0.0;
    J3b = 0.0;
    agent_site = 0;
    agent_vertex = 1;
    init_site = -1;
    agent_spin = 0;
    defect_site = -1;
    defect_spin = 0;
    csite = -1;
    difference_counter = 0;
    dE = 0.0;
    num_S0_updates = 0;

    std::cout << "MCSim initialized.\n";
};

void SquareIceGame::setMagneticField(double h1, double h2, double h3) {
    magnetic_field.resize(3);
	h1_t = h1 * 9.27400968 / 1.3806488;
	h2_t = h2 * 9.27400968 / 1.3806488; // unit conversion
	h3_t = h3 * 9.27400968 / 1.3806488;
};

void SquareIceGame::setCoupling(double J1) {
    model.set_J1(J1);
    J2  = 0.0;
    J3a = 0.0; 
    J3b = 0.0;
};

void SquareIceGame::setTemperature(double T) {
    if (T >= 0.0)
        kT = T;
    std::cout << " ------- Set temperature kT = " << kT << "\n";
}

void SquareIceGame::InitConfigs() {
    // We force only one replica exists
	magnetic_field[0] = h1_t;
	magnetic_field[1] = h2_t;
	magnetic_field[2] = h3_t;

    std::cout << " ------- Set h1_t = " << h1_t 
        << ", h2_t = " << h2_t 
        << ", h3_t = " << h3_t << "\n";

    // Construct Ising replica
	model.initialization(&S0, &latt, 1);
	model.initialize_observable(&S0, &latt, kT, magnetic_field);
    std::cout << " ------- The configuration of each replica is initialized." << std::endl;
}

/*
void SquareIceGame::reset_configs() {}
*/

double SquareIceGame::Start(int site) {
    CreateDefect(site);
    InitAgent(site);
    return getStDefectDensity();
}

int SquareIceGame::CreateDefect(int site) {
    // Create defect on St only, put the agent on it (first)
    // Defect live on S0, and then St should sync to it
    if (site < MC_info.Num_sites && site > 0) {
        defect_site = site;
        // Create defect on S0, flip
        St.Ising[defect_site] *= -1;
        FlipSite(defect_site);
        // Save the defect sign
        defect_spin = getStSpin(defect_site);
        //defect_map[site] = 1;
    } else {
        return -1;
    }
    return defect_spin;
}

int SquareIceGame::RecoverDefect() {
    /*
    if (!other_defects.empty()) {
        for (unsigned int i = 0; i < other_defects.size(); ++i) {
            S0.Ising[other_defects[i]] *= -1;
        }
        other_defects.clear();
    }
    */
    //if (S0.Ising[defect_site] == defect_spin) {
    S0.Ising[defect_site] *= -1;
    //}
    return defect_spin;
}

int SquareIceGame::InitAgent(int site) {
    const int N = MC_info.Num_sites;
    if (site < N) {
        update_agent_site(site);
        init_site = site;
        trajectory_queue.push_back(site);
        trajectory_map[site] = 1;
        trajectory_counter_map[site] += 1;
        visited_queue.push_back(site);
        visited_map[site] = 1;
        visited_counter_map[site] += 1;
   } else {
       return -1;
   }
    return agent_site;
}

void SquareIceGame::go_to_S0() {
    // update spin variables
    St.Ising        = S0.Ising;
    St_next.Ising   = S0.Ising;
    // energies assignment
    St.energy_update       = S0.energy_update;
    St_next.energy_update  = S0.energy_update;
    // magnization not yet
    difference_counter = 0.0;
    dE = 0.0;
    defect_spin = getS0Spin(defect_site);
    agent_spin = getS0Spin(agent_site);
    init_site = -1;

    // search S0
    for (int i = 0; i < N; ++i) {
        int spin = S0.Ising[i];
        switch(spin) {
            case 1:
                up_spin[i] = 1;
            break;
            case -1:
                down_spin[i] = 1;
            break;
        }
    }
    reset_maps();
}

void SquareIceGame::go_to_St() {
    // accpet configurations
    S0.Ising = St.Ising;
    // update energy
    S0.energy_update = model.total_energy(&S0, &latt);
    // update magnetization 
    S0.magnetization_update[0] = model.magnetization(&S0, &latt);
    difference_counter = 0.0;
    num_S0_updates++;
    dE = 0.0;
    reset_maps();
}

/* Update all agent site related information */
void SquareIceGame::update_agent_site(int new_site) {
    int old_agent = agent_site;
    agent_site = new_site;
    agent_spin = getSpin(agent_site);
    // update agent position map
    agent_pos_map[old_agent] = 0;
    agent_pos_map[agent_site] = 1;

    // update subi-lattice map
    for (int i =0; i < 4; ++i) {
        sublatt_map[latt.NN[old_agent][i]] = 0;
        action_map[latt.NN[old_agent][i]] = 0;
    }
    if (1 == latt.sub[old_agent]) {
        sublatt_map[latt.NNN[old_agent][0]] = 0;
        sublatt_map[latt.NNN[old_agent][2]] = 0;
        action_map[latt.NNN[old_agent][0]] = 0;
        action_map[latt.NNN[old_agent][2]] = 0;
    } else {
        sublatt_map[latt.NNN[old_agent][1]] = 0;
        sublatt_map[latt.NNN[old_agent][3]] = 0;
        action_map[latt.NNN[old_agent][1]] = 0;
        action_map[latt.NNN[old_agent][3]] = 0;
    }
    for (int i =0; i < 4; ++i) {
        sublatt_map[latt.NN[agent_site][i]] = 1;
        int spin = getSpin(latt.NN[agent_site][i]);
        if (spin != agent_spin) {
            action_map[latt.NN[agent_site][i]] = 1;
        }
    }
    if (1 == latt.sub[agent_site]) {
        sublatt_map[latt.NNN[agent_site][0]] = 1;
        sublatt_map[latt.NNN[agent_site][2]] = 1;
        if (agent_spin != getSpin(latt.NNN[agent_site][0]))
            action_map[latt.NNN[agent_site][0]] = 1;
        if (agent_spin != getSpin(latt.NNN[agent_site][2]))
            action_map[latt.NNN[agent_site][2]] = 1;
    } else {
        sublatt_map[latt.NNN[agent_site][1]] = 1;
        sublatt_map[latt.NNN[agent_site][3]] = 1;
        if (agent_spin != getSpin(latt.NNN[agent_site][1]))
            action_map[latt.NNN[agent_site][1]] = 1;
        if (agent_spin != getSpin(latt.NNN[agent_site][3]))
            action_map[latt.NNN[agent_site][3]] = 1;
    }

}

void SquareIceGame::reset_maps() {
    // reset defect map
    std::fill(agent_pos_map.begin(), agent_pos_map.end(), 0);
    std::fill(action_map.begin(), action_map.end(), 0);
    std::fill(sublatt_map.begin(), sublatt_map.end(), 0);
    std::fill(defect_map.begin(), defect_map.end(), 0);
    std::fill(visited_map.begin(), visited_map.end(), 0);
    std::fill(visited_counter_map.begin(), visited_counter_map.end(), 0);
    std::fill(trajectory_map.begin(), trajectory_map.end(), 0);
    std::fill(trajectory_counter_map.begin(), trajectory_counter_map.end(), 0);
    trajectory_queue.clear();
    visited_queue.clear();
    action_queue.clear();
}

void SquareIceGame::update_to_newstate() {
    St.Ising = St_next.Ising;
    St.energy_update = St_next.energy_update;
}

void SquareIceGame::MCRun(int mcSteps) {
    TT.timer_begin();
    for (int i = 0; i < mcSteps; ++i) {
            model.MCstep(&S0, &latt);
    }
    std::cout << "Monte carlo runs " << mcSteps << " steps.\n";
    TT.timer_end();
    std::cout << " ------- Simulation time = " << TT.timer_duration() << " seconds. " << std::endl; 
    // Set all states to S0
    go_to_S0();

    std::cout << "Set all states to initial state\n";
}

std::vector<int> SquareIceGame::getdSConfigs() {
    return difference_map;
}

std::vector<int> SquareIceGame::getS0Configs() {
    return S0.get_Ising();
}

std::vector<int> SquareIceGame::getStConfigs() {
    return St.get_Ising();
}

std::vector<int> SquareIceGame::getStNextConfigs() {
    return St_next.get_Ising();
}

int SquareIceGame::getStSpin(int site) {
    return St.get_Ising()[site];
}

int SquareIceGame::getS0Spin(int site) {
    return S0.get_Ising()[site];
}

std::vector<int> SquareIceGame::getNeighborSites(int site) {
    // latt.NN[site][num_of_neighbor];
    std::vector<int> neighbor_sites;
    //neighbor_sites.push_back(site);
    neighbor_sites.push_back(latt.NN[site][0]);
    neighbor_sites.push_back(latt.NN[site][1]);
    neighbor_sites.push_back(latt.NN[site][2]);
    neighbor_sites.push_back(latt.NN[site][3]);
    if (1 == latt.sub[site]) {
        neighbor_sites.push_back(latt.NNN[site][0]);
        neighbor_sites.push_back(latt.NNN[site][2]);
    } else {
        neighbor_sites.push_back(latt.NNN[site][1]);
        neighbor_sites.push_back(latt.NNN[site][3]);
    }
    return neighbor_sites;
}

std::vector<int> SquareIceGame::getNeighborSpins(int site) {
    // latt.NN[site][num_of_neighbor];
    std::vector<int> neighbor_spins;
    //neighbor_spins.push_back(St.Ising[site]);
    neighbor_spins.push_back(St.Ising[latt.NN[site][0]]);
    neighbor_spins.push_back(St.Ising[latt.NN[site][1]]);
    neighbor_spins.push_back(St.Ising[latt.NN[site][2]]);
    neighbor_spins.push_back(St.Ising[latt.NN[site][3]]);
    if (1 == latt.sub[site]) {
        neighbor_spins.push_back(St.Ising[latt.NNN[site][0]]);
        neighbor_spins.push_back(St.Ising[latt.NNN[site][2]]);
        //
    } else {
        neighbor_spins.push_back(St.Ising[latt.NNN[site][1]]);
        neighbor_spins.push_back(St.Ising[latt.NNN[site][3]]);
    }
    return neighbor_spins;
}

std::vector<int> SquareIceGame::getNeighborSitesAndSpins(int site) {
    std::vector<int> concat;
    std::vector<int> sites = getNeighborSites(site);
    std::vector<int> spins = getNeighborSpins(site);
    concat.reserve(sites.size() + spins.size());
    concat.insert(concat.end(), sites.begin(), sites.end());
    concat.insert(concat.end(), spins.begin(), spins.end());
    return concat;
}

double SquareIceGame::MetropolisJudgement() {
  const double N = MC_info.Num_sites;
  double score = 0.0;

  // Metropolis algorithm
  double delta_E = getDiffEnergy();
  double weight = exp(-delta_E / S0.get_temperature() );
  double dice = uni01_sampler();
  int status;

  unsigned int state_change = 0;
  for (int site = 0; site < MC_info.Num_sites; ++site) {
      // calculate all St and S0
    bool same = (St.Ising[site] == S0.Ising[site]);
    if (!same) {
        ++state_change;
        difference_map[site] = 1;
     } else {
        difference_map[site] = 0;
     }
  }

  double state_change_ratio = state_change / double(action_queue.size());

  if (delta_E < 0) {
    go_to_St();
    status = 1;
  } else if (dice < weight) {
    go_to_St();
    // rewards
    status = 0;
  } else {
    status = -1;
  } 

  if (status >=0 && state_change > 1) {
      score = state_change;
      // but too small value
  }

  return score;
} // end of metropolis

void SquareIceGame::Reset() {
    if (getS0DefectDensity() != 0.0) {
        std::cout << "[WORNING]: Defect recover fails\n";
        std::cout << "  defect density = " << getS0DefectDensity() << "\n";
    }
    go_to_S0();
}

double SquareIceGame::Act(int action) {
    double reward = default_reward;
    /* move functions */
    switch (action) {
        case 0: 
            flip_right();
        break;
        case 1: 
            flip_down();
        break;
        case 2: 
            flip_left();
        break;
        case 3: 
            flip_up();
        break;
        case 4: 
            flip_next_up();
        break;
        case 5: 
            flip_next_down();
        break;

        case 6:
            reward = MetropolisJudgement();
        break;

        case 7:
            flip_next_upper_right();
        break;
        case 8:
            reward = MetropolisJudgement();
        break;

        case 9:
            flip_left();
        break;
        case 10:
            flip_up();
        break;
        case 11:
            flip_next_up();
        break;
        case 12:
            flip_next_down();
        break;

        case 13:
            go_to_S0();
        break;

        default:
        break;
    }
    // Add the action to the queue
    action_queue.push_back(action);
    return reward;
}

int SquareIceGame::whichActionToFlipSite(int site) {
    int action = -1;
    if (site == latt.NN[agent_site][0]) {
            action = 7;
    }else if (site == latt.NN[agent_site][1]) {
            action = 8;
    }else if (site == latt.NN[agent_site][2]) {
            action = 9;
    }else if (site == latt.NN[agent_site][3]) {
            action = 10;
    }else if (site == latt.NNN[agent_site][0]) {
            action = 11;
    }else if (site == latt.NNN[agent_site][1]) {
            action = 11;
    }else if (site == latt.NNN[agent_site][2]) {
            action = 12;
    }else if (site == latt.NNN[agent_site][3]) {
            action = 12;
    } else {
            std::cout << "Site is not in agent's neighborhoods\n";
    }
    return action;
}

int SquareIceGame::whichActionToWalkSite(int site) {
    int action = -1;
    if (site == latt.NN[agent_site][0]) {
            action = 0;
    }else if (site == latt.NN[agent_site][1]) {
            action = 1;
    }else if (site == latt.NN[agent_site][2]) {
            action = 2;
    }else if (site == latt.NN[agent_site][3]) {
            action = 3;
    }else if (site == latt.NNN[agent_site][0]) {
            action = 4;
    }else if (site == latt.NNN[agent_site][1]) {
            action = 4;
    }else if (site == latt.NNN[agent_site][2]) {
            action = 5;
    }else if (site == latt.NNN[agent_site][3]) {
            action = 5;
    } else {
            std::cout << "Site is not in agent's neighborhoods\n";
    }
    return action;
}

double SquareIceGame::getDiffEnergy() {
  const double N = MC_info.Num_sites;
  double E_St = model.total_energy(&St, &latt);
  double E_S0 = model.total_energy(&S0, &latt);
  // std::cout << "E_S0 = " << E_S0 << ", E_St = " << E_St << "\n";
  double delta_E = E_St - E_S0;
  return delta_E /= N;
}

double SquareIceGame::calculate_St_defect_density() {
    double defect = 0;
    for (int i = 0; i < MC_info.Num_sites; ++i) {
        if (latt.sub[i] == 1){
            double sum = St.Ising[i] + St.Ising[latt.NN[i][0]] + St.Ising[latt.NN[i][1]] + St.Ising[latt.NNN[i][0]];
            if (sum != 0){
                defect_map[i] = 1;
            } else {
                defect_map[i] = 0;
            }
            defect += abs(sum);
    }
  }
  defect /= double(MC_info.Num_sites);
  //defect *= 2.0;
  //defect = abs(defect);
  return defect;
}

double SquareIceGame::calculate_S0_defect_density() {
    double defect = 0;
    for (int i = 0; i < MC_info.Num_sites; ++i) {
        if (latt.sub[i] == 1){
            double sum = S0.Ising[i] + S0.Ising[latt.NN[i][0]] + S0.Ising[latt.NN[i][1]] + S0.Ising[latt.NNN[i][0]];
            if (sum != 0) {
                if (i != defect_site) {
                    other_defects.push_back(i);
                    // Check the correcetness of the code
                }
            }
            defect += sum;
    }
  }
  defect /= double(MC_info.Num_sites);
  defect *= 2.0;
  return abs(defect);
}

// check on St (for less confusion)
void SquareIceGame::seperate_up_or_down(int site) {
    // get spin from St
    int spin = getSpin(site);
    if (spin == 1) {
        up_spin[site] = 1;
        down_spin[site] = 0;
    } else if(spin == -1) {
        up_spin[site] = 0;
        down_spin[site] = 1;
    }
}

int SquareIceGame::FlipOnSite() {
    ++visited_counter_map[agent_site];
    // add to trajectory
    if (0 == visited_map[agent_site]) {
        visited_map[agent_site] = 1;
    }
    FlipSite(agent_site);
    visited_queue.push_back(agent_site);
    return agent_spin;
}

int SquareIceGame::StayOnSite() {
    ++trajectory_counter_map[agent_site];
    if (0 == trajectory_map[agent_site]) {
        trajectory_map[agent_site] = 1;
    }
    trajectory_queue.push_back(agent_site);

    return agent_spin;
}

int SquareIceGame::getAgentSpin() {
    return getSpin(agent_site);
}

int SquareIceGame::getSpin(int site) {
    return St.Ising[site];
}

// Flip one site on St_next, then update to St
// 
int SquareIceGame::FlipSite(int site) {
    const int N = MC_info.Num_sites;
    // error
    if (site > N) {
        std::cout << "flip site is out of bound!\n";
        return -1;
    }
    // check if the same, true -> counter++ else counter--
    bool same = (S0.Ising[site] == St_next.Ising[site]);
    if (same) {
        ++difference_counter;
    } else {
        --difference_counter;
    }

    //double delta = (-2.0) * model.SE(site, &St_next, &latt);
    St_next.Ising[site] *= -1;
    int spin = St_next.Ising[site];
    //St_next.energy_update += delta;
    //dE += (delta/double(N));

    // update St <- St_next
    update_to_newstate();
    // Update defect map and density
    defect_density = calculate_St_defect_density();
    // then, classifiy the 'on-site' spin to up or down
    seperate_up_or_down(site);
    return spin;
}

// Not follow Ice rule (do not check spin sign)
bool SquareIceGame::go_up() {
    bool isMove = false;
    int newsite;
    if (1 == latt.sub[agent_site]) {
        newsite = latt.NN[agent_site][3];
    } else {
        newsite = latt.NN[agent_site][3];
    }
    if (getSpin(newsite) != agent_spin) isMove = true;
    update_agent_site(newsite);
    StayOnSite();
    return isMove;
}

// Not follow Ice rule (do not check spin sign)
bool SquareIceGame::flip_up() {
    bool isMove = false;
    int newsite;
    if (1 == latt.sub[agent_site]) {
        newsite = latt.NN[agent_site][3];
    } else {
        newsite = latt.NN[agent_site][3];
    }
    if (getSpin(newsite) != agent_spin) isMove = true;
    update_agent_site(newsite);
    StayOnSite();
    FlipOnSite();
    return isMove;
}

// Not follow Ice rule (do not check spin sign)
bool SquareIceGame::go_down() {
    bool isMove = false;
    int newsite;
    if (1 == latt.sub[agent_site]) {
        newsite = latt.NN[agent_site][1];
    } else {
        newsite = latt.NN[agent_site][1];
    }
    if (getSpin(newsite) != agent_spin) isMove = true;
    update_agent_site(newsite);
    StayOnSite();
    return isMove;
}

// Not follow Ice rule (do not check spin sign)
bool SquareIceGame::flip_down() {
    bool isMove = false;
    int newsite;
    if (1 == latt.sub[agent_site]) {
        newsite = latt.NN[agent_site][1];
    } else {
        newsite = latt.NN[agent_site][1];
    }
    if (getSpin(newsite) != agent_spin) isMove = true;
    update_agent_site(newsite);
    StayOnSite();
    FlipOnSite();
    return isMove;
}

// Not follow Ice rule (do not check spin sign)
bool SquareIceGame::go_right() {
    bool isMove = false;
    int newsite;
    if (1 == latt.sub[agent_site]) {
        newsite = latt.NN[agent_site][0];
    } else {
        newsite = latt.NN[agent_site][0];
    }
    if (getSpin(newsite) != agent_spin) isMove = true;
    update_agent_site(newsite);
    StayOnSite();
    return isMove;
}

// Not follow Ice rule (do not check spin sign)
bool SquareIceGame::flip_right() {
    bool isMove = false;
    int newsite;
    if (1 == latt.sub[agent_site]) {
        newsite = latt.NN[agent_site][0];
    } else {
        newsite = latt.NN[agent_site][0];
    }
    if (getSpin(newsite) != agent_spin) isMove = true;
    update_agent_site(newsite);
    StayOnSite();
    FlipOnSite();
    return isMove;
}

// Not follow Ice rule (do not check spin sign)
bool SquareIceGame::go_left() {
    bool isMove = false;
    int newsite;
    if (1 == latt.sub[agent_site]) {
        newsite = latt.NN[agent_site][2];
    } else {
        newsite = latt.NN[agent_site][2];
    }
    if (getSpin(newsite) != agent_spin) isMove = true;
    update_agent_site(newsite);
    StayOnSite();
    return isMove;
}

// Not follow Ice rule (do not check spin sign)
bool SquareIceGame::flip_left() {
    bool isMove = false;
    int newsite;
    if (1 == latt.sub[agent_site]) {
        newsite = latt.NN[agent_site][2];
    } else {
        newsite = latt.NN[agent_site][2];
    }
    if (getSpin(newsite) != agent_spin) isMove = true;
    update_agent_site(newsite);
    StayOnSite();
    FlipOnSite();
    return isMove;
}

// Not follow Ice rule (do not check spin sign)
bool SquareIceGame::go_next_up() {
    bool isMove = false;
    int newsite;
    if (1 == latt.sub[agent_site]) {
        newsite = latt.NNN[agent_site][2];
    } else {
        newsite = latt.NNN[agent_site][3];
    }
    if (getSpin(newsite) != agent_spin) isMove = true;
    update_agent_site(newsite);
    StayOnSite();
    return isMove;
}

// Not follow Ice rule (do not check spin sign)
bool SquareIceGame::flip_next_up() {
    bool isMove = false;
    int newsite;
    if (1 == latt.sub[agent_site]) {
        newsite = latt.NNN[agent_site][2];
    } else {
        newsite = latt.NNN[agent_site][3];
    }
    if (getSpin(newsite) != agent_spin) isMove = true;
    update_agent_site(newsite);
    StayOnSite();
    FlipOnSite();
    return isMove;
}

// Not follow Ice rule (do not check spin sign)
bool SquareIceGame::go_next_down() {
    bool isMove = false;
    int newsite;
    if (1 == latt.sub[agent_site]) {
        newsite = latt.NNN[agent_site][0];
    } else {
        newsite = latt.NNN[agent_site][1];
    }
    if (getSpin(newsite) != agent_spin) isMove = true;
    update_agent_site(newsite);
    StayOnSite();
    return isMove;
}

// Not follow Ice rule (do not check spin sign)
bool SquareIceGame::flip_next_down() {
    bool isMove = false;
    int newsite;
    if (1 == latt.sub[agent_site]) {
        newsite = latt.NNN[agent_site][0];
    } else {
        newsite = latt.NNN[agent_site][1];
    }
    if (getSpin(newsite) != agent_spin) isMove = true;
    update_agent_site(newsite);
    StayOnSite();
    FlipOnSite();
    return isMove;
}

bool SquareIceGame::flip_next_upper_right() {
    bool isMove = false;
    int newsite = latt.NNN[agent_site][3];
    if (getSpin(newsite) != agent_spin) isMove = true;
    update_agent_site(newsite);
    StayOnSite();
    FlipOnSite();
    return isMove;
}

bool SquareIceGame::flip_next_upper_left() {
    bool isMove = false;
    int newsite = latt.NNN[agent_site][2];
    if (getSpin(newsite) != agent_spin) isMove = true;
    update_agent_site(newsite);
    StayOnSite();
    FlipOnSite();
    return isMove;
}

bool SquareIceGame::flip_next_lower_left() {
    bool isMove = false;
    int newsite = latt.NNN[agent_site][1];
    if (getSpin(newsite) != agent_spin) isMove = true;
    update_agent_site(newsite);
    StayOnSite();
    FlipOnSite();
    return isMove;
}

bool SquareIceGame::flip_next_lower_right() {
    bool isMove = false;
    int newsite = latt.NNN[agent_site][0];
    if (getSpin(newsite) != agent_spin) isMove = true;
    update_agent_site(newsite);
    StayOnSite();
    FlipOnSite();
    return isMove;
}

bool SquareIceGame::isShortClosed() {
    bool is_closed = false;
    if (std::find(visited_counter_map.begin(), visited_counter_map.end(), 2) != visited_counter_map.end()) {
        is_closed = true;
        csite = std::find(visited_counter_map.begin(), visited_counter_map.end(), 2) - visited_counter_map.begin();
        std::cout << "Crossed site is " << csite << "\n";
    }
    return is_closed;
}

bool SquareIceGame::isLongClosed() {
    bool is_closed = false;
    int head = visited_queue[0];
    if (head >= 0 && head < MC_info.Num_sites) {
        if (visited_counter_map[head] >= 2 ) {
            is_closed = true;
        }
    }
    return is_closed;
}

int SquareIceGame::popVisitedQueue() {
    if (!visited_queue.empty()) {
        int site = visited_queue.back();
        visited_queue.pop_back();
        return site;
    } else {
        return -1;
    }
}

int SquareIceGame::popTrajectoryQueue() {
    if (!trajectory_queue.empty()) {
        int site = trajectory_queue.back();
        trajectory_queue.pop_back();
        return site;
    } else {
        return -1;
    }
}

int SquareIceGame::popActionQueue() {
    if (!action_queue.empty()) {
        int act = action_queue.back();
        action_queue.pop_back();
        return act;
    } else {
        return -1;
    }
}

int SquareIceGame::dual_action_flip(int action) {
    int dual_action = -1;
    switch (action) {
        case 7: 
            dual_action = 9;
        break;
        case 8: 
            dual_action = 10;
        break;
        case 9: 
            dual_action = 7;
        break;
        case 10: 
            dual_action = 8;
        break;
        case 11: 
            dual_action = 12;
        break;
        case 12: 
            dual_action = 11;
        break;
        default:
            dual_action = action;
        break;
    }
    return dual_action;
}

int SquareIceGame::dual_action_walk(int action) {
    int dual_action = -1;
    switch (action) {
        case 7: 
            dual_action = 2;
        break;
        case 8: 
            dual_action = 3;
        break;
        case 9: 
            dual_action = 0;
        break;
        case 10: 
            dual_action = 1;
        break;
        case 11: 
            dual_action = 5;
        break;
        case 12: 
            dual_action = 4;
        break;
    }
    return dual_action;
}

bool SquareIceGame::ShortLoopUpdate() {
    bool is_updated = false;
    if(isShortClosed()) {
        std::cout << "Short loop is closed!\n";
        // crossed site
        int csite = std::find(visited_counter_map.begin(), visited_counter_map.end(), 2) - visited_counter_map.begin();
        // flip back those not in loop 
        while(!visited_queue.empty()) {
            int site = visited_queue.back();
            visited_queue.pop_back();
            std::cout << "pop out " << site << " sites\n";
            if (visited_map[csite] == 0) {
                FlipSite(site);
                std::cout << "flip back " << site << " site's spin\n";
            }
            visited_map[site]--;
        }
        is_updated = true;
    }
    return is_updated;
}

int SquareIceGame::popShortLoopActionQueue() {
    int action = -1;
    //int dual_action = -1;
    if (!trajectory_queue.empty()) {
        //action = popActionQueue();
        trajectory_queue.pop_back();
        int prev_site = trajectory_queue.back();
        trajectory_queue.pop_back();
        int site = agent_site;
        std::cout << "\tprevious site = " << prev_site << " and current site = " << site << "\n";
        std::cout << "\tvisited[csite] = " << visited_counter_map[csite] << "\n";
        if (visited_counter_map[csite] == 0) {
            //dual_action = dual_action_flip(action);
            action = whichActionToFlipSite(prev_site);
            std::cout << "\tFlip back " << site << " site's spin\n";
        } else {
            //dual_action = dual_action_walk(action);
            action = whichActionToWalkSite(prev_site);
            std::cout << "\tWalk through " << site << " site's spin\n";
        }
        visited_counter_map[site]--;
    }
    return action;
}
