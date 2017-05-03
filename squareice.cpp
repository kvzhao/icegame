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
    difference_map.resize(N, 0);
    trajectory_map.resize(N, 0);

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
    init_site = -1;
    agent_spin = 0;
    defect_site = -1;
    defect_spin = 0;
    difference_counter = 0;
    dE = 0.0;
    num_S0_updates = 0;
    num_shortloop_occurs = 0;

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

int SquareIceGame::CreateDefect(int site) {
    // Create defect on St only, put the agent on it (first)
    // Defect live on S0, and then St should sync to it
    if (site < MC_info.Num_sites && site > 0) {
        defect_site = site;
        // Create defect on S0, flip
        // St.Ising[defect_site] *= -1;
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
    std::fill(trajectory_map.begin(), trajectory_map.end(), 0);
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
    neighbor_sites.push_back(site);
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

std::vector<double> SquareIceGame::getNeighborNormalizedSites(int site) {
    // latt.NN[site][num_of_neighbor];
    double denum = 1.0 / double(MC_info.Num_sites);
    std::vector<double> neighbor_sites;
    neighbor_sites.push_back(site * denum);
    neighbor_sites.push_back(latt.NN[site][0] * denum);
    neighbor_sites.push_back(latt.NN[site][1] * denum);
    neighbor_sites.push_back(latt.NN[site][2] * denum);
    neighbor_sites.push_back(latt.NN[site][3] * denum);
    if (1 == latt.sub[site]) {
        neighbor_sites.push_back(latt.NNN[site][0] * denum);
        neighbor_sites.push_back(latt.NNN[site][2] * denum);
    } else {
        neighbor_sites.push_back(latt.NNN[site][1] * denum);
        neighbor_sites.push_back(latt.NNN[site][3] * denum);
    }
    return neighbor_sites;
}

std::vector<int> SquareIceGame::getNeighborSpins(int site) {
    // latt.NN[site][num_of_neighbor];
    std::vector<int> neighbor_spins;
    neighbor_spins.push_back(St.Ising[site]);
    neighbor_spins.push_back(St.Ising[latt.NN[site][0]]);
    neighbor_spins.push_back(St.Ising[latt.NN[site][1]]);
    neighbor_spins.push_back(St.Ising[latt.NN[site][2]]);
    neighbor_spins.push_back(St.Ising[latt.NN[site][3]]);
    if (1 == latt.sub[site]) {
        neighbor_spins.push_back(St.Ising[latt.NNN[site][0]]);
        neighbor_spins.push_back(St.Ising[latt.NNN[site][2]]);
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

std::vector<double> SquareIceGame::getNeighborNormalizedSitesAndSpins(int site) {
    std::vector<double> concat;
    std::vector<double> sites = getNeighborNormalizedSites(site);
    std::vector<int> intspin = getNeighborSpins(site);
    std::vector<double> spins(intspin.begin(), intspin.end());
    concat.reserve(sites.size() + spins.size());
    concat.insert(concat.end(), sites.begin(), sites.end());
    concat.insert(concat.end(), spins.begin(), spins.end());
    return concat;
}

std::vector<int> SquareIceGame::getNeighborSitesWithSpins(int site) {
    std::vector<int> concat;
    concat.push_back(site);
    concat.push_back(St.Ising[site]);
    concat.push_back(latt.NN[site][0]);
    concat.push_back(St.Ising[latt.NN[site][0]]);
    concat.push_back(latt.NN[site][1]);
    concat.push_back(St.Ising[latt.NN[site][1]]);
    concat.push_back(latt.NN[site][2]);
    concat.push_back(St.Ising[latt.NN[site][2]]);
    concat.push_back(latt.NN[site][3]);
    concat.push_back(St.Ising[latt.NN[site][3]]);
    if (1 == latt.sub[site]) {
        concat.push_back(latt.NNN[site][0]);
        concat.push_back(St.Ising[latt.NNN[site][0]]);
        concat.push_back(latt.NNN[site][2]);
        concat.push_back(St.Ising[latt.NNN[site][2]]);
    } else {
        concat.push_back(latt.NNN[site][1]);
        concat.push_back(St.Ising[latt.NNN[site][1]]);
        concat.push_back(latt.NNN[site][3]);
        concat.push_back(St.Ising[latt.NNN[site][3]]);
    }
    return concat;
}

std::vector<double> SquareIceGame::getNeighborNormalizedSitesWithSpins(int site) {
    double denum = 1.0 / double(MC_info.Num_sites);
    std::vector<double> concat;
    concat.push_back(site * denum);
    concat.push_back(double(St.Ising[site]));
    concat.push_back(latt.NN[site][0] * denum);
    concat.push_back(double(St.Ising[latt.NN[site][0]]));
    concat.push_back(latt.NN[site][1] * denum);
    concat.push_back(double(St.Ising[latt.NN[site][1]]));
    concat.push_back(latt.NN[site][2] * denum);
    concat.push_back(double(St.Ising[latt.NN[site][2]]));
    concat.push_back(latt.NN[site][3] * denum);
    concat.push_back(St.Ising[latt.NN[site][3]]);
    if (1 == latt.sub[site]) {
        concat.push_back(latt.NNN[site][0] * denum);
        concat.push_back(double(St.Ising[latt.NNN[site][0]]));
        concat.push_back(latt.NNN[site][2] * denum);
        concat.push_back(double(St.Ising[latt.NNN[site][2]]));
    } else {
        concat.push_back(latt.NNN[site][1] * denum);
        concat.push_back(double(St.Ising[latt.NNN[site][1]]));
        concat.push_back(latt.NNN[site][3] * denum);
        concat.push_back(double(St.Ising[latt.NNN[site][3]]));
    }
    return concat;
}

std::vector<double> SquareIceGame::MetropolisJudgement() {
  std::vector<double> reward;
  reward.resize(num_reward_types);
  const double N = MC_info.Num_sites;

  // Metropolis algorithm
  double delta_E = getDiffEnergy();
  double weight = exp(-delta_E / S0.get_temperature() );
  double dice = uni01_sampler();

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

  reward[1] = state_change / N;
  reward[2] = delta_E;
  reward[3] = weight;

  if (delta_E < 0) {
    go_to_St();
    reward[0] = 1.0;
  } else if (dice < weight) {
    go_to_St();
    // rewards
    reward[0] = 0.0;
  } else {
    reward[0] = -1.0;
  } 
  return reward;
} // end of metropolis

void SquareIceGame::Reset() {
    if (getS0DefectDensity() != 0.0) {
        std::cout << "[WORNING]: Defect recover fails\n";
        std::cout << "  defect density = " << getS0DefectDensity() << "\n";
    }
    go_to_S0();
}

bool SquareIceGame::ActDirection(int action) {
    bool isIceRule = false;
    /* move functions */
    switch (action) {
        case 0: 
            isIceRule = go_right();
        break;
        case 1: 
            isIceRule = go_down();
        break;
        case 2: 
            isIceRule = go_left();
        break;
        case 3: 
            isIceRule = go_up();
        break;
        default:
        break;
    }
    /* evaluate this step */
    // rewards vector + isMove 
    // r[0]: acceptance status
    // r[1]: fraction of difference
    // r[2]: weights of boltzman factor
    // r[3]: delta energy
    // r[4]: 1 ice rule move
    //
    return isIceRule;
}

bool SquareIceGame::ActDirectionIce(int action) {
    bool isIceRule = false;
    /* move functions */
    switch (action) {
        case 0: 
            isIceRule = go_nn_right();
        break;
        case 1: 
            isIceRule = go_nn_down();
        break;
        case 2: 
            isIceRule = go_nn_left();
        break;
        case 3: 
            isIceRule = go_nn_up();
        break;
        case 4: 
            isIceRule = go_nnn_up();
        break;
        case 5:
            isIceRule = go_nnn_down();
        break;
        default:
        break;
    }
    /* evaluate this step */
    // rewards vector + isMove 
    // r[0]: acceptance status
    // r[1]: fraction of difference
    // r[2]: weights of boltzman factor
    // r[3]: delta energy
    // r[4]: 1 ice rule move
    //
    return isIceRule;
}

bool SquareIceGame::ActSpinSign(int action) {
    bool isIceRule = false;
    switch (action) {
        case 0 :
            isIceRule = go_same_spin_up();
            break;
        case 1:
            isIceRule = go_same_spin_down();
            break;
        default:
        break;
    }
    return isIceRule;
}

bool SquareIceGame::ActOneWay(int action) {
    bool isIceRule = false;
    /* move functions */
    switch (action) {
        case 0: 
            isIceRule = go_nn_up();
        break;
        case 1: 
            isIceRule = go_nn_right();
        break;
        case 2: 
        // This action maybe crucial
            isIceRule = go_nnn_down();
        break;
        default:
        break;
    }
    /* evaluate this step */
    // rewards vector + isMove 
    // r[0]: acceptance status
    // r[1]: fraction of difference
    // r[2]: weights of boltzman factor
    // r[3]: delta energy
    // r[4]: 1 ice rule move
    //
    return isIceRule;
}

bool SquareIceGame::ActFlipSite(int site, int action) {
    bool isFlip = false;
    //if (action == 0) return isFlip;
    //int onsite_spin = getSpin(site);
    update_agent_site(site);
    if (action == 1) {
        // Flip the site
        FlipOnSite();
        isFlip = true;
    }
    return isFlip;
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
    FlipSite(agent_site);
    walked_sites.push_back(agent_site);
    ++visited_map[agent_site];
    // add to trajectory
    if (0 == trajectory_map[agent_site]) {
        trajectory_map[agent_site] = 1;
    }
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

// Performance issues
// seems legacy code
void SquareIceGame::set_virt_configs(boost::python::list &ns) {
    vector<int> newising;
    for (int i=0; i < len(ns); ++i) {
        newising.push_back(boost::python::extract<int>(ns[i]));
    }
    St.Ising = newising;
}

// Not follow Ice rule (do not check spin sign)
bool SquareIceGame::go_right() {
    bool isMove = false;
    int x = floor(agent_site % L);
    int y = floor(agent_site / L);
    x = PBD(x, 1);
    int site = int(y * L + x);
    int new_spin = getSpin(site);
    if (new_spin == getAgentSpin()) {
        isMove = true;
    }
    update_agent_site(site);
    FlipOnSite();
    return isMove;
}

// Not follow Ice rule (do not check spin sign)
bool SquareIceGame::go_down() {
    bool isMove = false;
    int x = floor(agent_site % L);
    int y = floor(agent_site / L);
    y = PBD(y, 1);
    int site = int(y * L + x);
    int new_spin = getSpin(site);
    if (new_spin == getAgentSpin()) {
        isMove = true;
    }
    update_agent_site(site);
    FlipOnSite();
    return isMove;
}

// Not follow Ice rule (do not check spin sign)
bool SquareIceGame::go_up() {
    bool isMove = false;
    int x = floor(agent_site % L);
    int y = floor(agent_site / L);
    y = PBD(y, -1);
    int site = int(y * L + x);
    int new_spin = getSpin(site);
    if (new_spin == getAgentSpin()) {
        isMove = true;
    }
    update_agent_site(site);
    FlipOnSite();
    return isMove;
}

// Not follow Ice rule (do not check spin sign)
bool SquareIceGame::go_left() {
    bool isMove = false;
    int x = floor(agent_site % L);
    int y = floor(agent_site / L);
    x = PBD(x, -1);
    int site = int(y * L + x);
    int new_spin = getSpin(site);
    if (new_spin == getAgentSpin()) {
        isMove = true;
    }
    update_agent_site(site);
    FlipOnSite();
    return isMove;
}

/*bool SquareIceGame::go_nn_left() {
    bool isMove = false;
    int x = floor(agent_site % L);
    int y = floor(agent_site / L);
    x = PBD(x, -1);
    int site = int(y * L + x);
    int new_spin = getSpin(site);
    if (new_spin == getAgentSpin()) {
        isMove = true;
    }
    if(isMove) {
        update_agent_site(site);
        FlipOnSite();
    }
    return isMove;
}*/

/*bool SquareIceGame::go_nn_right() {
    bool isMove = false;
    int x = floor(agent_site % L);
    int y = floor(agent_site / L);
    x = PBD(x, 1);
    int site = int(y * L + x);
    int new_spin = getSpin(site);
    if (new_spin == getAgentSpin()) {
        isMove = true;
    }
    if (isMove) {
        update_agent_site(site);
        FlipOnSite();
    }
    return isMove;
}*/

bool SquareIceGame::go_nnn_down() {
    bool isMove = false;
    int site; 
    if (1 == latt.sub[agent_site]) {
        site = latt.NNN[agent_site][2];
    } else {
        site = latt.NNN[agent_site][3];
    }
    int new_spin = getSpin(site);
    if (new_spin == getAgentSpin()) {
        isMove = true;
    }
    if (isMove) {
        update_agent_site(site);
        FlipOnSite();
    }
    return isMove;
}

bool SquareIceGame::go_nnn_up() {
    bool isMove = false;
    int site; 
    if (1 == latt.sub[agent_site]) {
        site = latt.NNN[agent_site][0];
    } else {
        site = latt.NNN[agent_site][1];
    }
    int new_spin = getSpin(site);
    if (new_spin == getAgentSpin()) {
        isMove = true;
    }
    if (isMove) {
        update_agent_site(site);
        FlipOnSite();
    }
    return isMove;
}

/*bool SquareIceGame::go_nn_down() {
    bool isMove = false;
    int x = floor(agent_site % L);
    int y = floor(agent_site / L);
    y = PBD(y, 1);
    int site = int(y * L + x);
    int new_spin = getSpin(site);
    if (new_spin == getAgentSpin()) {
        isMove = true;
    }
    if (isMove) {
        update_agent_site(site);
        FlipOnSite();
    }
    return isMove;
} */

/*bool SquareIceGame::go_nn_up() {
    bool isMove = false;
    int x = floor(agent_site % L);
    int y = floor(agent_site / L);
    y = PBD(y, -1);
    int site = int(y * L + x);
    int new_spin = getSpin(site);
    if (new_spin == getAgentSpin()) {
        isMove = true;
    }
    if (isMove) {
        update_agent_site(site);
        FlipOnSite();
    }
    return isMove;
}*/

bool SquareIceGame::go_nn_up() {
    bool isMove = false;
    int up_site = latt.NN[agent_site][3]; // up site
    int new_spin = getSpin(up_site);
    if (new_spin == getAgentSpin()) {
        isMove = true;
    }
    if (isMove) {
        update_agent_site(up_site);
        FlipOnSite();
    }
    return isMove;
}

bool SquareIceGame::go_nn_down() {
    bool isMove = false;
    int down_site = latt.NN[agent_site][1]; // down site
    int new_spin = getSpin(down_site);
    if (new_spin == getAgentSpin()) {
        isMove = true;
    }
    if (isMove) {
        update_agent_site(down_site);
        FlipOnSite();
    }
    return isMove;
}

bool SquareIceGame::go_nn_left() {
    bool isMove = false;
    int left_site = latt.NN[agent_site][2]; // left  site
    int new_spin = getSpin(left_site);
    if (new_spin == getAgentSpin()) {
        isMove = true;
    }
    if (isMove) {
        update_agent_site(left_site);
        FlipOnSite();
    }
    return isMove;
}

bool SquareIceGame::go_nn_right() {
    bool isMove = false;
    int right_site = latt.NN[agent_site][0]; // right site
    int new_spin = getSpin(right_site);
    if (new_spin == getAgentSpin()) {
        isMove = true;
    }
    if (isMove) {
        update_agent_site(right_site);
        FlipOnSite();
    }
    return isMove;
}

bool SquareIceGame::go_same_spin_up() {
    bool isMove = false;
    // Find out nieghboring candidates
    // 1. get directional up and nnn up site and spin
    int upsite = latt.NN[agent_site][3];
    int upspin = getSpin(upsite);
    int onsite_spin = getSpin(agent_site);
    int nnnupsite;
    // depend on the tpye of sublattice, choose nnn up site
    if (1 == latt.sub[agent_site]) {
        nnnupsite = latt.NNN[agent_site][0];
    } else {
        nnnupsite = latt.NNN[agent_site][1];
    }
    int nnnupspin = getSpin(nnnupsite);
    // 2. 
    if (upspin == onsite_spin) {
        // update position map
        update_agent_site(upsite);
        FlipOnSite();
        isMove = true;
    } else if (nnnupspin == onsite_spin){
        // update position map
        update_agent_site(nnnupsite);
        FlipOnSite();
        isMove = true;
    }
    return isMove;
}

bool SquareIceGame::go_same_spin_down() {
    bool isMove = false;
    // Find out nieghboring candidates
    // 1. get directional up and nnn up site and spin
    int downsite = latt.NN[agent_site][1];
    int nnndownsite;
    // depend on the tpye of sublattice, choose nnn up site
    if (1 == latt.sub[agent_site]) {
        nnndownsite = latt.NNN[agent_site][2];
    } else {
        nnndownsite = latt.NNN[agent_site][3];
    }
    int onsite_spin = getSpin(agent_site);
    int downspin    = getSpin(downsite);
    int nnndownspin = getSpin(nnndownsite);
    // 2. 
    if (downspin == onsite_spin) {
        // update position map
        update_agent_site(downsite);
        FlipOnSite();
        isMove = true;
    } else if (nnndownspin == onsite_spin){
        // update position map
        update_agent_site(nnndownsite);
        FlipOnSite();
        isMove = true;
    }
    return isMove;
}

bool SquareIceGame::go_spinup() {
    bool isMove = false;
    std::vector<int> nspins = getNeighborSpins(agent_site);
    std::vector<int> nsites = getNeighborSites(agent_site);
    for (std::vector<int>::size_type i = 1; i != nspins.size(); ++i) {
        if (nspins[i] > 0) {
            int newspin = getSpin(nsites[i]);
            if (newspin == getAgentSpin()) {
                isMove = true;
            }
            update_agent_site(nsites[i]);
            FlipOnSite();
            break;
        }
    }
    return isMove;
}

bool SquareIceGame::go_spindown() {
    bool isMove = false;
    std::vector<int> nspins = getNeighborSpins(agent_site);
    std::vector<int> nsites = getNeighborSites(agent_site);
    for (std::vector<int>::size_type i = 1; i != nspins.size(); ++i) {
        if (nspins[i] < 0) {
            int newspin = getSpin(nsites[i]);
            if (newspin == getAgentSpin()) {
                isMove = true;
            }
            update_agent_site(nsites[i]);
            FlipOnSite();
            break;
        }
    }
    return isMove;
}

bool SquareIceGame::is_short_closed() {
    bool is_closed = false;
    if (std::find(visited_map.begin(), visited_map.end(), 2) != visited_map.end()) {
        is_closed = true;
    }
    return is_closed;
}

bool SquareIceGame::is_long_closed() {
    bool is_closed = false;
    if (init_site >= 0 && init_site < MC_info.Num_sites) {
        if (visited_map[init_site] >= 2 ) {
            is_closed = true;
        }
    }
    return is_closed;
}

bool SquareIceGame::ShortLoopUpdate() {
    bool is_updated = false;
    if(is_short_closed()) {
        std::cout << "Short loop is closed!\n";
        // crossed site
        int csite = std::find(visited_map.begin(), visited_map.end(), 2) - visited_map.begin();
        // flip back those not in loop 
        while(!walked_sites.empty()) {
            int site = walked_sites.back();
            walked_sites.pop_back();
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
