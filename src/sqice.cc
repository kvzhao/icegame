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
    agent_site = -1;
    init_agent_site = -1;
    steps_counter = 0;
    updated_counter = 0;

    mag_fields.push_back(h1_t);
    mag_fields.push_back(h2_t);
    mag_fields.push_back(h3_t);

    // initialze all the member data
    state_0.resize(N, 0);
    state_t.resize(N, 0);
    state_tp1.resize(N, 0);
    sites_counter.resize(N, 0);

    agent_map.resize(N, 0);
    canvas_traj_map.resize(N, 0);
    canvas_spin_map.resize(N, 0);
    energy_map.resize(N, 0.0);
    defect_map.resize(N, 0.0);
    diff_map.resize(N, 0);

    std::cout << "[GAME] Square Ice Game is created.\n";
}

void SQIceGame::reset_maps() {
    std::fill(sites_counter.begin(),
              sites_counter.end(), 0);
    std::fill(canvas_traj_map.begin(),
              canvas_traj_map.end(), 0);
    std::fill(agent_map.begin(),
              agent_map.end(), 0);
    std::fill(canvas_spin_map.begin(),
              canvas_spin_map.end(), 0);
    std::fill(energy_map.begin(),
              energy_map.end(), 0.0);
    std::fill(defect_map.begin(),
              defect_map.end(), 0.0);
    traj_sites.clear();
    traj_spins.clear();
    traj_norepeat.clear();

    steps_counter = 0;
}

void SQIceGame::update_state_to_config() {
    int diff = _cal_config_difference();
    vector<int> backup = ice_config.Ising;
    ice_config.Ising = state_t;
    state_0 = state_t;
    state_tp1 = state_t;
    // state_0 = state_t 
    // ... Sanity check!
    if (   _cal_defect_density_of_state(ice_config.Ising) == 0.0 \
        && _cal_energy_of_state(ice_config.Ising) == -1.0) {
        std::cout << "[GAME] Updated Succesfully!\n";
        updated_counter++; 
        accepted_looplength.push_back(diff);
   } else {
        std::cout << "[GAME] Ice Config is RUINED. Restore.\n";
        ice_config.Ising = backup;
        restore_config_to_state();
   }
    // reset maps
}

void SQIceGame::UpdateConfig() {
    update_state_to_config();
    //?
}

void SQIceGame::restore_config_to_state() {
    state_t = ice_config.Ising;
    state_tp1 = ice_config.Ising;
    state_0 = ice_config.Ising;
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
    state_tp1 = ice_config.Ising;
    config_mean = _cal_mean(state_0);
    config_stdev = _cal_stdev(state_0);

    std::cout << "[GAME] Average Energy E = " << _cal_energy_of_state(state_0) << "\n";
    std::cout << "[GAME] Defect Density D = " << _cal_defect_density_of_state(state_0) << "\n";
    std::cout << "[GAME] Config mean = " << config_mean << " , and std = " << config_stdev << "\n";
}

int SQIceGame::Start(int init_site) {
    set_agent_site(init_site);
    return agent_site;
}

void SQIceGame::Reset() {
    reset_maps();
    init_agent_site = agent_site;
    restore_config_to_state();
}

// member functions 
vector<double> SQIceGame::Metropolis() {
    vector<double> rets(4);
    bool is_accept = false;
    double E0 = _cal_energy_of_state(state_0);
    double Et = _cal_energy_of_state(state_t);
    double dE = Et - E0;
    double dd = _cal_defect_density_of_state(state_t);
    int diff_counts = _cal_config_difference();
    double diff_ratio = diff_counts / double(N);
    if (dE == 0.0) {
        is_accept = true;
        rets[0] = 1.0;
    } else {
        rets[0] = -1.0;
    }
    rets[1] = dE;
    rets[2] = dd;
    rets[3] = diff_ratio;
    return rets;
}

void SQIceGame::flip_site(int site) {
    // maybe we need cal some info when flipping
    if(site >= 0 && site < N) {
        state_t[site] *= -1;
    }
}

void SQIceGame::flip_along_traj(const vector<int> &traj) {
    // check traj is not empty
    // check traj is cont. or not
    // flip along traj on state_t
    // done, return ?
    if(!traj.empty()) {
        for (auto const & site : traj) {
            flip_site(site);
            #ifdef DEBUG
            std::cout << "Flip site " << site 
                    << " and dE = " << _cal_energy_of_state(state_t) - _cal_energy_of_state(state_0)
                    << ", dd = " << _cal_defect_density_of_state(state_t) << endl;
            #endif
        }
    }
}

void SQIceGame::set_agent_site(int site) {
    if (site >= 0 && site < N) {
        agent_site = site;
        if (traj_sites.size() == 0) {
            init_agent_site = site;
            traj_sites.push_back(site);
            traj_norepeat.push_back(site);
        }
    } else {
        std::cout << "[GAME] WORNING, Set Agent on Illegal Site!\n";
    }
}

int SQIceGame::get_agent_site() {
    return agent_site;
}

int SQIceGame::get_spin(int site) {
    int spin = 0;
    if(site >= 0 && site < N) {
        spin = state_t[site];
    }
    return spin;
}

int SQIceGame::get_agent_spin() {
    return get_spin(agent_site);
}

vector<double> SQIceGame::Draw(int dir_idx) {
    int site = -1;
    int curt_spin = get_agent_spin();
    int next_spin = 0;
    vector<double> rets(3);
    // move agent
    switch (dir_idx) {
        case RIGHT:
            next_spin = get_spin(get_neighbor_site_by_direction(RIGHT));
            site = go(RIGHT);
        break;
        case DOWN:
            next_spin = get_spin(get_neighbor_site_by_direction(DOWN));
            site = go(DOWN);
        break;
        case LEFT:
            next_spin = get_spin(get_neighbor_site_by_direction(LEFT));
            site = go(LEFT);
        break;
        case UP:
            next_spin = get_spin(get_neighbor_site_by_direction(UP));
            site = go(UP);
        break;
        case UPPER_RIGHT:
            next_spin = get_spin(get_neighbor_site_by_direction(UPPER_RIGHT));
            site = go(UPPER_RIGHT);
        break;
        case LOWER_RIGHT:
            next_spin = get_spin(get_neighbor_site_by_direction(LOWER_RIGHT));
            site = go(LOWER_RIGHT);
        break;
        case LOWER_LEFT:
            next_spin = get_spin(get_neighbor_site_by_direction(LOWER_LEFT));
            site = go(LOWER_LEFT);
        break;
        case UPPER_LEFT:
            next_spin = get_spin(get_neighbor_site_by_direction(UPPER_LEFT));
            site = go(UPPER_LEFT);
        break;
        case UPPER_NEXT:
            next_spin = get_spin(get_neighbor_site_by_direction(UPPER_NEXT));
            site = go(UPPER_NEXT);
        break;
        case LOWER_NEXT:
            next_spin = get_spin(get_neighbor_site_by_direction(LOWER_NEXT));
            site = go(LOWER_NEXT);
        break;
    }
    // draw on canvas
    if (canvas_traj_map[site] == 1) {
        canvas_traj_map[site] = 0;
    } else {
        canvas_traj_map[site] = 1;
    }
    canvas_spin_map[site] = state_t[site];

    double dE = _cal_energy_of_state(state_tp1) - _cal_energy_of_state(state_t);
    double dd = _cal_defect_density_of_state(state_tp1);

    if (curt_spin == next_spin ) {
        rets[0] = -1.0;
    } else {
        rets[0] = +1.0;
    }
    rets[1] = dE; rets[2] = dd;

    #ifdef DEBUG 
    std::cout << "  current spin " << curt_spin << " , next spin " << next_spin << "\n";
    std::cout << " Draw dE = " << dE
                << ", dd = " << dd << "\n";
    #endif

    return rets;
}

void SQIceGame::IceMove(int act_idx) {
    int dir_idx; 
    int agent_spin = get_agent_spin();

    switch (act_idx) {
        // go spin up
        case 0:
            if (agent_spin > 0) {
                dir_idx = how_to_go(icemove(true));
            } else {
                dir_idx = how_to_go(icemove(false));
            }
            Draw(dir_idx);
        break;
        // go spin down
        case 1:
            if (agent_spin < 0) {
                dir_idx = how_to_go(icemove(true));
            } else {
                dir_idx = how_to_go(icemove(false));
            }
            Draw(dir_idx);
        break;
    }
}

int SQIceGame::go(DIR dir) {
    int old_site = agent_site;
    int new_site = get_neighbor_site_by_direction(dir);
    set_agent_site(new_site);
    if (!_is_visited(new_site)) {
        traj_norepeat.push_back(new_site);
    }
    traj_spins.push_back(get_spin(new_site));
    traj_sites.push_back(new_site);
    action_list.push_back(dir);
    // no repeated trajectory

    // agent map
    agent_map[old_site] = 0;
    agent_map[new_site] = 1; 

    if (sites_counter[agent_site] < 1) {
        state_tp1[agent_site] *= -1;
    }

    steps_counter++;
    sites_counter[agent_site]++;

    return agent_site;
}

DIR SQIceGame::how_to_go(int site) {
    DIR dir = NOWAY;
    dir = get_direction_by_sites(agent_site, site);
    return dir;
}

int SQIceGame::icemove(bool same_spin) {
    int site = -1;
    vector<int> cands = get_neighbor_candidates(same_spin);
    if (!cands.empty()) {
        site = cands[0];
    }
    return site;
}

vector<int> SQIceGame::get_neighbor_sites() {
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

vector<int> SQIceGame::get_neighbor_spins() {
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

vector<int> SQIceGame::get_neighbor_candidates(bool same_spin) {
    int inverse_agent_spin = 0;
    if (same_spin) {
        inverse_agent_spin = get_agent_spin();
    } else {
        inverse_agent_spin = -1 * get_agent_spin();
    }
    vector<int> nspins = get_neighbor_spins();
    vector<int> nsites = get_neighbor_sites();
    vector<int> idxes;
    for (std::vector<int>::size_type i =0; i < nspins.size(); i++) {
        if (inverse_agent_spin == nspins[i]) {
            idxes.push_back(i);
        }
    }
    vector<int> candidates;
    for (auto const & idx : idxes) {
        candidates.push_back(nsites[idx]);
    }
    std::random_shuffle(candidates.begin(), candidates.end());
    // should i avoid repeating?
    return candidates;
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
            site = latt.NNN[site][1];
    } else if (dir == UPPER_LEFT) {
            site = latt.NNN[site][2];
    } else if (dir == UPPER_RIGHT) {
            site = latt.NNN[site][3];
    } else if (dir == LOWER_NEXT) {
        if (latt.sub[site] == 1) {
            site = latt.NNN[site][0];
        } else {
            site = latt.NNN[site][1];
        }
    } else if (dir == UPPER_NEXT) {
        if (latt.sub[site] == 1) {
            site = latt.NNN[site][2];
        } else {
            site = latt.NNN[site][3];
        }
    }
    return site;
}

DIR SQIceGame::get_direction_by_sites(int site, int next_site) {
    DIR dir = NOWAY;
    int right_site = latt.NN[site][0];
    int left_site = latt.NN[site][2];
    int up_site = latt.NN[site][3];
    int down_site = latt.NN[site][1];
    int lower_right_site = latt.NNN[site][0];
    int lower_left_site  = latt.NNN[site][1];
    int upper_left_site  = latt.NNN[site][2];
    int upper_right_site = latt.NNN[site][3];

    // check next state is in its neighbots
    if (next_site == right_site) {
        dir = RIGHT;
    } else if (next_site == left_site) {
        dir = LEFT;
    } else if (next_site == up_site) {
        dir = UP;
    } else if (next_site == down_site) {
        dir  = DOWN;
    } else if (next_site == upper_right_site) {
        dir = UPPER_RIGHT;
    } else if (next_site == lower_right_site) {
        dir = LOWER_RIGHT;
    } else if (next_site == lower_left_site) {
        dir = LOWER_LEFT;
    } else if (next_site == upper_left_site) {
        dir = UPPER_LEFT;
    }

    #ifdef DEBUG
    std::cout << "right site = " << right_site << "\n"
              << "left site = " << left_site << "\n"
              << "up site = " << up_site << "\n"
              << "down site = " << down_site << "\n"
              << "upper right site = " << upper_right_site << "\n"
              << "lower right site = " << lower_right_site << "\n"
              << "lower left  site = " << lower_left_site << "\n"
              << "upper left  site = " << upper_left_site << "\n";
    #endif

    return dir;
}

vector<double> SQIceGame::GetEnergyMap() {
    for (uint i = 0; i < N; i++) {
        double se = 0.0;
        se = _cal_energy_of_site(state_t, i);
        energy_map[i] = se/6.0;
    }
    return energy_map;
}

vector<double> SQIceGame::GetStateTMap() {
    vector<double> map(state_t.begin(), state_t.end());
    return map;
}

vector<double> SQIceGame::GetCanvasMap() {
    vector<double> map(canvas_traj_map.begin(),
                        canvas_traj_map.end());
    return map;
}

vector<double> SQIceGame::GetDefectMap() {
    for (uint i =0; i < N; i++) {
        double dd = 0.0;
        if (latt.sub[i] == 1) {
            dd = state_t[i] + state_t[latt.NN[i][0]] + state_t[latt.NN[i][1]] + state_t[latt.NNN[i][0]];
            dd /= 4.0;
        }
        defect_map[i] = dd;
    }
    return defect_map;
}

void SQIceGame::TEST(){
    std::cout << "/=========================================/\n";

    set_agent_site(5000);
    std::cout << "Set Agent site at 5000, agent = " << get_agent_site() << "\n";
    set_agent_site(100);
    std::cout << "Set Agent site at 100, agent = " << get_agent_site() << "\n";
    std::cout << "Its neighbors are ";
    _print_vector(get_neighbor_sites());
    std::cout << ".... with corresponding spins are ";
    _print_vector(get_neighbor_spins());

    std::cout << " init agent site = " << init_agent_site << "\n";

    std::cout << "\twhose up site is " << get_neighbor_site_by_direction(UP) << "\n";
    std::cout << "\twhose down site is " << get_neighbor_site_by_direction(DOWN) << "\n";
    std::cout << "\twhose left site is " << get_neighbor_site_by_direction(LEFT) << "\n";
    std::cout << "\twhose right site is " << get_neighbor_site_by_direction(RIGHT) << "\n";
    std::cout << "\twhose upper right site is " << get_neighbor_site_by_direction(UPPER_RIGHT) << "\n";
    std::cout << "\twhose lower right site is " << get_neighbor_site_by_direction(LOWER_RIGHT) << "\n";
    std::cout << "\twhose lower left site is " << get_neighbor_site_by_direction(LOWER_LEFT) << "\n";
    std::cout << "\twhose upper left site is " << get_neighbor_site_by_direction(UPPER_LEFT) << "\n";
    std::cout << "\twhose upper next site is " << get_neighbor_site_by_direction(UPPER_NEXT) << "\n";
    std::cout << "\twhose lower next site is " << get_neighbor_site_by_direction(LOWER_NEXT) << "\n";
    std::cout << "/=========================================/\n";
    long_loop_algorithm();

    /*
    Draw(0);
    _print_vector(get_neighbor_candidates());
    Draw(0);
    _print_vector(get_neighbor_candidates());
    Draw(3);
    Draw(3);
    Draw(2);
    Draw(2);
    Draw(1);
    Draw(1);
    std::cout << " init agent site = " << init_agent_site << "\n";
    DIR d = get_direction_by_sites(traj_sites[0], traj_sites[1]);
    std::cout << "From " << traj_sites[0] << " to " << traj_sites[1] 
                << " is done by action " << d << "\n"; 
    if(_is_traj_continuous()) {
        std::cout << "Loop is cont\n";
    } else {
        std::cout << "Loop is NOT cont\n";
    }
    if(_is_start_end_meets(agent_site)) {
        std::cout << "End meets!\n";
    }
    std::cout << "traj: ";
    _print_vector(traj_sites);
    std::cout << "no repetition traj: ";
    _print_vector(traj_norepeat);

    //std::cout << "Agent site after ice move : " << icemove() << "\n";
    //std::cout << "Agent site after ice move : " << icemove() << "\n";
    //std::cout << "Agent site after ice move : " << icemove() << "\n";
    //std::cout << "Agent site after ice move : " << icemove() << "\n";
    //std::cout << "Agent site after ice move : " << icemove() << "\n";
    flip_along_traj(traj_norepeat);
    */
}

// ###### private function ######
double SQIceGame::_cal_energy_of_state(const vector<int> & s) {
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

double SQIceGame::_cal_energy_of_site(const vector<int> &s, int site) {
    double se = 0.0;
    int i = site;
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
    se = J1 * se;
    return se;
}

double SQIceGame::_cal_defect_density_of_state(const vector<int> & s) {
    double dd = 0.0;
    for (uint i = 0; i < N; ++i) {
        if (latt.sub[i] == 1) {
            dd += s[i] + s[latt.NN[i][0]] + s[latt.NN[i][1]] + s[latt.NNN[i][0]];
        }
    }
    dd /= 2.0;
    dd /= N;
    return abs(dd);
}

double SQIceGame::_cal_mean(const vector<int> &s) {
    return 1.0* std::accumulate(s.begin(), s.end(), 0LL) / s.size();
}

double SQIceGame::_cal_stdev(const vector<int> &s) {
    double mean = _cal_mean(s);
    vector<double> diff(s.size());
    std::transform(s.begin(), s.end(), diff.begin(), [mean](double x) { return x - mean; });
    double sq_sum = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
    double stdev = std::sqrt(sq_sum/ s.size());
    return stdev;
}

void SQIceGame::_print_vector(const vector<int> &v) {
    std::cout << "[";
    for (const auto &i: v)
        std::cout << i << ", ";
    std::cout << "]\n";
}

bool SQIceGame::_is_visited(int site) {
    bool visited = false;
    if (std::find(traj_sites.begin(), 
                    traj_sites.end(), site) != traj_sites.end()) {
        visited = true;
    }
    return visited;
}

bool SQIceGame::_is_start_end_meets(int site) {
    return (site == init_agent_site) ? true : false;
}

bool SQIceGame::_is_traj_continuous() {
    bool cont = true;
    // try walk through traj with operation
    // check next site is in its neighbors
    for (std::vector<int>::size_type i = 1;  i < traj_sites.size(); i++) {
        std::cout << "step " << i << endl;
        DIR dir = get_direction_by_sites(traj_sites[i-1], traj_sites[i]); 
        std::cout << "From " << traj_sites[i-1] << " to " << traj_sites[i] 
                << " is done by action " << dir << "\n"; 
        if (dir == NOWAY) {
            cont = false;
            break;
        }
    }
    return cont;
}

int SQIceGame::_cal_config_difference() {
    int diff_counter = 0;
    for (size_t i = 0 ; i < N; i++) {
        if (state_0[i] != state_t[i]) {
            diff_counter++;
            diff_map[i] = 1;
        }
    }
    return diff_counter;
}

bool SQIceGame::_is_long_loop() {
    bool loop_occur = false;
    return loop_occur;
}

bool SQIceGame::_is_short_loop() {
    bool loop_occur = false;
    return loop_occur;
}

void SQIceGame::long_loop_algorithm() {
    set_agent_site(100);
    int next_site = icemove(false);
    while (!_is_start_end_meets(next_site)) {
        Draw(how_to_go(next_site));
        next_site = icemove(false);
    }
    Draw(how_to_go(next_site));
    _print_vector(traj_norepeat);
    flip_along_traj(traj_sites);
    //flip_along_traj(traj_norepeat);
}