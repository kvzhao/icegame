#include "saw.hpp"

SAWMaze::SAWMaze(int L) {
    N = L*L;
    maze.resize(N, 0);
    sites_counter.resize(N, 0);
    Reset();
}

void SAWMaze::Reset() {
    std::fill(maze.begin(), maze.end(), -1);
    std::fill(sites_counter.begin(), sites_counter.end(), 0);
    traj.clear();
    traj_action.clear();
    step_counter = 0;
}

int SAWMaze::Start(int init_site) {
    set_agent_site(init_site);
    bool insane  = flip_agent_site();
    if (insane) {
        std::cout << "Problematic Initialzation\n";
    }
    return agent_site;
}

void SAWMaze::set_agent_site(int site) {
    if (site >= 0 && site < N) {
        agent_site = site;
    }
}

bool SAWMaze::flip_agent_site() {
    bool visited = false;
    int site = agent_site;
    // flip from -1 -> 1
    // then count
    if (site >= 0 && site < N) {
        maze[site] *= -1;
        if (sites_counter[site] >= 1) {
            visited = true;
        }
        sites_counter[site]++;
    }
    return visited;
}