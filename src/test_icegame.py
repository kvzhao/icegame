from icegame import SQIceGame, INFO
import numpy as np

# physical parameters
L = 32
kT = 0.0001
J = 1
N = L**2

num_neighbors = 2
num_replicas = 1
num_mcsteps = 2000
num_bins = 1
num_thermalization = num_mcsteps
tempering_period = 1

mc_info = INFO(L, N, num_neighbors, num_replicas, num_bins, num_mcsteps, tempering_period, num_thermalization)

# initalize the system, lattice config
sim = SQIceGame(mc_info)
sim.set_temperature (kT)
sim.init_model()
sim.mc_run(num_mcsteps)

sim.start(100)

eng_map = sim.get_energy_map()
print(eng_map)
print(type(eng_map))

for i in range(10):
    print (sim.draw(np.random.randint(6)))
    print (sim.get_trajectory())
sites = sim.get_trajectory()
site_diffs = [j-i for i, j in zip(sites[:-1], sites[1:])]
sim.flip_trajectory()
print(sim.metropolis())

print ('start point = {}'.format(sim.get_start_point()))
print ('Execute {} steps in episode'.format(sim.get_ep_step_counter()))
print ('Action counter in this ep = {}'.format(sim.get_ep_action_counters()))
print ('Action list in this ep = {}'.format(sim.get_ep_action_list()))
print ('Action Statistics = {}'.format(sim.get_action_statistics()))

