from libicegame import SQIceGame, INFO

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
sim.TEST()