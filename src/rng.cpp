// Try to construct the simplest MC simulation code with OOP
// initialization -> update -> measure -> output -> binning/bootstrap
// measure energy only
// without magnetic field
#include "rng.hpp"


boost::mt19937 rng(time(NULL));
boost::uniform_01<boost::mt19937> uni01_sampler(rng);
