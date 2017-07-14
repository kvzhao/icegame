// Try to construct the simplest MC simulation code with OOP
// initialization -> update -> measure -> output -> binning/bootstrap
// measure energy only
// without magnetic field

#ifndef _RNG_H_INCLUDED
#define _RNG_H_INCLUDED

#include <iostream>
#include <time.h>
#include <boost/random.hpp> // for mt19937
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_01.hpp>

extern boost::mt19937 rng;
extern boost::uniform_01<boost::mt19937> uni01_sampler;
#endif
