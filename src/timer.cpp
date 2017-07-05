// Try to construct the simplest MC simulation code with OOP
// initialization -> update -> measure -> output -> binning/bootstrap
// measure energy only
// without magnetic field

#include "timer.hpp"

void Timer::timer_begin(){
      start = clock();
}

void Timer::timer_end(){
      finish = clock();
}

double Timer::timer_duration(){
      duration = (double)(finish-start) / CLOCKS_PER_SEC;
      return duration;
}

