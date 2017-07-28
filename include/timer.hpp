// Try to construct the simplest MC simulation code with OOP
// initialization -> update -> measure -> output -> binning/bootstrap
// measure energy only
// without magnetic field

#ifndef _TIMER_H_INCLUDED
#define _TIMER_H_INCLUDED

#include <time.h>

class Timer{
  private:
    clock_t start, finish;
    double duration;

  public:
    void timer_begin();
    void timer_end();
    double timer_duration();
};

#endif

