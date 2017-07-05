// Try to construct the simplest MC simulation code with OOP
// initialization -> update -> measure -> output -> binning/bootstrap
// measure energy only
// without magnetic field

//#include "sample.hpp"
//#include "hamiltonian.cpp"
#include "observable.hpp"


void Observable::output_open(){
      ostringstream filename;
      filename << name << counter_replica << ".data";
      
      output.open( filename.str().c_str() );
}

void Observable::output_temperature(Sample* ss){
      output << ss->get_temperature() << '\n';
}

void Observable::output_field(Sample* ss, int n){
      output << ss->get_field(n) << '\n';
}

void Observable::output_data(double a){
      output << a << '\n';
}

void Observable::output_close(){
      output.close();
}

void Observable::output_dataline(vector<double> &a){
      vector<double>::iterator it;
      for (it = a.begin(); it != a.end(); it++)
	output << (*it) << '\t';
      output << '\n';
}

void Observable::set_counter(int n){
      counter_replica = n;
}

single_value_obs::single_value_obs(string given_name){ name = given_name; } 

single_array_obs::single_array_obs(string given_name){
  name = given_name;
  //data.resize(size);
}

void single_value_obs::measure(Sample* ss, Lattice* la){
  double result = (*obs_ptr)(ss, la);
  output_data(result);

  return;
}

void single_array_obs::measure(Sample* ss, Lattice* la){
  double result = (*obs_aptr)(ss, la, data);
  output_dataline(data); 

  return;
}
