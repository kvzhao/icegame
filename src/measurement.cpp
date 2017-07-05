// Try to construct the simplest MC simulation code with OOP
// initialization -> update -> measure -> output -> binning/bootstrap
// measure energy only
// without magnetic field

#include "measurement.hpp"

void Measurement::go(Sample* ss, Lattice* la){
      //for (int i = 0; i < Num_replicas; i++)
      //for (int i = 0; i < machine_size; i++)
	   // machine[i]->measure(ss);
      for (iter = machine.begin(); iter != machine.end(); iter++)
	  (*iter)->measure(ss, la);
}
    
void Measurement::set_counter(int n){
      //for (int i = 0; i < Num_replicas; i++)
	for (iter = machine.begin(); iter != machine.end(); iter++)
	     //for (int i = 0; i < machine_size; i++)
		  //machine[i]->set_counter(n);
		   (*iter)->set_counter(n);
}    

void Measurement::set_files(){
      //for (int i = 0; i < machine_size; i++)
	   //machine[i]->output_open(); 
      for (iter = machine.begin(); iter != machine.end(); iter++)
	(*iter)->output_open();
}
    
void Measurement::set_files_closed(){
      //for (int i = 0; i < machine_size; i++)
	    //machine[i]->output_close();
      for (iter = machine.begin(); iter != machine.end(); iter++)
	(*iter)->output_close();
}
    
void Measurement::set_files_info(Sample* ss, int index){
      //for (int i = 0; i < Num_replicas; i++){
	for (iter = machine.begin(); iter != machine.end(); iter++){
	  (*iter)->output_temperature(ss);
	  (*iter)->output_field(ss, index);
	}
	//for (int i = 0; i < machine_size; i++){
	  //    machine[i]->output_temperature(ss);
	    //  machine[i]->output_field(ss);
	//}
      //}
}

void Measurement::insert_single_obs(string obs_name){
      machine.resize(machine_size);
      machine[machine_size-1] = new single_value_obs(obs_name);
      machine_size++;
      latest++;
}

void Measurement::insert_array_obs(string obs_name){
      machine.resize(machine_size);
      machine[machine_size-1] = new single_array_obs(obs_name);
      machine_size++;
      latest++;
}

Measurement::Measurement(){ // constructor
      machine_size = 1;
      latest = -1;
}


Measurement::~Measurement(){
      //delete e;
      //delete m;
      //delete s;
      //delete a;
      //delete sc;
      machine.clear();
}

