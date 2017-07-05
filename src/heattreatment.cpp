// Try to construct the simplest MC simulation code with OOP
// initialization -> update -> measure -> output -> binning/bootstrap
// measure energy only
// without magnetic field

#include "heattreatment.hpp"


void HeatTreatment::import_info(const INFO &MC_info){
      Num_replicas = MC_info.Num_replicas;
      Num_thermalization = MC_info.Num_thermalization;
      temperature_period = MC_info.tempering_period;
}


parallel_tempering::parallel_tempering(){
}

parallel_tempering::parallel_tempering(const INFO &MC_info){
      import_info(MC_info);
      acceptance.resize(Num_replicas, 0);
}  

void parallel_tempering::init (INFO MC_info){
      import_info(MC_info);
      acceptance.resize(Num_replicas, 0);
}  
    
void parallel_tempering::reset_accptance(){
      acceptance.clear();
      acceptance.resize(Num_replicas, 0);
}

void parallel_tempering::thermalization(Hamiltonian* model, Lattice* la, vector<Sample *> vs){
      for (int i = 0; i < Num_thermalization; i++){
	for (int j = 0; j < Num_replicas; j++){	  
	  model->MCstep(vs[j], la);
	}
	tempering(i, vs);
      }

      reset_accptance();
}

void parallel_tempering::thermalization2(Hamiltonian* model, Lattice* la, vector<Sample *> vs){
      for (int i = 0; i < Num_thermalization; i++){
	for (int j = 0; j < Num_replicas; j++){	  
	  model->MCstep(vs[j], la);
	}
      }

      reset_accptance();
}



void parallel_tempering::tempering(int n, vector<Sample *> vs){
      
      if (n % temperature_period == 0){
      
	for (int i = 0; i < (Num_replicas - 1); i++){		
	  double delta_e = vs[i]->energy_update - vs[i+1]->energy_update;
	  double delta_t = ( 1.0 / vs[i]->get_temperature() ) - ( 1.0 / vs[i+1]->get_temperature() );
	  double delta = delta_e * delta_t;
	  double weight = exp(delta);
	  double dice = uni01_sampler();

	  if (dice < weight){
	    vs[i]->Ising.swap(vs[i+1]->Ising);
	    swap(vs[i]->energy_update, vs[i+1]->energy_update);
	    swap(vs[i]->magnetization_update, vs[i+1]->magnetization_update);
	    acceptance[i] ++;
	  }
      
	}
      }
} 

simulated_annealing::simulated_annealing(INFO MC_info){
      import_info(MC_info);
      current_replica = (Num_replicas - 1);
}  

void simulated_annealing::thermalization(Hamiltonian* model, Lattice* la, vector<Sample *> vs){
      for (int i = 0; i < Num_thermalization; i++){
	  model->MCstep(vs[current_replica], la);
	}
}

void simulated_annealing::annealing(int n, int current, vector<Sample*> vs){
      if (n % temperature_period == 0){
	if ( current != 0){
	  vs[current-1]->Ising = vs[current]->Ising;
	  vs[current-1]->energy_update = vs[current]->energy_update;
	  vs[current-1]->magnetization_update = vs[current]->magnetization_update;
	}
      }
}


field_tempering::field_tempering(INFO MC_info, int n){
      import_info(MC_info);
      acceptance.resize(Num_replicas, 0);
      index = n; // vary which magnetic field component
}  
    
void field_tempering::reset_accptance(){
      acceptance.clear();
      acceptance.resize(Num_replicas, 0);
}

void field_tempering::thermalization(Hamiltonian* model, Lattice* la, vector<Sample *> vs){
      for (int i = 0; i < Num_thermalization; i++){
	for (int j = 0; j < Num_replicas; j++){	  
	  model->MCstep(vs[j], la);
	}
	tempering(i, vs);
      }

      reset_accptance();
}

void field_tempering::thermalization2(Hamiltonian* model, Lattice* la, vector<Sample *> vs){
      for (int i = 0; i < Num_thermalization; i++){
	for (int j = 0; j < Num_replicas; j++){	  
	  model->MCstep(vs[j], la);
	}
      }

      reset_accptance();
}



void field_tempering::tempering(int n, vector<Sample *> vs){
      
      if (n % temperature_period == 0){
      
	for (int i = 0; i < (Num_replicas - 1); i++){		
          double delta_m = vs[i]->magnetization_update[index] - vs[i+1]->magnetization_update[index];
          double delta_h = vs[i]->get_field(index) - vs[i+1]->get_field(index);
          double delta = (-1.0) * delta_m * delta_h / vs[i]->get_temperature();
	  double weight = exp(delta);
	  double dice = uni01_sampler();

	  if (dice < weight){
	    vs[i]->Ising.swap(vs[i+1]->Ising);
	    swap(vs[i]->energy_update, vs[i+1]->energy_update);
	    swap(vs[i]->magnetization_update, vs[i+1]->magnetization_update);
	    acceptance[i] ++;
	  }
      
	}
      }
} 
