

#include "mc.hpp"



int main(){
	
	// initialization setting  
	cout << " ------- Classical Monte Carlo Simulation ------- " << endl;	
	const int L = 32;
	int N;
	int initial = 0;
	int layer_neighbors = 2; // NN: 1, NNN: 2
	int Num_replicas;
	int Num_bins = 20;
	int Num_MCsteps = 10000;
	int Num_thermalization = 10000;
	int tempering_period = 10;	
	bool data_analysis = 1;
	vector<double> temperature;
	vector< vector<double> > magnetic_field;
	int index_vary_field = 0; // change if using field-tempering
	double t = 0.4;
	double h1 = 0;  // field strength in [111] direction (Tesla)
	double h2 = 0;   // field strength in [-1-12] direction (Tesla)
	double h3 = 0;  //field strength in [-110] direction (Tesla)
	
	double h1_t = h1 * 9.27400968 / 1.3806488;
	double h2_t = h2 * 9.27400968 / 1.3806488; // unit conversion
	double h3_t = h3 * 9.27400968 / 1.3806488;
		
	Timer TT;
	
	// import PT and setting the Monte Carlo simulation

	ifstream PTtemp("PT");
	PTtemp >> Num_replicas;	
	temperature.resize(Num_replicas);
	magnetic_field.resize(Num_replicas);
	for (int i = 0; i < Num_replicas; i++){
	  PTtemp >> temperature[i];
	  magnetic_field[i].resize(3);
	  magnetic_field[i][0] = h1_t;
	  magnetic_field[i][1] = h2_t;
	  magnetic_field[i][2] = h3_t;
	}
	PTtemp.close();

	Square lattice(L, N, layer_neighbors);	
    lattice.show_lattice();
	
	INFO MC_info(L, N, layer_neighbors, Num_replicas, Num_bins, Num_MCsteps, tempering_period, Num_thermalization);	
	
	Square_ice model(MC_info, layer_neighbors); // 1: NN only 2: up to NNN 	
	
	model.set_J1(1.0);
	
	
	
	//lattice.show_lattice();
//return 0;

// -------------------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------------------
// DESIGN the measurement
// -------------------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------------------

	// Construct and initialize the measuring machine.
	std::vector<Measurement*> Measure_sample(Num_replicas);

	// Decide which observables are needed, and link them to functions in Hamiltonian class by function pointers

	for (int i = 0; i < Num_replicas; i++){
	  Measure_sample[i] = new Measurement;
		
		Measure_sample[i]->insert_array_obs("energy");
		Measure_sample[i]->machine[Measure_sample[i]->latest]->obs_aptr = model.obs_energy;
		
		//Measure_sample[i]->insert_array_obs("magnetization");
		//Measure_sample[i]->machine[Measure_sample[i]->latest]->obs_aptr = model.magnetization;

		Measure_sample[i]->insert_array_obs("magnetization");
		Measure_sample[i]->machine[Measure_sample[i]->latest]->obs_aptr = model.obs_magnetization;
		//Measure_sample[i]->machine[Measure_sample[i]->latest]->obs_ptr = model.stagger_magnetization;
	 	
		//Measure_sample[i]->insert_array_obs("spincorrelation");
		//Measure_sample[i]->machine[Measure_sample[i]->latest]->obs_aptr = model.spin_correlation;
		
        	Measure_sample[i]->insert_array_obs("acceptance");
		Measure_sample[i]->machine[Measure_sample[i]->latest]->obs_aptr = model.obs_acceptanceSSF;
		
		Measure_sample[i]->insert_array_obs("defect");
		Measure_sample[i]->machine[Measure_sample[i]->latest]->obs_aptr = model.defect_density;

		//Measure_sample[i]->insert_array_obs("polarization");
		//Measure_sample[i]->machine[Measure_sample[i]->latest]->obs_aptr = model.FP_order; 

		//Measure_sample[i]->insert_array_obs("domainwall");
		//Measure_sample[i]->machine[Measure_sample[i]->latest]->obs_aptr = model.domain_order;

		//Measure_sample[i]->insert_array_obs("vortex");
		//Measure_sample[i]->machine[Measure_sample[i]->latest]->obs_aptr = model.H_order;
		//Measure_sample[i]->insert_array_obs("neutronmap");
		//Measure_sample[i]->machine[Measure_sample[i]->latest]->obs_aptr = model.neutron_map;

		Measure_sample[i]->set_counter(i);
		Measure_sample[i]->set_files();
	}
	
	cout << " ------- The measuring machine is constructed." << endl;


// -------------------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------------------
// INITIALIZE every sample
// -------------------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------------------

	// Construct and initialize the Ising replicas.
	std::vector<Sample*> Magnet(Num_replicas);
	for (int i = 0; i < Num_replicas; i++){
	  Magnet[i] = new Sample(MC_info);
	  model.initialization(Magnet[i], &lattice, initial);
	  model.initialize_observable(Magnet[i], &lattice, temperature[i], magnetic_field[i]);
	  Measure_sample[i]->set_files_info(Magnet[i], index_vary_field);
	}
	cout << " ------- The configuration of each replica is initialized." << endl;

// -------------------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------------------
// START the Monte Carlo Simulation
// -------------------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------------------

	TT.timer_begin();

	// thermalization
	parallel_tempering PT_heating(MC_info); 	
	PT_heating.thermalization(&model, &lattice, Magnet); // including tempering swap
	cout << " ------- The thermalization process is completed." << endl;

	//Magnet[0]->output_configuration("snapshot", 0);

	// MC update and data collect
	for (int i = 0; i < Num_MCsteps; i++){	  
		for (int j = 0; j < Num_replicas; j++){
	    model.MCstep(Magnet[j], &lattice);     // one MC step
	    
		/*	if ( (i % 50) == 0){
	    	Magnet[j]->output_configuration("snapshot", j);	
			}*/
			Measure_sample[j]->go(Magnet[j], &lattice);   // measure and store the observables
//cout << " GO!" << endl;	  
		}
	    //PT_heating.tempering(i, Magnet);    // PT swap
	}
	
	cout << " ------- The raw data production process is completed." << endl;
	
	for (int pt = 0; pt < Num_replicas; pt++){
	  PT_heating.acceptance[pt] /= double(Num_MCsteps / tempering_period); // averaging PT acceptance 
		Magnet[pt]->get_configuration(pt);     // output the final configuration of each replica
	}
	
	TT.timer_end();
	cout << " ------- Simulation time = " << TT.timer_duration() << " seconds. " << endl; 
	
	for (int i = 0; i < Num_replicas; i++){
	  Measure_sample[i]->set_files_closed();
	}

//return 0;

// -------------------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------------------
// END of Monte Carlo Simulation
// -------------------------------------------------------------------------------------------------------------

// -------------------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------------------
// DATA ANALYSIS : binning, bootstrap, output summary
// -------------------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------------------

    if (data_analysis == 1){
	Analysis analyze_energy("energy", 1, 1, MC_info);
	Analysis analyze_magnetization("magnetization", 1, 2, MC_info);
	Analysis analyze_acceptance("acceptance", 1, 0, MC_info);
	//Analysis analyze_spincorrelation("spincorrelation", L/2, 3, MC_info);
	Analysis analyze_defect("defect", 1, 0, MC_info);
	//Analysis analyze_domainwall("domainwall", 3, 0, MC_info);
	//Analysis analyze_polarization("polarization", 4, 0, MC_info);
	//Analysis analyze_scattering("neutronmap", Num_Qpoints, 3, MC_info);
	//Analysis analyze_vortex("vortex", 1, 0, MC_info);
    }
	return 0;
}


