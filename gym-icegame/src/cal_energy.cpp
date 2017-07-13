

#include "cal_energy.hpp"



int main(){
	
	// initialization setting  
	cout << " ------- Classical Monte Carlo Simulation ------- " << endl;	
	const int L = 60;
	int N;
	int initial = 0;
	int layer_neighbors = 3; // NN: 1, NNN: 2
	int Num_replicas;
	//int Num_data = 0;
	int Num_bins = 20;
	int Num_MCsteps = 10000;
	int Num_thermalization = 10000;
	int tempering_period = 10;	
	bool data_analysis = 1;
	//double field = 0;
	// 1.75, -0.25
	double h1 = 1.5;  // field strength in [111] direction (Tesla)
	double h2 = -0.25;   // field strength in [-1-12] direction (Tesla)
	double h3 = 0;  //field strength in [-110] direction (Tesla)
	double h1_t = h1 * 9.27400968 / 1.3806488;
	double h2_t = h2 * 9.27400968 / 1.3806488; // unit conversion
	double h3_t = h3 * 9.27400968 / 1.3806488;
 
	//double h_frequency = 100.0;
	//h_frequency *= 2.0 * 3.14159265358979323846;
		
	Timer TT;
	// import PT and setting the Monte Carlo simulation
	ifstream PTtemp("PT");
	PTtemp >> Num_replicas;	
	
	Kagome lattice(L, N, layer_neighbors);	
	
	MC_info.set_INFO(L, N, Num_replicas, Num_bins, Num_MCsteps, tempering_period, Num_thermalization);	
	
	Kagome_ice model(MC_info, layer_neighbors); // 1: NN only 2: up to NNN 	
	model.set_J1(1.0);
	model.set_J2(-1.0/3.0);
	//model.set_J3a(0.05);
	//model.set_J3b(0.05);	
	
	
	//lattice.Hstate_Bloch_phase();
	
	//lattice.show_lattice();
//return 0;

// -------------------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------------------
// DESIGN the measurement
// -------------------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------------------
/*
	// Construct and initialize the measuring machine.
	std::vector<Measurement*> Measure_sample(Num_replicas);

	// Decide which observables are needed, and link them to functions in Hamiltonian class by function pointers

	for (int i = 0; i < Num_replicas; i++){
	  Measure_sample[i] = new Measurement;
		
		Measure_sample[i]->insert_array_obs("energy");
		Measure_sample[i]->machine[Measure_sample[i]->latest]->obs_aptr = model.obs_energy;
		
		Measure_sample[i]->insert_array_obs("magnetization");
		Measure_sample[i]->machine[Measure_sample[i]->latest]->obs_aptr = model.magnetization;

		//Measure_sample[i]->insert_single_obs("magnetization");
		//Measure_sample[i]->machine[Measure_sample[i]->latest]->obs_ptr = model.obs_magnetization;
		//Measure_sample[i]->machine[Measure_sample[i]->latest]->obs_ptr = model.stagger_magnetization;
		

		Measure_sample[i]->insert_array_obs("acceptance");
		Measure_sample[i]->machine[Measure_sample[i]->latest]->obs_aptr = model.obs_acceptanceSSF;
		
		Measure_sample[i]->insert_array_obs("defect");
		Measure_sample[i]->machine[Measure_sample[i]->latest]->obs_aptr = model.defect_density;

		//Measure_sample[i]->insert_array_obs("polarization");
		//Measure_sample[i]->machine[Measure_sample[i]->latest]->obs_aptr = model.FP_order; 

		Measure_sample[i]->insert_array_obs("domainwall");
		Measure_sample[i]->machine[Measure_sample[i]->latest]->obs_aptr = model.domain_order;

		Measure_sample[i]->insert_array_obs("vortex");
		Measure_sample[i]->machine[Measure_sample[i]->latest]->obs_aptr = model.H_order;

		Measure_sample[i]->set_counter(i);
		Measure_sample[i]->set_files();
	}
	
	cout << " ------- The measuring machine is constructed." << endl;
*/

// -------------------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------------------
// INITIALIZE every sample
// -------------------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------------------
/*
	// Construct and initialize the Ising replicas.
	std::vector<Sample*> Magnet(Num_replicas);
	for (int i = 0; i < Num_replicas; i++){
	  Magnet[i] = new Sample(MC_info);
	  model.initialization(Magnet[i], &lattice, initial);
	  double t;
	  PTtemp >> t;
	  Magnet[i]->set_temperature(t);
	  model.set_field(Magnet[i], h1_t, h2_t, h3_t);
	  Measure_sample[i]->set_files_info(Magnet[i]);
	  model.initialize_observable(Magnet[i], &lattice);
	}
	PTtemp.close();	
	cout << " ------- The configuration of each replica is initialized." << endl;
*/
// -------------------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------------------
// OUTPUT some special configuration for checking (not necessary)
// -------------------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------------------
	/*
	int Num_states = 1;
	std::vector<Sample*> Candidate(Num_states);
	for (int i = 0; i < Num_states; i++){
		Candidate[i] = new Sample(MC_info);
		model.initialization(Candidate[0], &lattice, 2);
		Candidate[0]->set_temperature(0);
		model.set_field(Candidate[0], h1_t, h2_t, h3_t);
	}
	for (double ja = -0.2; ja <= 0.2; ja+=0.01){
	  for (double jb = -0.2; jb <= 0.2; jb+=0.01){
			model.set_J3a(ja);
	    model.set_J3b(jb);
	    double eng = model.total_energy(Candidate[0], &lattice);
	    eng /= double(N);
	    cout << "J3a = " << ja << ", J3b = " << jb << ", E = " << eng << endl;
	  }
	}
	
	*/
	
	/*
	Sample qx_file(MC_info);
	model.initialization(&qx_file, &lattice, 2);
	qx_file.output_configuration("qx");
	//model.defect_configuration(Magnet[0], &lattice, "monopole", 0);
	//Magnet[0]->output_configuration("CheckConfiguration", 0);
	
	model.degenerate_qx(Magnet[0], &lattice, 0);	
	model.defect_configuration(Magnet[0], &lattice, "monopole", 0);
	Magnet[0]->output_configuration("CheckConfiguration", 0);
	
	model.degenerate_qx(Magnet[0], &lattice, 1);	
	model.defect_configuration(Magnet[0], &lattice, "monopole", 0);
	Magnet[0]->output_configuration("CheckConfiguration", 0);
	

	model.degenerate_vortex(Magnet[0], &lattice, 0);	
	model.defect_configuration(Magnet[0], &lattice, "monopole", 0);
	Magnet[0]->output_configuration("CheckConfiguration", 0);
	
	model.degenerate_vortex(Magnet[0], &lattice, 1);	
	model.defect_configuration(Magnet[0], &lattice, "monopole", 0);
	Magnet[0]->output_configuration("CheckConfiguration",0);

	model.degenerate_vortex(Magnet[0], &lattice, 2);	
	model.defect_configuration(Magnet[0], &lattice, "monopole", 0);
	Magnet[0]->output_configuration("CheckConfiguration",0);
		
	//return 0;
*/
// -------------------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------------------
// START the Monte Carlo Simulation
// -------------------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------------------
/*
	TT.timer_begin();

	// thermalization
	parallel_tempering PT_heating(MC_info); 	
	PT_heating.thermalization(&model, &lattice, Magnet); // including tempering swap
	cout << " ------- The thermalization process is completed." << endl;

	//model.defect_configuration(Magnet[0], &lattice, "monopole", 0);
	//Magnet[0]->output_configuration("snapshot", 0);

	// MC update and data collect
	for (int i = 0; i < Num_MCsteps; i++){	  
	  //double h_update = h1_t * cos(h_frequency * i);
		for (int j = 0; j < Num_replicas; j++){
			//model.set_field(Magnet[j], h_update, h2_t, h3_t);
	    model.MCstep(Magnet[j], &lattice);     // one MC step
			
			Measure_sample[j]->go(Magnet[j], &lattice);   // measure and store the observables
	  }
	    //PT_heating.tempering(i, Magnet);    // PT swap
	}
	
	cout << " ------- The raw data production process is completed." << endl;
	
	for (int pt = 0; pt < Num_replicas; pt++){
	  PT_heating.acceptance[pt] /= double(Num_MCsteps / tempering_period); // averaging PT acceptance 
		model.defect_configuration(Magnet[pt], &lattice, "monopole", pt);
		Magnet[pt]->get_configuration(pt);     // output the final configuration of each replica
	}
	
	TT.timer_end();
	cout << " ------- Simulation time = " << TT.timer_duration() << " seconds. " << endl; 
	
	for (int i = 0; i < Num_replicas; i++){
	  Measure_sample[i]->set_files_closed();
	}

//return 0;
*/
// -------------------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------------------
// END of Monte Carlo Simulation
// -------------------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------------------
// CHECK the final configuration at lowest temperature, compare to known states
// -------------------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------------------
	TT.timer_begin();

	vector<Sample*> check_state;
	check_state.resize(5);
	for (int i = 0; i < 5; i++){
		check_state[i] = new Sample(MC_info);
		check_state[i]->set_temperature(1);
		model.set_product_hs(check_state[i], h1_t, h2_t, h3_t);
	}


	model.degenerate_vortex(check_state[0], &lattice, 0);
	//model.degenerate_vortex(check_state[1], &lattice, 1);
	//model.degenerate_vortex(check_state[2], &lattice, 2);
	model.degenerate_qx(check_state[1], &lattice, 0);
	//model.degenerate_qx(check_state[4], &lattice, 1);
	model.initialization(check_state[2], &lattice, 3);
	model.monopole_vortex(check_state[3], &lattice, 0);
	model.initialization(check_state[4], &lattice, 4);
	//model.monopole_vortex(check_state[7], &lattice, 1);
	//model.monopole_vortex(check_state[8], &lattice, 2);
	//model.anti_vortex(check_state[9], &lattice, 0);
	//model.anti_vortex(check_state[10], &lattice, 1);
	//model.anti_vortex(check_state[11], &lattice, 2);

	ofstream energy_phase;
	energy_phase.open("Summary_phase");
	energy_phase << "ja" << '\t'
		           << "jb" << '\t'
							 << "phase" << '\n';

	/*
	energy_phase << "ja" << '\t'
		           << "jb" << '\t'
							 << "v1" << '\t'
							 << "v2" << '\t'
							 << "v3" << '\t'
							 << "qx1" << '\t'
							 << "qx2" << '\t'
							 << "sat" << '\t'
							 << "dv1" << '\t'
							 << "dv2" << '\t'
							 << "dv3" << '\t'
							 << "v4" << '\t'
							 << "v5" << '\t'
							 << "v6" << '\n';
*/
	vector<double> energy_compare(5);
	double e_min = 0;
	int emin_index;
	for (double jb = -0.2; jb <= 0.2; jb+=0.005){
	  for (double ja = -0.2; ja <= 0.2; ja+=0.005){
			model.set_J3a(ja);
			model.set_J3b(jb);
			e_min = 0;
			energy_phase << ja << '\t';
			energy_phase << jb << '\t';
			for (int i = 0; i < 5; i++){
	    	energy_compare[i] = model.total_energy(check_state[i], &lattice) / double(N);
				if (e_min > energy_compare[i]){
					e_min = energy_compare[i];
					emin_index = i;
				}
				//double eng = model.total_energy(check_state[i], &lattice);
				//eng /= double(N);
				//energy_phase << eng << '\t';
			}
			energy_phase << emin_index << '\n';
		}
	}


	energy_phase.close();

	TT.timer_end();
	cout << " ------- Calculation time = " << TT.timer_duration() << " seconds. " << endl; 


	//vector<double> domain_order;
	
	/*
	cout << " ****** vortex state (1) = " << model.compare_state(check_state[0], Magnet[0]) << endl;
	cout << " ****** vortex state (2) = " << model.compare_state(check_state[1], Magnet[0]) << endl;
	cout << " ****** vortex state (3) = " << model.compare_state(check_state[2], Magnet[0]) << endl;
	
	cout << " ****** vortex state (4) = " << model.compare_state(check_state[9], Magnet[0]) << endl;
	cout << " ****** vortex state (5) = " << model.compare_state(check_state[10], Magnet[0]) << endl;
	cout << " ****** vortex state (6) = " << model.compare_state(check_state[11], Magnet[0]) << endl;

	cout << " ****** q = X state  (1) = " << model.compare_state(check_state[3], Magnet[0]) << endl;
	cout << " ****** q = X state  (2) = " << model.compare_state(check_state[4], Magnet[0]) << endl;
	cout << " ****** saturated state  = " << model.compare_state(check_state[5], Magnet[0]) << endl;
	cout << " ***** defect vortex (1) = " << model.compare_state(check_state[6], Magnet[0]) << endl;
	cout << " ***** defect vortex (2) = " << model.compare_state(check_state[7], Magnet[0]) << endl;
	cout << " ***** defect vortex (3) = " << model.compare_state(check_state[8], Magnet[0]) << endl;
	
  model.domain_order(Magnet[0], &lattice, domain_order);
	cout << " ***** domain spin order = " << domain_order[0] << endl;
	cout << " ***** alpha chain order = " << domain_order[1] << endl;
	cout << " ***** beta chain order  = " << domain_order[2] << endl;
*/
	check_state.clear();
//return 0;


	//cout << " ########## total energy = " << Magnet[0]->energy_update / double(MC_info.Num_sites) << endl;


// -------------------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------------------
// DATA ANALYSIS : binning, bootstrap, output summary
// -------------------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------------------
/*
    if (data_analysis == 1){
	Analysis analyze_energy("energy", 1, 1, MC_info);
	Analysis analyze_magnetization("magnetization", 4, 2, MC_info);
	Analysis analyze_acceptance("acceptance", 1, 0, MC_info);
	Analysis analyze_defect("defect", 1, 0, MC_info);
	Analysis analyze_domainwall("domainwall", 3, 0, MC_info);
	Analysis analyze_vortex("vortex", 1, 0, MC_info);
    }

*/
	return 0;
}


