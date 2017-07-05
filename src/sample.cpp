// Try to construct the simplest MC simulation code with OOP
// initialization -> update -> measure -> output -> binning/bootstrap
// measure energy only
// without magnetic field
#include "sample.hpp"

INFO::INFO(int L, int N, int nei, int r, int b, int s, int p, int t){
	lattice_size = L;
	Num_sites = N;
    Num_neighbors = nei;
	Num_replicas = r;
	Num_bins = b;
	Num_MCsteps = s;
	Num_steps_per_bin = s / b;
	tempering_period = p;
	Num_thermalization = t;
	//cout << " ------- Monte Carlo simulation is started. -------" << endl;
	cout << " ------- Number of sites in the lattice = " << Num_sites << endl;
	cout << " ------- Number of replicas in the simulation = " << Num_replicas << endl;
	cout << " ------- Number of Monte Carlo steps for each replica = " << Num_MCsteps << endl;
	cout << " ------- Number of bins for data analysis = " << Num_bins << endl;
	cout << " ------- Number of thermalization stpes for each replica = " << Num_thermalization << endl;
	cout << " ------- Number of Monte Carlo steps in each bin = " << Num_steps_per_bin << endl;	
	print_INFO();
} 

INFO::INFO(const INFO &in) {
	lattice_size = in.lattice_size;
	Num_sites = in.Num_sites;
	Num_replicas = in.Num_replicas;
    Num_neighbors = in.Num_neighbors;
	Num_bins = in.Num_bins;
	Num_MCsteps = in.Num_MCsteps ;
	Num_steps_per_bin = in.Num_steps_per_bin;
	tempering_period = in.tempering_period;
	Num_thermalization = in.Num_thermalization;
    print_INFO();
}

void INFO::print_INFO(){
      ofstream info_out;
      info_out.open( "info");
      info_out << "lattice_size" << " " << lattice_size << endl
	       << "Num_sites" << " " << Num_sites << endl
	       << "Num_replicas" << " " << Num_replicas << endl
	       << "Num_MCsteps" << " " << Num_MCsteps << endl
	       << "Num_bins" << " " << Num_bins << endl
	       << "Num_steps_per_bin" << " " << Num_steps_per_bin << endl
	       << "tempering_period" << " " << tempering_period << endl
	       << "Num_thermalization" << " " << Num_thermalization << endl;
      info_out.close();
}

Sample::Sample(INFO MC_info){
      N = MC_info.Num_sites; // get number of lattice sites, depends on size L;
      L = MC_info.lattice_size;
      //initial_state = initial; // define the initial configuration;
      Ising.resize(N); // adjust the vector size to fit lattice sites N; 
      ctr_output = 0;
      temperature = 0;
      field.resize(1);
      field[0] = 0;
      //correlation.resize(L, 0);
}

void Sample::init(INFO MC_info){
      N = MC_info.Num_sites; // get number of lattice sites, depends on size L;
      L = MC_info.lattice_size;
      Ising.resize(N); // adjust the vector size to fit lattice sites N; 
      ctr_output = 0;
      temperature = 0;
      field.resize(1);
      field[0] = 0;
}

Sample::~Sample(){
      Ising.clear();
}
	  
void Sample::set_energy(double e){
      energy_update = e;
}

void Sample::resize_order(int a){
  order_update.resize(a);
}

void Sample::resize_magnetization(int a){
  magnetization_update.resize(a);
}

void Sample::set_magnetization(int a, double m){
    magnetization_update[a] = m;
}

void Sample::set_temperature(double t){
      temperature = t;
}

void Sample::resize_field(int n){
      field.resize(n);
}

void Sample::set_field(double h, int n){
      field[n] = h;
}
double Sample::get_temperature(){
      return temperature;
}
double Sample::get_field(int n){
      return field[n];
}

void Sample::get_configuration(int num){
  ofstream confi;
  ostringstream filename;
  filename << "IsingConfiguration" << num << ".txt";
  confi.open( filename.str().c_str() );
  confi << "temperature" << '\t' << temperature << endl;
  confi << "field" << '\t' << field[0] << endl;
  for (int i = 0; i < N; i++)
    confi << Ising[i] << endl;
  confi.close();
}

void Sample::output_configuration(string name){
  ofstream confi;
  ostringstream filename;
  filename << name << ".txt";
  confi.open( filename.str().c_str() );
  confi << "temperature" << '\t' << temperature << endl;
  confi << "field" << '\t' << field[0] << endl;
  for (int i = 0; i < N; i++)
    confi << Ising[i] << endl;
  confi.close();
  
  //ctr_output++;
}


void Sample::output_configuration(string name, int num){
  ofstream confi;
  ostringstream filename;
  filename << name << num << "_" << ctr_output << ".txt";
  confi.open( filename.str().c_str() );
  confi << "temperature" << '\t' << temperature << endl;
  confi << "field" << '\t' << field[0] << endl;
  for (int i = 0; i < N; i++)
    confi << Ising[i] << endl;
  confi.close();
  
  ctr_output++;
}

void Sample::output_configuration(string name, int num1, int num2){
  ofstream confi;
  ostringstream filename;
  filename << name <<  num1 << "_" << num2 << ".txt";
  confi.open( filename.str().c_str() );
  confi << "temperature" << '\t' << temperature << endl;
  confi << "field" << '\t' << field[0] << endl;
  for (int i = 0; i < N; i++)
    confi << Ising[i] << endl;
  confi.close();
}
