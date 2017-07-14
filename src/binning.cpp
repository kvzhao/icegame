// Try to construct the simplest MC simulation code with OOP
// initialization -> update -> measure -> output -> binning/bootstrap
// measure energy only
// without magnetic field

#include "binning.hpp"
/*
Bin_analysis::Bin_analysis(string obs, int s, int i, INFO mc_info){
    
      input_open();
      set_parameters(obs, s, i ,mc_info);
}*/

void Bin_analysis::set_parameters(string obs, int s, int i, INFO mc_info){

      Num_sites = mc_info.Num_sites;
      Num_replicas = mc_info.Num_replicas;
      Num_bins = mc_info.Num_bins;
      size_bin = mc_info.Num_steps_per_bin;
      counter_replica = i;
      name = obs;
      size_data = s;
      Num_tests = 1000;
      
      data.resize(size_data, vector<double> (Num_bins, 0));
      data2.resize(size_data, vector<double> (Num_bins, 0));
      data4.resize(size_data, vector<double> (Num_bins, 0));
      data_response.resize(size_data, vector<double> (Num_bins, 0));
      data_Binder.resize(size_data, vector<double> (Num_bins, 0));
      data_abs.resize(size_data, vector<double> (Num_bins, 0));
      
      data_bootstrap.resize(size_data, vector<double> (Num_tests, 0));
      data_response_bootstrap.resize(size_data, vector<double> (Num_tests, 0));
      data_Binder_bootstrap.resize(size_data, vector<double> (Num_tests, 0));
}

Bin_analysis::~Bin_analysis(){
    //data.clear(); data2.clear(); data4.clear(); data_abs.clear();
    //data_response.clear(); data_Binder.clear(); data_bootstrap.clear();
    //data_response_bootstrap.clear(); data_Binder_bootstrap.clear();
}
 
void Bin_analysis::input_open(){  
      ostringstream filename;
      filename << name << counter_replica << ".data";
      input.open( filename.str().c_str() );
      input >> temperature;
      input >> field;
}
    
void Bin_analysis::input_close(){
      input.close();
}

double Bin_analysis::get_temperature(){
      return temperature;
}

double Bin_analysis::get_field(){
      return field;
}

void Bin_analysis::binning(){  //do average only, for dataline(list of data)
            
      for (int i = 0; i < Num_bins; i++){
	for (int j = 0; j < size_bin; j++){
	  for (int k = 0; k < size_data; k++){
	    double dd;
	    input >> dd;
	    data[k][i] += dd;
	    data2[k][i] += dd * dd;
	    data4[k][i] += pow(dd,4);
	    data_abs[k][i] += abs(dd);
	  }
	}
	for (int l = 0; l < size_data; l++){
	  data[l][i] /= double(size_bin);
	  data2[l][i] /= double(size_bin);
	  data4[l][i] /= double(size_bin);
	  data_abs[l][i] /= double(size_bin);
	  data_Binder[l][i] = 1.0 - 1.0 / 3.0 * data4[l][i] / data2[l][i] / data2[l][i];
	}
    }
}

double Bin_analysis::get_mean_value(int i){
    double result = Avg_vector(Num_bins, data[i]);
    return result;
}


double Bin_analysis::get_squared_value(int i){
      double result = Avg_vector(Num_bins, data2[i]);
      return result;
}

double Bin_analysis::get_abs_mean(int i){
      double result = Avg_vector(Num_bins, data_abs[i]);
      return result;
}
    
double Bin_analysis::get_Binder_cumulant(int i){
      double result = Avg_vector(Num_bins, data_Binder[i]);
      return result;
}    

void Bin_analysis::bootstrap(int i){
      //int Num_test = 1000;
      //data_bootstrap[i].resize(Num_tests, 0);
      //data_response_bootstrap[i].resize(Num_tests, 0);
      //data_Binder_bootstrap[i].resize(Num_tests, 0);
      
      for (int k = 0; k < Num_tests; k++){    
        for (int j = 0; j < Num_bins; j++){
            int random = int (uni01_sampler() * Num_bins);
            data_bootstrap[i][k] += data[i][random];
	    data_response_bootstrap[i][k] += data_response[i][random];
	    data_Binder_bootstrap[i][k] += data_Binder[i][random];
        }
        data_bootstrap[i][k] /= double(Num_bins);
	data_response_bootstrap[i][k] /= double(Num_bins);
	data_Binder_bootstrap[i][k] /= double(Num_bins);	
      }    
}

double Bin_analysis::get_mean_bootstrap(int i){
      return Avg_vector(Num_tests, data_bootstrap[i]);
}

double Bin_analysis::get_mean_error(int i){
      return Std_vector(Num_tests, data_bootstrap[i]);
}

double Bin_analysis::get_reponse_bootstrap(int i){
      return Avg_vector(Num_tests, data_response_bootstrap[i]);
}

double Bin_analysis::get_response_error(int i){
      return Std_vector(Num_tests, data_response_bootstrap[i]);
}

double Bin_analysis::get_Binder_bootstrap(int i){
      return Avg_vector(Num_tests, data_Binder_bootstrap[i]);
}

double Bin_analysis::get_Binder_error(int i){
      return Std_vector(Num_tests, data_Binder_bootstrap[i]);
}

// ------------------------------------------------------------------------------------------------------------------

// ------------------------------------------------------------------------------------------------------------------

Simple_bin::Simple_bin(string obs, int s, int i, INFO mc_info){
    set_parameters(obs, s, i, mc_info);
    input_open();
}

double Simple_bin::get_response_function(int i){
      double result = 0;
      return result;
}

Energy_bin::Energy_bin(string obs, int s, int i, INFO mc_info){
    set_parameters(obs, s, i ,mc_info);
    input_open();
}

double Energy_bin::get_response_function(int i){
      for (int j = 0; j < Num_bins; j++){
	data_response[i][j] = ( data2[i][j] - data[i][j] * data[i][j] ) / pow(temperature, 2) * double(Num_sites); 
      }
      double result = Avg_vector(Num_bins, data_response[i]);
      return result;
}


Magnetization_bin::Magnetization_bin(string obs, int s, int i, INFO mc_info){
    set_parameters(obs, s, i, mc_info);
    input_open();
}

double Magnetization_bin::get_response_function(int i){
      for (int j = 0; j < Num_bins; j++){
	data_response[i][j] = ( data2[i][j] - data_abs[i][j] * data_abs[i][j] ) / temperature * double(Num_sites); 
      }
      double result = Avg_vector(Num_bins, data_response[i]);
      return result; 
}


double Avg_vector(int num, vector<double> array){
    
    double summation = 0;
    for (int i = 0; i < num; i++){
        summation += array[i];
    }
    summation /= num;
    return summation;
}

double Var_vector(int num, vector<double> array){
    
    double variance = 0;
    for (int i = 0; i < num; i++){
        variance += array[i] * array[i];
    }
    variance /= num;
    return variance;
}

double Std_vector(int num, vector<double> array){

    double mean = Avg_vector(num, array);
    double standard_deviation = 0; 
    for (int i = 0; i < num; i++){
        standard_deviation += (array[i] - mean) * (array[i] - mean);
    }
    standard_deviation /= (num - 1.0);
    
    return sqrt(standard_deviation);           
}

