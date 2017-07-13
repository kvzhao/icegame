
#include "analysis.hpp"

Analysis::Analysis(string obs, int n, int m, INFO MC_info){
  
  set_parameters(obs, n, m, MC_info);
  
  switch(mode){
    case 0:
      simple_analysis(MC_info);
      output_open();
      output_write();
      output_close();
      break;
    case 1:
      energy_analysis(MC_info);
      output_open();
      output_write();
      output_close();
      break;
    case 2:
      magnetization_analysis(MC_info);
      output_open();
      output_write();
      output_close();
      break;
    case 3:
      correlation_analysis(MC_info);
      output_open();
      output_write_no_bootstrap();
      output_close();
      break;
  } 
}

void Analysis::set_parameters(string obs, int n, int m, INFO MC_info){

  name = obs;
  Num_observables = n;
  Num_replicas = MC_info.Num_replicas;
  mode = m;

  temperature.resize(Num_replicas);
  field.resize(Num_replicas);

  summ_mean.resize(Num_replicas);
  summ_error.resize(Num_replicas);
  summ_response.resize(Num_replicas);
  summ_Binder.resize(Num_replicas);

  summ_bootstrap.resize(Num_replicas);
  summ_bootstrap_error.resize(Num_replicas);
  summ_bootstrap_response.resize(Num_replicas);
  summ_bootstrap_response_error.resize(Num_replicas);
  summ_bootstrap_Binder.resize(Num_replicas);
  summ_bootstrap_Binder_error.resize(Num_replicas);
    
  for (int i = 0; i < Num_replicas; i++){
    summ_mean[i].resize(Num_observables, 0);
    summ_error[i].resize(Num_observables, 0);
    summ_response[i].resize(Num_observables, 0);
    summ_Binder[i].resize(Num_observables, 0);

    summ_bootstrap[i].resize(Num_observables, 0);
    summ_bootstrap_error[i].resize(Num_observables, 0);
    summ_bootstrap_response[i].resize(Num_observables, 0);
    summ_bootstrap_response_error[i].resize(Num_observables, 0);
    summ_bootstrap_Binder[i].resize(Num_observables, 0);
    summ_bootstrap_Binder_error[i].resize(Num_observables, 0);
  }
}

void Analysis::simple_analysis(INFO MC_info){

  for (int i = 0; i < Num_replicas; i++){  
    Simple_bin data_bin(name, Num_observables, i, MC_info);        
    data_bin.binning();
    temperature[i] = data_bin.get_temperature();
    //cout << "temp = " << temperature[i] << endl;
    field[i] = data_bin.get_field();
   
    for (int j = 0; j < Num_observables; j++){
      summ_mean[i][j] = data_bin.get_mean_value(j);
      //cout << " summ_mean" << summ_mean[i][j] << endl;
      summ_response[i][j] = data_bin.get_response_function(j);
      summ_Binder[i][j] = data_bin.get_Binder_cumulant(j);

      data_bin.bootstrap(j);
      summ_bootstrap[i][j] = data_bin.get_mean_bootstrap(j);
      summ_bootstrap_error[i][j] = data_bin.get_mean_error(j);
      summ_bootstrap_response[i][j] = data_bin.get_reponse_bootstrap(j);
      summ_bootstrap_response_error[i][j] = data_bin.get_response_error(j);
      summ_bootstrap_Binder[i][j] = data_bin.get_Binder_bootstrap(j);
      summ_bootstrap_Binder_error[i][j] = data_bin.get_Binder_error(j);
    }
  data_bin.input_close();  
  }
}

void Analysis::correlation_analysis(INFO MC_info){

  for (int i = 0; i < Num_replicas; i++){  
    Simple_bin data_bin(name, Num_observables, i, MC_info);        
    data_bin.binning();
    temperature[i] = data_bin.get_temperature();
    field[i] = data_bin.get_field();
   
    for (int j = 0; j < Num_observables; j++){
      summ_mean[i][j] = data_bin.get_mean_value(j);
      summ_response[i][j] = data_bin.get_response_function(j);
      summ_Binder[i][j] = data_bin.get_Binder_cumulant(j);

    }
  data_bin.input_close();  
  }
}



void Analysis::energy_analysis(INFO MC_info){

  for (int i = 0; i < Num_replicas; i++){  
    Energy_bin data_bin(name, Num_observables, i, MC_info);        
    data_bin.binning();
    temperature[i] = data_bin.get_temperature();
    field[i] = data_bin.get_field();
   
    for (int j = 0; j < Num_observables; j++){
      summ_mean[i][j] = data_bin.get_mean_value(j);
      summ_response[i][j] = data_bin.get_response_function(j);
      summ_Binder[i][j] = data_bin.get_Binder_cumulant(j);

      data_bin.bootstrap(j);
      summ_bootstrap[i][j] = data_bin.get_mean_bootstrap(j);
      summ_bootstrap_error[i][j] = data_bin.get_mean_error(j);
      summ_bootstrap_response[i][j] = data_bin.get_reponse_bootstrap(j);
      summ_bootstrap_response_error[i][j] = data_bin.get_response_error(j);
      summ_bootstrap_Binder[i][j] = data_bin.get_Binder_bootstrap(j);
      summ_bootstrap_Binder_error[i][j] = data_bin.get_Binder_error(j);
    }
   data_bin.input_close(); 
  }
}

void Analysis::magnetization_analysis(INFO MC_info){

  for (int i = 0; i < Num_replicas; i++){  
    Magnetization_bin data_bin(name, Num_observables, i, MC_info);        
    data_bin.binning();
    temperature[i] = data_bin.get_temperature();
    field[i] = data_bin.get_field();
   
    for (int j = 0; j < Num_observables; j++){
      // ************************************************************************************************
      //summ_mean[i][j] = data_bin.get_abs_mean(j);
      summ_mean[i][j] = data_bin.get_mean_value(j);
      // ************************************************************************************************
      summ_response[i][j] = data_bin.get_response_function(j);
      summ_Binder[i][j] = data_bin.get_Binder_cumulant(j);

      data_bin.bootstrap(j);
      summ_bootstrap[i][j] = data_bin.get_mean_bootstrap(j);
      summ_bootstrap_error[i][j] = data_bin.get_mean_error(j);
      summ_bootstrap_response[i][j] = data_bin.get_reponse_bootstrap(j);
      summ_bootstrap_response_error[i][j] = data_bin.get_response_error(j);
      summ_bootstrap_Binder[i][j] = data_bin.get_Binder_bootstrap(j);
      summ_bootstrap_Binder_error[i][j] = data_bin.get_Binder_error(j);
    }
   data_bin.input_close(); 
  }
}




Analysis::~Analysis(){
  summ_mean.clear();
  summ_error.clear();
  summ_response.clear();
  summ_bootstrap.clear();
}
    
void Analysis::output_open(){
  ostringstream filename;
  filename << "Summary_" << name;
  output.open( filename.str().c_str() );
} 

void Analysis::output_close(){ output.close(); }

void Analysis::output_write(){
  
  output << "index" << '\t'
         << "temperature" << '\t'
         << "field" << '\t';
  for (int i = 0; i < Num_observables; i++){
    output << "mean_" << i << '\t'
           << "response_" << i << '\t'
           << "Binder_" << i << '\t'
           << "bt_mean_" << i << '\t'
           << "bt_error_" << i << '\t'
           << "bt_resp_" << i << '\t'
           << "bt_resp_err_" << i << '\t'
           << "bt_Binder_" << i << '\t'
           << "bt_Binder_err" << i << '\t';
  }
  output << '\n'; 
  

  for (int i = 0; i < Num_replicas; i++){
    output << i << '\t'
           << temperature[i] << '\t'
           << field[i] << '\t';
    for (int j = 0; j < Num_observables; j++){
      output << summ_mean[i][j] << '\t'
             << summ_response[i][j] << '\t'
             << summ_Binder[i][j] << '\t'
             << summ_bootstrap[i][j] << '\t'
             << summ_bootstrap_error[i][j] << '\t'
             << summ_bootstrap_response[i][j] << '\t'
             << summ_bootstrap_response_error[i][j] << '\t'
             << summ_bootstrap_Binder[i][j] << '\t'
             << summ_bootstrap_Binder_error[i][j] << '\t';
    }
    output << '\n';
  }

}

void Analysis::output_write_no_bootstrap(){
  
  output << "index" << '\t'
         << "temperature" << '\t'
         << "field" << '\t';
  for (int i = 0; i < Num_observables; i++){
    output << "mean_" << i << '\t';
  }
  output << '\n'; 
  

  for (int i = 0; i < Num_replicas; i++){
    output << i << '\t'
           << temperature[i] << '\t'
           << field[i] << '\t';
    for (int j = 0; j < Num_observables; j++){
      output << summ_mean[i][j] << '\t';
    }
    output << '\n';
  }

}
