// Try to construct the simplest MC simulation code with OOP
// initialization -> update -> measure -> output -> binning/bootstrap
// measure energy only
// without magnetic field

#include "hamiltonian.hpp"


// -----------------------------------------------------------------------------------------------------------
// Ising model with 2nd nearest-neighbor interaction
// update algorithm: SSF, SSF2
// measured quantity: energy, magnetization, stagger_magnetization(square), acceptance
// 		      energy_update, magnetization_update
// initial state: random configuration, ferromagnetic orderded state 
// -----------------------------------------------------------------------------------------------------------

double Hamiltonian::compare_state(Sample* s1, Sample* s2){
  double result = 0;
  for (int i = 0; i < N; i++){
    if (s1->Ising[i] == s2->Ising[i])
      result += 1;
  }
  result /= double(N);
  
  return result;
}

Ising::Ising(INFO MC, int nei){
      L = MC.lattice_size;
      N = MC.Num_sites;
      R = MC.Num_replicas;
      Num_neighbors = nei;
      J1 = -1.0; //default value, can be changed by other functions
      J2 = 0;
      J3a = 0;
      J3b = 0;
      cout << " ------- Hamiltonian: Ising model with layer of neighbors = " << nei << endl;
}
    
Ising::~Ising(){
}

void Ising::set_J1(double J){
  J1 = J;
  cout << " ------- Set parameter J1 = " << J1 << endl;
}

void Ising::set_J2(double J){
  J2 = J;
  if (Num_neighbors >1)
    cout << " ------- Set parameter J2 = " << J2 << endl;
}

void Ising::set_J3a(double J){
  J3a = J;
  if (Num_neighbors > 2)
    cout << " ------- Set parameter J3a = " << J3a << endl;
}

void Ising::set_J3b(double J){
  J3b = J;
  if (Num_neighbors > 2)
    cout << " ------- Set parameter J3b = " << J3b << endl;
}

double Ising::SE(int p, Sample* ss, Lattice* la){
      double se = 0;
      vector<int>::iterator NN_itr;
      for (NN_itr = la->NN[p].begin(); NN_itr != la->NN[p].end(); NN_itr++)
	se += ss->Ising[p] * ss->Ising[(*NN_itr)] * J1;
      return se;
}
 
double Ising::SE2(int p, Sample* ss, Lattice* la){
      double se = 0;
      vector<int>::iterator NNN_itr;
      for (NNN_itr = la->NNN[p].begin(); NNN_itr != la->NNN[p].end(); NNN_itr++)
	se += ss->Ising[p] * ss->Ising[(*NNN_itr)] * J2;
      return se;
}

double Ising::SE3(int p, Sample* ss, Lattice* la){
      double se = 0;
      vector<int>::iterator NNNNa_itr;
      vector<int>::iterator NNNNb_itr;
      for (NNNNa_itr = la->NNNN[p].begin(); NNNNa_itr != la->NNNN[p].end(); NNNNa_itr++)
        se += ss->Ising[p] * ss->Ising[(*NNNNa_itr)] * J3a;
      for (NNNNb_itr = la->NNNN2[p].begin(); NNNNb_itr != la->NNNN2[p].end(); NNNNb_itr++)
        se += ss->Ising[p] * ss->Ising[(*NNNNb_itr)] * J3b;
      return se;
}      

double Ising::ME(int p, Sample* ss){
      double me = (-1.0) * ss->Ising[p] * ss->get_field(0);
      return me;
}

void Ising::MCstep(Sample* ss, Lattice* la){
      
      double acceptance = 0;
      
      switch(Num_neighbors){

	case 1:{
	    for (int j = 0 ; j < N; j++){
	      if ( SSF(ss, la) )
	      acceptance++;
	    }
	    ss->acceptance_SSF = acceptance;
	    break;
	}
	case 2:{
	    for (int j = 0 ; j < N; j++){
	      if ( SSF2(ss, la) )
	      acceptance++;
	    }
	    ss->acceptance_SSF = acceptance;
	    break;
	}
      }
}

bool Ising::SSF3(Sample* ss, Lattice* la){
      int random = (uni01_sampler() * N);
      double delta_m = (-2.0) * ss->Ising[random]; 
      double delta_se = (-2.0) * SE(random, ss, la) + (-2.0) * SE2(random, ss, la) + (-2.0) * SE3(random, ss, la);
      double delta_me = (-2.0) * ME(random, ss);
      double delta = delta_se + delta_me;
      
      double weight = exp(-delta / ss->get_temperature() );
      double dice = uni01_sampler();
      
      if (dice < weight){
	ss->Ising[random] *= -1;
	ss->energy_update += delta;
	ss->magnetization_update[0] += delta_m;
	return true;
      }
      else{
	return false;
      }
 
}


bool Ising::SSF2(Sample* ss, Lattice* la){
      int random = (uni01_sampler() * N);
      double delta_m = (-2.0) * ss->Ising[random];
      double delta_se = (-2.0) * SE(random, ss, la) + (-2.0) * SE2(random, ss, la);
      double delta_me = (-2.0) * ME(random, ss);
      double delta = delta_se + delta_me;
      
      double weight = exp(-delta / ss->get_temperature() );
      double dice = uni01_sampler();
      
      if (dice < weight){
	ss->Ising[random] *= -1;
	ss->energy_update += delta;
	ss->magnetization_update[0] += delta_m;
	return true;
      }
      else{
	return false;
      }
 
}

bool Ising::SSF(Sample* ss, Lattice* la){
      int random = (uni01_sampler() * N);
      double delta_m = (-2.0) * ss->Ising[random];
      double delta_se = (-2.0) * SE(random, ss, la);
      double delta_me = (-2.0) * ME(random, ss);
      double delta = delta_se + delta_me;
      
      double weight = exp(-delta / ss->get_temperature() );
      double dice = uni01_sampler();
      
      if (dice < weight){
	ss->Ising[random] *= -1;
	ss->energy_update += delta;
	ss->magnetization_update[0] += delta_m;
	return true;
      }
      else{
	return false;
      }      
}

double Ising::total_energy(Sample* ss, Lattice* la){
  double eng = 0;
  
  switch(Num_neighbors){
    case 1:{
      for (int i = 0; i < N; i++)
	eng += SE(i,ss,la) / double(2.0) + ME(i,ss);
      break;
    }
    case 2:{
      for (int i = 0; i < N; i++)
	eng += SE(i, ss, la) / double(2.0) + SE2(i, ss, la) / double(2.0) + ME(i, ss);
      break;
    }
    case 3:{
      for (int i = 0; i < N; i++)
	eng += SE(i,ss,la) / double(2.0) + SE2(i,ss,la) / double(2.0) + SE3(i,ss,la) / double(2.0) + ME(i,ss);
      break;
    }
  }

  return eng;
} 

double Ising::magnetization(Sample* ss, Lattice* la){
  double mag = 0;
  for (int i = 0; i < N; i++)
    mag += ss->Ising[i];
  return mag;
}

double Ising::stagger_magnetization(Sample* ss, Lattice* la){
  double mag = 0;
  for (int i = 0; i < ss->N; i++)
    mag += ss->Ising[i] * la->sub[i];
  mag /= double(ss->N);
  return mag;
}



double Ising::spin_correlation(Sample* ss, Lattice* la,vector<double> &data){
  
  int LL = la->L / 2;
  int NN = la->N;

  data.clear();
  data.resize(LL);
  
  for (int i = 0; i < LL; i++){
    double avg_sc = 0;
    for (int j = 0; j < NN; j++){
    double sc = 0;
    int xx = la->site2d[j].x;
    int yy = la->site2d[j].y;
    //int xp = xx + i;
    //int yp = yy + i;
    //int xm = xx - i;
    //int ym = yy - i;
    
    int xp = la->periodic_bd(xx, i);
    int yp = la->periodic_bd(yy, i);
    int xm = la->periodic_bd(xx, -i);
    int ym = la->periodic_bd(yy, -i);
/*
    if (xp >= LL)
      xp = xp % LL;
    if (yp >= LL)
      yp = yp % LL;
    if (xm < 0)
      xm = (xm % LL + LL) % LL;
    if (ym < 0)
      ym = (ym % LL + LL) % LL;
*/    
    //cout << "xm = " << xm << ", ym = " << ym << endl;

    int p1 = la->site_index2d[xp][yy];
    int p2 = la->site_index2d[xx][yp];
    int p3 = la->site_index2d[xm][yy];
    int p4 = la->site_index2d[xx][ym]; 
    //cout << "p1 = " << model->site[p1].x << ", " << model->site[p1].y << endl;
    //cout << "p2 = " << model->site[p2].x << ", " << model->site[p2].y << endl;
    //cout << "p3 = " << model->site[p3].x << ", " << model->site[p3].y << endl;
    //cout << "p4 = " << model->site[p4].x << ", " << model->site[p4].y << endl; 
    sc += ss->Ising[j] * (ss->Ising[p1] + ss->Ising[p2] + ss->Ising[p3] + ss->Ising[p4]);
    sc /= double(4);
    avg_sc += sc;
    }
    data[i] = avg_sc / double(NN);
    //ss->correlation[i] += sc;
  }
  return 0;
}


double Ising::obs_energy(Sample* ss, Lattice* la, vector<double> &data){
  data.clear();
  data.resize(1);
  double e = ss->energy_update;
  e /= double(ss->N);
  data[0] = e;
  return e;
}

double Ising::obs_magnetization(Sample* ss, Lattice* la, vector<double> &data){
  data.clear();
  data.resize(2);
  double m = ss->magnetization_update[0];
  m /= double(ss->N);
  data[0] = m;
  data[1] = abs(m);
  return m;
}

double Ising::obs_acceptanceSSF(Sample* ss, Lattice* la, vector<double> &data){
  data.clear();
  data.resize(1);
  double a = ss->acceptance_SSF;
  a /= double(ss->N);
  data[0] = a;
  return a;
}


void Ising::initialization(Sample* ss, Lattice* la, int initial_state){
  
  switch(initial_state){
    case 1:  // all spins up
	     for(int i = 0; i < N; i++)
	       ss->Ising[i] = 1;
	     break;
    case 0:  // random initial state
	     for (int i = 0; i < N; i++)
	       ss->Ising[i] = 2 * int( uni01_sampler() * 2) - 1;
	   
  }
}

void Ising::initialize_observable(Sample* ss, Lattice* la, double t, vector<double> h){
  ss->set_temperature(t);
  ss->resize_field(1);
  ss->set_field(h[0], 0);

  ss->resize_magnetization(1);
  ss->resize_order(0);

  ss->set_energy(total_energy(ss, la));
  ss->set_magnetization(0, magnetization(ss, la));

}

// -----------------------------------------------------------------------------------------------------------
// 3D spin ice model with 2nd nearest-neighbor interaction
// update algorithm: SSF, SSF2
// measured quantity: energy, magnetization, stagger_magnetization(square), acceptance
// 		      energy_update, magnetization_update
// initial state: random configuration, ferromagnetic orderded state 
// -----------------------------------------------------------------------------------------------------------

Kagome_ice::Kagome_ice(INFO MC, int nei){
      L = MC.lattice_size;
      N = MC.Num_sites;
      R = MC.Num_replicas;
      Num_neighbors = nei;
      J1 = -1.0; //default value, can be changed by other functions
      J2 = 0;
      J3a = 0;
      J3b = 0;
      //easyaxis = new spin_easyaxis;
      //set_easyaxis();
      
     cout << " ------- Hamiltonian: Spin ice model with layer of neighbors = " << nei << endl;
}
    
Kagome_ice::~Kagome_ice(){
}

void Kagome_ice::get_easyaxis(int a, vector<double> &vtr){
  double sq2 = sqrt(2);
  double sq3 = sqrt(3);
  switch(a){
    case 0:
      vtr[0] = -1.0 * sq2 / sq3;
      vtr[1] = -1.0 * sq2 / 3.0;
      vtr[2] = 1.0 / 3.0;
      break;
    case 1:
      vtr[0] = sq2 / sq3;
      vtr[1] = -1.0 * sq2 / 3.0;
      vtr[2] = 1.0 / 3.0;
      break;
    case 2:
      vtr[0] = 0;
      vtr[1] = 2.0 * sq2 / 3.0;
      vtr[2] = 1.0 / 3.0;
      break;
  }
} 
/*
Kagome_ice::spin_easyaxis::spin_easyaxis(){
      sq2 = sqrt(2);
      sq3 = sqrt(3);
      sublattice.resize(3);
      for (int i = 0; i < 3; i++)
       sublattice[i].resize(3);	

      sublattice[0][0] = -1.0 * sq2 / sq3;
      sublattice[0][1] = -1.0 * sq2 / 3.0;
      sublattice[0][2] = 1.0 / 3.0;
      sublattice[1][0] = sq2 / sq3;
      sublattice[1][1] = -1.0 * sq2 / 3.0;
      sublattice[1][2] = 1.0 / 3.0;
      sublattice[2][0] = 0;
      sublattice[2][1] = 2.0 * sq2 / 3.0;
      sublattice[2][2] = 1.0 / 3.0;
}

static void Kagome_ice::set_easyaxis(){
  double sq2 = sqrt(2);
  double sq3 = sqrt(3);

  easyaxis[0].x = -1.0 * sq2 / sq3;
  easyaxis[0].y = -1.0 * sq2 / 3.0;
  easyaxis[0].z = 1.0 / 3.0;
  easyaxis[1].x = sq2 / sq3;
  easyaxis[1].y = -1.0 * sq2 / 3.0;
  easyaxis[1].z = 1.0 / 3.0;
  easyaxis[2].x = 0;
  easyaxis[2].y = 2.0 * sq2 / 3.0;
  easyaxis[2].z = 1.0 / 3.0;

}*/

void Kagome_ice::set_J1(double J){
  J1 = J;
  //cout << " ------- Set parameter J1 = " << J1 << endl;
}

void Kagome_ice::set_J2(double J){
  J2 = J;
  //if (Num_neighbors > 1)
    //cout << " ------- Set parameter J2 = " << J2 << endl;
}

void Kagome_ice::set_J3a(double J){
  J3a = J;
  //if (Num_neighbors > 2)
    //cout << " ------- Set parameter J3a = " << J3a << endl;
}

void Kagome_ice::set_J3b(double J){
  J3b = J;
  //if (Num_neighbors > 2)
    //cout << " ------- Set parameter J3b = " << J3b << endl;
}

void Kagome_ice::set_product_hs(Sample* ss, double h1, double h2, double h3){
    
    ss->product_hs.resize(3);
    // 3-dim vector of magnetic field, {x, y, z} = {[111], [-1-12], [-110]}
    double h_vector[3] = {h1, h2, h3};
    
    // easy axis in the kagome coordinates, {x, y, z} = {[111], [-1-12], [-110]}
    double s0[3] = {-1.0/3.0, -sqrt(2)/3.0, -sqrt(2.0/3.0)};
    double s1[3] = {-1.0/3.0, -sqrt(2)/3.0, sqrt(2.0/3.0)};
    double s2[3] = {-1.0/3.0 , 2.0*sqrt(2)/3.0, 0};
    
    // product of h ans s for the four sublattices
    ss->product_hs[0] = std::inner_product(s0, s0+3, h_vector, 0.0);
    ss->product_hs[1] = std::inner_product(s1, s1+3, h_vector, 0.0);
    ss->product_hs[2] = std::inner_product(s2, s2+3, h_vector, 0.0);
}


double Kagome_ice::SE(int p, Sample* ss, Lattice* la){
      double se = 0;
      vector<int>::iterator NN_itr;
      for (NN_itr = la->NN[p].begin(); NN_itr != la->NN[p].end(); NN_itr++)
	se += ss->Ising[p] * ss->Ising[(*NN_itr)] * J1;
      return se;
}
 
double Kagome_ice::SE2(int p, Sample* ss, Lattice* la){
      double se = 0;
      vector<int>::iterator NNN_itr;
      for (NNN_itr = la->NNN[p].begin(); NNN_itr != la->NNN[p].end(); NNN_itr++)
	se += ss->Ising[p] * ss->Ising[(*NNN_itr)] * J2;
      return se;
}

double Kagome_ice::SE3(int p, Sample* ss, Lattice* la){
      double se = 0;
      vector<int>::iterator NNNNa_itr;
      vector<int>::iterator NNNNb_itr;
      for (NNNNa_itr = la->NNNN[p].begin(); NNNNa_itr != la->NNNN[p].end(); NNNNa_itr++)
        se += ss->Ising[p] * ss->Ising[(*NNNNa_itr)] * J3a;
      for (NNNNb_itr = la->NNNN2[p].begin(); NNNNb_itr != la->NNNN2[p].end(); NNNNb_itr++)
        se += ss->Ising[p] * ss->Ising[(*NNNNb_itr)] * J3b;
      return se;
}      


double Kagome_ice::ME(int p, Sample* ss, Lattice* la){

    double result = 0;

  switch(la->sub[p]){
    case 0: 
	     result = (-1.0) * ss->product_hs[0] * ss->Ising[p];
	     break;
    case 1:
	     result = (-1.0) * ss->product_hs[1] * ss->Ising[p];
	     break;
    case 2:
	     result = (-1.0) * ss->product_hs[2] * ss->Ising[p];
	     break;
  }
  return result;
}


double Kagome_ice::M_111(int p, Sample* ss, Lattice* la){ //2d kagome ice
	double result;
	result = (-1.0) * (1.0/3.0) * ss->Ising[p];
	return result;
}

double Kagome_ice::M_112(int p, Sample* ss, Lattice* la){ //2d kagome ice
    double result;
    switch(la->sub[p]){
        case 0:
            result = (-1.0) * ss->Ising[p] * sqrt(2)/3.0;
            break; 
        case 1:
            result = (-1.0) * ss->Ising[p] * sqrt(2)/3.0;
            break;
        case 2:
            result = ss->Ising[p] * 2.0 * sqrt(2)/3.0;
            break; 
    }
    return result;
}

double Kagome_ice::M_110(int p, Sample* ss, Lattice* la){ //2d kagome ice
    double result;
    switch(la->sub[p]){
        case 0:
            result = (-1.0) * ss->Ising[p] * sqrt(2.0/3.0);
            break;
        case 1:
            result = ss->Ising[p] * sqrt(2.0/3.0);
            break;
        case 2:
            result = 0;
            break;
    }
    return result;
}

double Kagome_ice::magnetization_111(Sample* ss, Lattice* la){ //2d kagome
    
    double sum_M = 0;
    for (int i = 0; i < N; i++)
        sum_M += M_111(i, ss, la);
    return (sum_M);
}

double Kagome_ice::magnetization_112(Sample* ss, Lattice* la){ //2d kagome
    
    double sum_M112 = 0;
    for (int i = 0; i < N; i++)
      sum_M112 += M_112(i, ss, la);    
    return (sum_M112);
}

double Kagome_ice::magnetization_110(Sample* ss, Lattice* la){ //2d kagome

    double sum_M110 = 0;
    for (int i = 0; i < N; i++)
      sum_M110 += M_110(i, ss, la);
    return (sum_M110);
}


void Kagome_ice::MCstep(Sample* ss, Lattice* la){
      
      double acceptance = 0;
      
      switch(Num_neighbors){

	case 1:
	    for (int j = 0 ; j < N; j++){
	      if ( SSF(ss, la) )
	      acceptance++;
	    }
	    ss->acceptance_SSF = acceptance;
	    break;
	case 2:
	    for (int j = 0 ; j < N; j++){
	      if ( SSF2(ss, la) )
	      acceptance++;
	    }
	    ss->acceptance_SSF = acceptance;
	    break;
        case 3:
	    for (int j = 0; j < N; j++){
	     if (SSF3(ss, la) ) 
	       acceptance++;
	    }
	    ss->acceptance_SSF = acceptance;
	    break;
      }
}

bool Kagome_ice::SSF3(Sample* ss, Lattice* la){
      int random = (uni01_sampler() * N);
      double delta_m111 = (-2.0) * M_111(random , ss, la);
      double delta_m112 = (-2.0) * M_112(random, ss, la);
      double delta_m110 = (-2.0) * M_110(random, ss, la);

      double delta_se = (-2.0) * SE(random, ss, la) + (-2.0) * SE2(random, ss, la) + (-2.0) * SE3(random, ss, la);
      double delta_me = (-2.0) * ME(random, ss, la);
      double delta = delta_se + delta_me;
      
      double weight = exp(-delta / ss->get_temperature() );
      double dice = uni01_sampler();
      
      if (dice < weight){
	ss->Ising[random] *= -1;
	ss->energy_update += delta;
	ss->magnetization_update[0] += delta_m111;
	ss->magnetization_update[1] += delta_m112;
	ss->magnetization_update[2] += delta_m110;
	return true;
      }
      else{
	return false;
      }
 
}


bool Kagome_ice::SSF2(Sample* ss, Lattice* la){
      int random = (uni01_sampler() * N);
      double delta_m111 = (-2.0) * M_111(random , ss, la);
      double delta_m112 = (-2.0) * M_112(random, ss, la);
      double delta_m110 = (-2.0) * M_110(random, ss, la);

      double delta_se = (-2.0) * SE(random, ss, la) + (-2.0) * SE2(random, ss, la);
      double delta_me = (-2.0) * ME(random, ss, la);
      double delta = delta_se + delta_me;
      
      double weight = exp(-delta / ss->get_temperature() );
      double dice = uni01_sampler();
      
      if (dice < weight){
	ss->Ising[random] *= -1;
	ss->energy_update += delta;
	ss->magnetization_update[0] += delta_m111;
	ss->magnetization_update[1] += delta_m112;
	ss->magnetization_update[2] += delta_m110;
	return true;
      }
      else{
	return false;
      }
 
}

bool Kagome_ice::SSF(Sample* ss, Lattice* la){
      int random = (uni01_sampler() * N);
      double delta_m111 = (-2.0) * M_111(random, ss, la);
      double delta_m112 = (-2.0) * M_112(random, ss, la);
      double delta_m110 = (-2.0) * M_110(random, ss, la);
       
      double delta_se = (-2.0) * SE(random, ss, la);
      double delta_me = (-2.0) * ME(random, ss, la);
      double delta = delta_se + delta_me;
      
      double weight = exp(-delta / ss->get_temperature() );
      double dice = uni01_sampler();
      
      if (dice < weight){
	ss->Ising[random] *= -1;
	ss->energy_update += delta;
	ss->magnetization_update[0] += delta_m111;
	ss->magnetization_update[1] += delta_m112;
	ss->magnetization_update[2] += delta_m110;
	return true;
      }
      else{
	return false;
      }      
}


double Kagome_ice::total_energy(Sample* ss, Lattice* la){
  double eng = 0;
  
  switch(Num_neighbors){
    case 1:
      for (int i = 0; i < N; i++)
	eng += SE(i,ss,la) / double(2.0) + ME(i,ss, la);
      break;
    case 2:
      for (int i = 0; i < N; i++)
	eng += SE(i, ss, la) / double(2.0) + SE2(i, ss, la) / double(2.0) + ME(i, ss, la);
      break;
    case 3:
      for (int i = 0; i < N; i++)
	eng += SE(i, ss, la) / double(2.0) + SE2(i, ss, la) / double(2.0) + SE3(i, ss, la) / double(2.0) + ME(i, ss, la);
      break;
  }
  return eng;
} 

double Kagome_ice::spin_correlation(Sample* ss, Lattice* la,vector<double> &sc_list){
  
  /*
  for (int i = 0; i < L; i++){
    double sc = 0;
    int xp = i;
    int yp = i;
    int xm = -i;
    int ym = -i;

    if (xp >= L)
      xp = xp % L;
    if (yp >= L)
      yp = yp % L;
    if (xm < 0)
      xm = (xm % L + L) % L;
    if (ym < 0)
      ym = (ym % L + L) % L;
    
    //cout << "xm = " << xm << ", ym = " << ym << endl;

    int p1 = la->site_index2[xp][0];
    int p2 = la->site_index2[0][yp];
    int p3 = la->site_index2[xm][0];
    int p4 = la->site_index2[0][ym]; 
  
    //cout << "p1 = " << model->site[p1].x << ", " << model->site[p1].y << endl;
    //cout << "p2 = " << model->site[p2].x << ", " << model->site[p2].y << endl;
    //cout << "p3 = " << model->site[p3].x << ", " << model->site[p3].y << endl;
    //cout << "p4 = " << model->site[p4].x << ", " << model->site[p4].y << endl; 
    sc += ss->Ising[0] * (ss->Ising[p1] + ss->Ising[p2] + ss->Ising[p3] + ss->Ising[p4]);
    sc /= double(4);
    sc_list[i] = sc;
    ss->correlation[i] += sc;
  }*/
  return 0;
}


double Kagome_ice::obs_energy(Sample* ss, Lattice* la, vector<double> &data){
  data.clear();
  data.resize(1);
  double e = ss->energy_update;
  e /= double(ss->N);
  data[0] = e;
  return 0;
}

double Kagome_ice::obs_magnetization111(Sample* ss, Lattice* la){
  double m = ss->magnetization_update[0];
  m /= double(ss->N);
  return m;
}


double Kagome_ice::obs_magnetization112(Sample* ss, Lattice* la){
  double m = ss->magnetization_update[1];
  m /= double(ss->N);
  return m;
}


double Kagome_ice::obs_magnetization110(Sample* ss, Lattice* la){
  double m = ss->magnetization_update[2];
  m /= double(ss->N);
  return m;
}


double Kagome_ice::obs_acceptanceSSF(Sample* ss, Lattice* la, vector<double> &data){
  data.clear();
  data.resize(1);
  double a = ss->acceptance_SSF;
  a /= double(ss->N);
  data[0] = a;
  return 0;
}
/*
void Kagome_ice::make_configuration(vector<Sample*> &vs, INFO MC_info){

  vs.resize(9);
  vector<Sample*>::iterator iter;
  for (iter = vs.begin(); iter != vs.end(); iter++){
    (*iter) = new Sample(MC_info);
  }
  

}
*/
void Kagome_ice::monopole_vortex(Sample* ss, Lattice* la, int n){

  switch(n){
    case 0:
      degenerate_vortex(ss, la, 0);
      break;
    case 1:
      degenerate_vortex(ss, la, 1);
      break;
    case 2:
      degenerate_vortex(ss, la, 2);
      break;
  }
  for (int i = 0; i < N; i++)
    if (la->sub[i] == 2)
      ss->Ising[i] = -1;
}

void Kagome_ice::anti_vortex(Sample* ss, Lattice* la, int n){

  switch(n){
    case 0:
      degenerate_vortex(ss, la, 0);
      break;
    case 1:
      degenerate_vortex(ss, la, 1);
      break;
    case 2:
      degenerate_vortex(ss, la, 2);
      break;
  }
  for (int i = 0; i < N; i++)
    ss->Ising[i] *= -1;
}

void Kagome_ice::degenerate_vortex(Sample* ss, Lattice* la, int n){

  initialization(ss, la, 1);
  switch(n){
    case 0:
      break;
    case 1:
      for (int i = 0; i < L; i++){
	for (int j = 0; j < L; j++){
	  int n0 = ss->Ising[la->site_index3d[i][j][0]];
	  int n1 = ss->Ising[la->site_index3d[i][j][1]];
	  int n2 = ss->Ising[la->site_index3d[i][j][2]];
	  ss->Ising[la->site_index3d[i][j][0]] = n1;
	  ss->Ising[la->site_index3d[i][j][1]] = n2;
	  ss->Ising[la->site_index3d[i][j][2]] = n0;
	}
      }
      break;
    case 2:
      for (int i = 0; i < L; i++){
	for (int j = 0; j < L; j++){
	  int n0 = ss->Ising[la->site_index3d[i][j][0]];
	  int n1 = ss->Ising[la->site_index3d[i][j][1]];
	  int n2 = ss->Ising[la->site_index3d[i][j][2]];
	  ss->Ising[la->site_index3d[i][j][0]] = n2;
	  ss->Ising[la->site_index3d[i][j][1]] = n0;
	  ss->Ising[la->site_index3d[i][j][2]] = n1;
	}
      }
      break;

  }
}

void Kagome_ice::degenerate_qx(Sample* ss, Lattice* la, int n){

  initialization(ss, la, 2);
  switch(n){
    case 0:
      break;
    case 1:
      for (int i = 0; i < N; i++){
	if (la->sub[i] != 2)
	  ss->Ising[i] *= -1;
      }
      break;
  }

}

void Kagome_ice::initialization(Sample* ss, Lattice* la, int initial_state){
  
  switch(initial_state){
    case 1:  // vortex state  
	    for(int j = 0; j < L; j++){ // y-axis
	      for(int i = 0; i < L; i++){	// x-axis
		if (j == 0){
		  int key = i % 3;
		  switch(key){
		    case 0:
		      ss->Ising[la->site_index3d[i][j][0]] = 1;
		      ss->Ising[la->site_index3d[i][j][1]] = -1;
		      ss->Ising[la->site_index3d[i][j][2]] = -1;
		      break;
		    case 1:
		      ss->Ising[la->site_index3d[i][j][0]] = -1;
		      ss->Ising[la->site_index3d[i][j][1]] = 1;
		      ss->Ising[la->site_index3d[i][j][2]] = -1;
		      break;
		    case 2:
		      ss->Ising[la->site_index3d[i][j][0]] = -1;
		      ss->Ising[la->site_index3d[i][j][1]] = -1;
		      ss->Ising[la->site_index3d[i][j][2]] = 1;
		      break;
		  }
	       }
	      else{
		if (i == 0){
		  ss->Ising[la->site_index3d[i][j][0]] = ss->Ising[la->site_index3d[L - 1][j - 1][0]];
		  ss->Ising[la->site_index3d[i][j][1]] = ss->Ising[la->site_index3d[L - 1][j - 1][1]];
		  ss->Ising[la->site_index3d[i][j][2]] = ss->Ising[la->site_index3d[L - 1][j - 1][2]];
		}
		else{
		  ss->Ising[la->site_index3d[i][j][0]] = ss->Ising[la->site_index3d[i - 1][j - 1][0]];
		  ss->Ising[la->site_index3d[i][j][1]] = ss->Ising[la->site_index3d[i - 1][j - 1][1]];
		  ss->Ising[la->site_index3d[i][j][2]] = ss->Ising[la->site_index3d[i - 1][j - 1][2]];
		}
	      }
	      }
	    }
	    break;

    case 0:  // random initial state
	     for (int i = 0; i < N; i++)
	       ss->Ising[i] = 2 * int( uni01_sampler() * 2) - 1;
	     break;
    case 2: // q = X state
	     for (int i = 0; i < N; i++){
		if (la->sub[i] == 2)
		  ss->Ising[i] = -1;
		else{
		  int key = (la->site3d[i].y + la->sub[i]) % 2;
		  switch(key){
		    case 0:
		      ss->Ising[i] = 1;
		      break;
		    case 1:
		      ss->Ising[i] = -1;
		      break;
		  }
		}
	     }
	     break;
    case 3: // saturated state
	     for (int i = 0; i < N; i++){
	      ss->Ising[i] = -1;
	     }
	    break; 
    case 4: // two q = X domain state
	     for (int i = 0; i < N; i++){
		if (la->sub[i] == 2)
		  ss->Ising[i] = -1;
		else{
		  int key = (la->site3d[i].y + la->sub[i]) % 2;
		  int boundary = L / 2;
		  if (la->site3d[i].x >= boundary){
		    switch(key){
		      case 0:
			ss->Ising[i] = -1;
			break;
		      case 1:
			ss->Ising[i] = 1;
			break;
		    }
		  }
		  else{
		    switch(key){
		      case 0:
			ss->Ising[i] = 1;
			break;
		      case 1:
			ss->Ising[i] = -1;
			break;
		    }
		  }
		}
	     }
	     break;
	   
  }
}


void Kagome_ice::initialize_observable(Sample* ss, Lattice* la, double t, vector<double> h){
  
  ss->set_temperature(t);
  ss->resize_field(3);  
  ss->set_field(h[0], 0);
  ss->set_field(h[1], 1);
  ss->set_field(h[2], 2);
  
  set_product_hs(ss, h[0], h[1], h[2]);

  ss->resize_magnetization(3);
  //ss->resize_order(4); // for 3 fully-polarized order

  ss->set_energy(total_energy(ss, la));
  ss->set_magnetization(0, magnetization_111(ss, la));
  ss->set_magnetization(1, magnetization_112(ss, la));
  ss->set_magnetization(2, magnetization_110(ss, la));
}

double Kagome_ice::defect_density(Sample* ss, Lattice* la, vector<double> &data){ //2d kagome
    
    data.clear();
    data.resize(1);

    double defect = 0;
    //int count = 0;
    int Num_triangles = la->N / 3;
    for (int i = 0; i < la->N; i++){    
        if (la->sub[i] == 0){ // target each tetrahedron by 0-sublattice sites
            
            int triangle = ss->Ising[i];
            for(int j = 0; j < 2; j++){ // sum up basis tetrahedrons
                triangle += ss->Ising[la->NN[i][j]];
            }
            //cout << " tetra = " << tetrahedron << endl;
            defect += triangle;
            //defect += abs(tetrahedron); // (2,2) = 0, (1,3) = 1, (0,4) = 2
            //count++;
        }
    }
    //cout << " number of tetra = " << count << endl;
    defect /= double(Num_triangles);
    data[0] = abs(defect);
    return 0;
}

void Kagome_ice::defect_configuration(Sample* ss, Lattice* la, string name, int num){
  ofstream confi;
  ostringstream filename;
  filename << name << num << "_" << ss->ctr_output << ".txt";
  confi.open( filename.str().c_str() );
  //int Num_triangles = la->Num_sites / 3;
  for(int i = 0; i < la->N; i++){
    if (la->sub[i] == 0){
      int up_triangle = ss->Ising[i];
      up_triangle += ss->Ising[la->NN[i][0]] + ss->Ising[la->NN[i][1]];
      int down_triangle = ss->Ising[i];
      down_triangle += ss->Ising[la->NN[i][2]] + ss->Ising[la->NN[i][3]];
      down_triangle *= (-1);
      confi << la->coor[i].x << '\t' << la->coor[i].y << '\t' << up_triangle << '\t' << down_triangle << endl;
    }
  }
  confi.close();
}

double Kagome_ice::magnetization(Sample* ss, Lattice* la, vector<double> &data){

  data.clear();
  data.resize(4);

  data[0] = ss->magnetization_update[0] / double(la->N);
  data[1] = ss->magnetization_update[1] / double(la->N);
  data[2] = ss->magnetization_update[2] / double(la->N);

  double M = sqrt(data[0] * data[0] + data[1] * data[1] + data[2] * data[2]);
  data[3] = M;
  return 0;
} 

double Kagome_ice::FP_order(Sample* ss, Lattice* la, vector<double>& data){
    
    data.clear();
    data.resize(4);

    for (int i = 0; i < ss->N; i++){

      switch(la->sub[i]){
	case 0:{
		 data[0] += ss->Ising[i];
		 data[1] += (-1.0) * ss->Ising[i];
		 data[2] += (-1.0) * ss->Ising[i];
		 break;
	       }
	case 1:{
		 data[0] += (-1.0) * ss->Ising[i];
		 data[1] += ss->Ising[i];
		 data[2] += (-1.0) * ss->Ising[i];
		 break;
	       }
	case 2:{
		 data[0] += (-1.0) * ss->Ising[i];
		 data[1] += (-1.0) * ss->Ising[i];
		 data[2] += ss->Ising[i];
		 break;
	       }
      }
    }

    data[3] = data[0];

    if(data[3] < data[1])
        data[3] = data[1];
    if(data[3] < data[2])
        data[3] = data[2];


  for(int i = 0; i < 4; i++)
    data[i] /= double(ss->N);

  return 0;  
}


double Kagome_ice::domain_order(Sample* ss, Lattice* la, vector<double> &data){
    
    data.clear();
    data.resize(3,0);    // 0: wall order, 1: alpha chain order, 2: beta chain staggered order

    double *single_chain_order;
    single_chain_order = new double[ss->L];
    for (int k = 0; k < la->L; k++)
        single_chain_order[k] = 0;
    
    
    for (int i = 0; i < ss->N; i++){
        
        if (la->sub[i] == 2)
            data[0] += ss->Ising[i];
        
        else{
            //int chain_index = site[i].y / 2;
            int key = ( la->site3d[i].y + la->sub[i] ) % 2;
            switch(key){
                case 0:
                    single_chain_order[la->site3d[i].y] += ss->Ising[i];
                    break;
                case 1:
                    single_chain_order[la->site3d[i].y] += (-1) * ss->Ising[i];
                    break;
            }
        }

    }
    
    for (int j = 0; j < la->L; j++){
        
        data[1] += abs(single_chain_order[j]);
        
        data[2] += single_chain_order[j];
    }
    data[0] = abs(data[0]); 
    data[0] = data[0] / double(ss->N) * 3.0;
    data[1] = data[1] / double(ss->N) / 2.0 * 3.0;
    data[2] = abs(data[2]);
    data[2] = data[2] / double(ss->N) / 2.0 * 3.0;
    delete [] single_chain_order;    
    return 0;    
}


double Kagome_ice::H_order(Sample* ss, Lattice* la, vector<double>& data){
  
  data.clear();
  data.resize(1);
  
  double result = 0;
  double phase_cos = 0;
  double phase_sin = 0;
  for (int i = 0; i < la->N; i++){
    phase_cos += ss->Ising[i] * la->Bloch_phase_cos[i];
    phase_sin += ss->Ising[i] * la->Bloch_phase_sin[i];
  }
  result = phase_cos * phase_cos + phase_sin * phase_sin;
  result = sqrt(result);
  result /= double(la->N);
  data[0] = result;
  return 0;
}

double Kagome_ice::neutron_map(Sample* ss, Lattice* la, vector<double>& data){
  
  int Num_Qpoints = la->neutron_Q.size();
  data.clear();
  data.resize(Num_Qpoints);

  vector< vector<double> > M_Q(Num_Qpoints, vector<double> (6));
  for (int i = 0; i < Num_Qpoints; i++){
   for (int j = 0; j < 6; j++){
      M_Q[i][j] = 0;
   }
  }
  
  for (int  q = 0; q < Num_Qpoints; q++){
    double Qx = la->neutron_Q[q][0];
    double Qy = la->neutron_Q[q][1];

    for (int i = 0; i < la->N; i++){
      vector<double> spin_easyaxis(3);
      get_easyaxis(la->sub[i], spin_easyaxis);
      double Ex = spin_easyaxis[0];
      double Ey = spin_easyaxis[1];
      double Ez = spin_easyaxis[2];
      //double Ex = 0;
      //double Ey = 0;
      //double Ez = 0;
      double proj = Ex * Qx + Ey * Qy;
      proj /= (Qx * Qx + Qy * Qy);
      double Sx = Ex - proj * Qx;
      double Sy = Ey - proj * Qy;
      double Sz = Ez;

      M_Q[q][0] += Sx * ss->Ising[i] * la->neutron_cos[q][i];
      M_Q[q][1] += Sx * ss->Ising[i] * la->neutron_sin[q][i];
      M_Q[q][2] += Sy * ss->Ising[i] * la->neutron_cos[q][i];
      M_Q[q][3] += Sy * ss->Ising[i] * la->neutron_sin[q][i];
      M_Q[q][4] += Sz * ss->Ising[i] * la->neutron_cos[q][i];
      M_Q[q][5] += Sz * ss->Ising[i] * la->neutron_sin[q][i];
    }
    double MM = pow(M_Q[q][0],2) + pow(M_Q[q][1],2) + pow(M_Q[q][2],2) + pow(M_Q[q][3],2) + pow(M_Q[q][4],2) + pow(M_Q[q][5],2);
    data[q] = MM / double(la->N);
  }

  return 0;
}

// -----------------------------------------------------------------------------------------------------------
// 2D square ice model with nearest-neighbor interaction
// update algorithm: SSF
// measured quantity: energy, magnetization, acceptance
// 		      energy_update, magnetization_update
// initial state: random configuration, ferromagnetic orderded state 
// -----------------------------------------------------------------------------------------------------------
//Square_ice::Square_ice() {}

Square_ice::Square_ice(INFO MC, int nei){
      L = MC.lattice_size;
      N = MC.Num_sites;
      R = MC.Num_replicas;
      Num_neighbors = 2; // for ice model we need 2nd neighbors of square lattice
      J1 = -1.0; //default value, can be changed by other functions
      cout << " ------- Hamiltonian: Ising model with layer of neighbors = " << nei << endl;
}

void Square_ice::init(INFO info) {
    L = info.lattice_size;
    N = info.Num_sites;
    R = info.Num_replicas;
    Num_neighbors = info.Num_neighbors;
    J1 = -1.0;
}
    
Square_ice::~Square_ice(){
}

void Square_ice::set_J1(double J){
  J1 = J;
  cout << " ------- Set parameter J1 = " << J1 << endl;
}

double Square_ice::SE(int p, Sample* ss, Lattice* la){
      //double se = ss->Ising[p];
      double se = 0;
      vector<int>::iterator NN_itr;
      for (NN_itr = la->NN[p].begin(); NN_itr != la->NN[p].end(); NN_itr++)
	se += ss->Ising[(*NN_itr)];
      
      switch(la->sub[p]){
	case 1:
	  se += ss->Ising[la->NNN[p][0]] + ss->Ising[la->NNN[p][2]];
	  break;
	case -1:
	  se += ss->Ising[la->NNN[p][1]] + ss->Ising[la->NNN[p][3]];
	  break;
      }
      se *= ss->Ising[p];
      return se * J1;
}
 
double Square_ice::ME(int p, Sample* ss){
      double me = (-1.0) * ss->Ising[p] * ss->get_field(0);
      return me;
}

void Square_ice::MCstep(Sample* ss, Lattice* la){
      
      double acceptance = 0;
      

	    for (int j = 0 ; j < N; j++){
	      if ( SSF(ss, la) )
	      acceptance++;
	    }
	    ss->acceptance_SSF = acceptance;
}

bool Square_ice::SSF(Sample* ss, Lattice* la){
      int random = (uni01_sampler() * N);
      double delta_m = (-2.0) * ss->Ising[random];
      double delta_se = (-2.0) * SE(random, ss, la);
      //double delta_me = (-2.0) * ME(random, ss);
      //double delta = delta_se + delta_me;
      double delta = delta_se;

      double weight = exp(-delta / ss->get_temperature() );
      double dice = uni01_sampler();
      
      if (dice < weight){
	ss->Ising[random] *= -1;
	ss->energy_update += delta;
	ss->magnetization_update[0] += delta_m;
	return true;
      }
      else{
	return false;
      }      
}

double Square_ice::total_energy(Sample* ss, Lattice* la){
  double eng = 0;
  
  for (int i = 0; i < N; i++)
      eng += SE(i,ss,la);
  return eng / double(2.0);
} 

double Square_ice::magnetization(Sample* ss, Lattice* la){
  double mag = 0;
  for (int i = 0; i < N; i++)
    mag += ss->Ising[i];
  return mag;
}

double Square_ice::obs_energy(Sample* ss, Lattice* la, vector<double> &data){
  data.clear();
  data.resize(1);
  double e = ss->energy_update;
  e /= double(ss->N);
  //e *= 2.0;
  data[0] = e;
  return e;
}

double Square_ice::obs_magnetization(Sample* ss, Lattice* la, vector<double> &data){
  data.clear();
  data.resize(2);
  double m = ss->magnetization_update[0];
  m /= double(ss->N);
  data[0] = m;
  data[1] = abs(m);
  return m;
}

double Square_ice::obs_acceptanceSSF(Sample* ss, Lattice* la, vector<double> &data){
  data.clear();
  data.resize(1);
  double a = ss->acceptance_SSF;
  a /= double(ss->N);
  data[0] = a;
  return a;
}

double Square_ice::defect_density(Sample* ss, Lattice* la, vector<double> &data){
  data.clear();
  data.resize(1);
  double defect = 0;
  for (int i = 0; i < la->N; i++){
    if (la->sub[i] == 1){
      defect += ss->Ising[i] + ss->Ising[la->NN[i][0]] + ss->Ising[la->NN[i][1]] + ss->Ising[la->NNN[i][0]];
    }
  }
  defect /= double(la->N);
  defect *= 2.0;

  data[0] = abs(defect);
  return defect;
}


void Square_ice::initialization(Sample* ss, Lattice* la, int initial_state){
  
  switch(initial_state){
    case 1:  // all spins up
	     for(int i = 0; i < N; i++)
	       ss->Ising[i] = 1;
	     break;
    case 0:  // random initial state
	     for (int i = 0; i < N; i++)
	       ss->Ising[i] = 2 * int( uni01_sampler() * 2) - 1;
	   
  }
}

void Square_ice::initialize_observable(Sample* ss, Lattice* la, double t, vector<double> h){
  ss->set_temperature(t);
  ss->resize_field(1);
  ss->set_field(h[0], 0);

  ss->resize_magnetization(1);
  ss->resize_order(0);

  ss->set_energy(total_energy(ss, la));
  ss->set_magnetization(0, magnetization(ss, la));

}

// -------------------------------------------------------------------------------------------------------------
//
// Other useful global functions 
//
// -------------------------------------------------------------------------------------------------------------


