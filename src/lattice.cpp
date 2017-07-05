#include "lattice.hpp"

Square::Square() {
}

Square::Square(int lattice_size, int& Num_sites, int LN){
      name = "square2d";
      Num_neighbors = LN;
      L = lattice_size;
      N = L * L;
      Num_sites = N;
      construct_lattice();
      find_neighbors();
      if (LN > 1)
        find_neighbors2();
      output_lattice();
      cout << " ------- Lattice: two dimensional square lattice " << endl;
}

void Square::init_from_size(int lattice_size, int num_neighbors) {
      name = "square2d";
      Num_neighbors = num_neighbors;
      L = lattice_size;
      N = L * L;
      construct_lattice();
      find_neighbors();
      if (num_neighbors > 1)
        find_neighbors2();
      output_lattice();
      cout << " ------- Lattice: two dimensional square lattice " << endl;
}

void Square::init(INFO info) {
    int L = info.lattice_size;
    int nei = info.Num_neighbors;
    init_from_size(L, nei);
}

Square::~Square(){
    delete site2d;
    site_index2d.clear();
    NN.clear();
    NNN.clear();
    sub.clear();
}

int Square::periodic_bd(int p, int d){
  int result = ( (p + d) % L + L) % L;
  return result;
}

void Square::construct_lattice(){
	
	int k = 0;
	site2d = new Point2d[N];
	
	for (int i = 0; i < N; i++)
	  sub.push_back(0);

	for (int s= 0; s < L; s++)
	  site_index2d.push_back(vector<int> (L));
	

	for (int i = 0; i < L; i++){
	    for (int j = 0; j < L; j++){
		site2d[k].y = i;
		site2d[k].x = j;
		if ( (i+j) % 2 == 0)
		  sub[k] = 1;
		else if ( (i+j) % 2 == 1)
		  sub[k] = -1;
		k++;
		}
	}
	
	for (int a = 0; a < N; a++)
	  site_index2d[site2d[a].x][site2d[a].y] = a;

}

void Square::find_neighbors(){
  
  for (int i = 0; i < N; i++){
    NN.push_back(vector<int> (4));
  }

  for (int p = 0; p < N; p++){
	int x = site2d[p].x;
	int y = site2d[p].y;
	int xp = periodic_bd(site2d[p].x, 1);
	int xm = periodic_bd(site2d[p].x, -1);
	int yp = periodic_bd(site2d[p].y, 1);
	int ym = periodic_bd(site2d[p].y, -1);
	
	NN[p][0] = site_index2d[xp][y];
	NN[p][1] = site_index2d[x][yp];
	NN[p][2] = site_index2d[xm][y];
	NN[p][3] = site_index2d[x][ym];
  }
}

void Square::find_neighbors2(){
  
  for (int i = 0; i < N; i++){
    NNN.push_back(vector<int> (4));
  }

  for (int p = 0; p < N; p++){
	int xp = periodic_bd(site2d[p].x, 1);
	int xm = periodic_bd(site2d[p].x, -1);
	int yp = periodic_bd(site2d[p].y, 1);
	int ym = periodic_bd(site2d[p].y, -1);
	
	NNN[p][0] = site_index2d[xp][yp];
	NNN[p][1] = site_index2d[xm][yp];
	NNN[p][2] = site_index2d[xm][ym];
	NNN[p][3] = site_index2d[xp][ym];
  }
}

void Square::find_neighbors_ANNN(){
  for (int i = 0; i < N; i++){
    ANNN.push_back(vector<int> (2));
  }
  for (int p = 0; p < N; p++){
    int x = site2d[p].x;
    int ypp = periodic_bd(site2d[p].y, 2);
    int ymm = periodic_bd(site2d[p].y, -2);
  
    ANNN[p][0] = site_index2d[x][ypp];
    ANNN[p][1] = site_index2d[x][ymm];
  }
}

void Square::show_lattice(){
	cout << "Show all the lattice sites of size L = " << L << " system." << endl;
	for (int i = 0; i < N; i++){
		cout << " site[" << i << "] = " << "(" << site2d[i].x << "," << site2d[i].y << "), "
		     << " sublattice = " << sub[i] << endl;
	}
	for (int a = 0; a < L; a++)
	  for (int b = 0; b < L; b++)
	    cout << " index[" << a << "][" << b << "] = " << site_index2d[a][b] << endl;

	for (int j = 0; j < N; j++)
	  for (int k = 0; k < 4; k++)
	    cout << " NN[" << j << "][" << k << "] = " << NN[j][k] << endl;
}

void Square::output_lattice(){
    ofstream la;
  la.open("SquareLattice");
    for (int i = 0; i < N; i++){
      la << site2d[i].x << " " << site2d[i].y << endl;
    }
    la.close();
}

//------
//
//

Kagome::Kagome(int lattice_size, int& Num_sites, int LN){
      name = "kagome";
      Num_neighbors = LN;
      L = lattice_size;
      N = 3 * L * L;
      Num_sites = N;
      construct_lattice();
      find_neighbors();
      if (LN > 1)
        find_neighbors2();
      if (LN > 2)
        find_neighbors3();
      
      output_lattice();
      cout << " ------- Lattice: two dimensional Kagome lattice is constructed. " << endl;

}
    
Kagome::~Kagome(){
    //site_index.clear();
    NN.clear();
    NNN.clear();
    NNNN.clear();
    NNNN2.clear();
    sub.clear();
}




void Kagome::construct_lattice(){

  //site.resize(N);
  //coor.resize(N);
  site3d = new Point3d[N];
  coor = new Coordinates2d[N];
  sub.resize(N);


  site_index3d.resize(L); //site_index3d[L][L][3]
  vtr3d::iterator it3d;
  vtr2d::iterator it2d;
  vtr1d::iterator it1d;
  for (it3d = site_index3d.begin(); it3d != site_index3d.end(); it3d++){
    (*it3d).resize(L);
    for (it2d = (*it3d).begin(); it2d != (*it3d).end(); it2d++){
      (*it2d).resize(3);
    }
  }

  int counter = 0;
  for(int j = 0; j < L; j++){
    for(int i = 0; i < L; i++){
      for(int k = 0; k < 3; k++){
        site_index3d[i][j][k] = counter;
        sub[counter] = k;
        site3d[counter].x = i;
        site3d[counter].y = j;
        site3d[counter].z = k;
        counter++;
      }
    }
  }

  Coordinates2d Bravais_a, Bravais_b;
  Bravais_a.x = 1.0; Bravais_a.y = 0.0;
  Bravais_b.x = 0.5; Bravais_b.y = sqrt(3) / 2.0;


  for(int i = 0; i < N; i++){
    switch(sub[i]){
      case 0:{
	       coor[i].x = site3d[i].x * Bravais_a.x + site3d[i].y * Bravais_b.x;
	       coor[i].y = site3d[i].x * Bravais_a.y + site3d[i].y * Bravais_b.y;
	       break;
	     }
      case 1:{
	       coor[i].x = (site3d[i].x + 0.5) * Bravais_a.x + site3d[i].y * Bravais_b.x;
	       coor[i].y = (site3d[i].x + 0.5) * Bravais_a.y + site3d[i].y * Bravais_b.y;
	       break;
	     }
      case 2:{
	       coor[i].x = site3d[i].x * Bravais_a.x + (site3d[i].y + 0.5) * Bravais_b.x;
	       coor[i].y = site3d[i].x * Bravais_a.y + (site3d[i].y + 0.5) * Bravais_b.y;
	       break;
	     }
    }
  }
	       
        
}

int Kagome::periodic_bd(int p, int d){
  int result = ( (p + d) % L + L) % L;
  return result;
}

void Kagome::find_neighbors(){

  for (int i = 0; i < N; i++){
    NN.push_back(vector<int> (4));
  }

 for (int p = 0; p < N; p++){

    int x = site3d[p].x; int y = site3d[p].y; 
    int xp = periodic_bd(site3d[p].x, 1); int yp = periodic_bd(site3d[p].y, 1);
    int xm = periodic_bd(site3d[p].x, -1); int ym = periodic_bd(site3d[p].y, -1);
           
        // neighbor 0,1,2,3 = right, up, left, down
        if (sub[p] == 0){
                NN[p][0] = site_index3d[x][y][1];
                NN[p][1] = site_index3d[x][y][2];
                NN[p][2] = site_index3d[xm][y][1];
                NN[p][3] = site_index3d[x][ym][2];
        
        }        
        if (sub[p] == 1){
                NN[p][0] = site_index3d[x][y][2];
                NN[p][1] = site_index3d[x][y][0];
                NN[p][2] = site_index3d[xp][ym][2];
		NN[p][3] = site_index3d[xp][y][0];
        }
        if (sub[p] == 2){
                NN[p][0] = site_index3d[x][y][0];
                NN[p][1] = site_index3d[x][y][1];
                NN[p][2] = site_index3d[x][yp][0];
                NN[p][3] = site_index3d[xm][yp][1]; 
        }
 }

} 

void Kagome::find_neighbors2(){
 

  for (int i = 0; i < N; i++){
    NNN.push_back(vector<int> (4));
  }

 for (int p = 0; p < N; p++){

    int x = site3d[p].x; int y = site3d[p].y; 
    int xp = periodic_bd(site3d[p].x, 1); int yp = periodic_bd(site3d[p].y, 1);
    int xm = periodic_bd(site3d[p].x, -1); int ym = periodic_bd(site3d[p].y, -1);
           
        // neighbor 0,1,2,3 = right, up, left, down
        if (sub[p] == 0){
                NNN[p][0] = site_index3d[x][ym][1];
                NNN[p][1] = site_index3d[xp][ym][2];
                NNN[p][2] = site_index3d[xm][yp][1];
                NNN[p][3] = site_index3d[xm][y][2];
        }
        if (sub[p] == 1){
                NNN[p][0] = site_index3d[xp][ym][0];
                NNN[p][1] = site_index3d[xp][y][2];
                NNN[p][2] = site_index3d[x][yp][0];
                NNN[p][3] = site_index3d[x][ym][2];
        }
        if (sub[p] == 2){
                NNN[p][0] = site_index3d[xp][y][0];
                NNN[p][1] = site_index3d[x][yp][1];
                NNN[p][2] = site_index3d[xm][yp][0];
                NNN[p][3] = site_index3d[xm][y][1];
        }
 }
} 

void Kagome::find_neighbors3(){

  for (int i = 0; i < N; i++){
    NNNN.push_back(vector<int> (4));
    NNNN2.push_back(vector<int> (2));
  }

 for (int p = 0; p < N; p++){

    int x = site3d[p].x; int y = site3d[p].y;
    int xp = periodic_bd(site3d[p].x, 1); int yp = periodic_bd(site3d[p].y, 1);
    int xm = periodic_bd(site3d[p].x, -1); int ym = periodic_bd(site3d[p].y, -1);
    
        // neighbor 0,1,2,3 = right, up, left, down
        if (sub[p] == 0){
                NNNN[p][0] = site_index3d[xp][y][0];
                NNNN[p][1] = site_index3d[x][yp][0];
                NNNN[p][2] = site_index3d[xm][y][0];
                NNNN[p][3] = site_index3d[x][ym][0];
                NNNN2[p][0] = site_index3d[xp][ym][0];
                NNNN2[p][1] = site_index3d[xm][yp][0];
        }
        
        if (sub[p] == 1){
                NNNN[p][0] = site_index3d[xp][ym][1];
                NNNN[p][1] = site_index3d[xp][y][1];
                NNNN[p][2] = site_index3d[xm][yp][1];
                NNNN[p][3] = site_index3d[xm][y][1];
                NNNN2[p][0] = site_index3d[x][yp][1];
                NNNN2[p][1] = site_index3d[x][ym][1];
        }

        if (sub[p] == 2){
                NNNN[p][0] = site_index3d[xp][ym][2];
                NNNN[p][1] = site_index3d[x][yp][2];
                NNNN[p][2] = site_index3d[xm][yp][2];
                NNNN[p][3] = site_index3d[x][ym][2];
                NNNN2[p][0] = site_index3d[xp][y][2];
                NNNN2[p][1] = site_index3d[xm][y][2];
        }
 }
} 

void Kagome::show_lattice(){
	cout << "Show all the lattice sites of size L = " << L << " system." << endl;
	for (int i = 0; i < N; i++){
		cout << " site[" << i << "] = " << "(" << site3d[i].x << "," << site3d[i].y << "), "
		     << " sublattice = " << sub[i] << endl;
	}
	for (int b = 0; b < L; b++)
	  for (int a = 0; a < L; a++)
            for (int c = 0; c < 3; c++)
	    cout << " index[" << a << "][" << b << "][" << c << "] = " << site_index3d[a][b][c] << endl;

	for (int j = 0; j < N; j++)
	  for (int k = 0; k < 4; k++)
	    cout << " NN[" << j << "][" << k << "] = " << NN[j][k] << endl;
}

void Kagome::output_lattice(){
    ofstream la;
  la.open("KagomeLattice");
    for (int i = 0; i < N; i++){
      la << coor[i].x << " " << coor[i].y << endl;
    }
    la.close();
}

void Kagome::Hstate_Bloch_phase(){
  
  double pi = 3.14159265358979323846;
  Bloch_phase_sin.resize(N);
  Bloch_phase_cos.resize(N);
  for (int i = 0; i < N; i++){
    double summ = 0;
    switch(sub[i]){
      case 0:
	summ += 2.0 * site3d[i].x + site3d[i].y + 0.5;
	break;
      case 1:
	summ += 2.0 * site3d[i].x + site3d[i].y;
	break;
      case 2:
	summ += 2.0 * site3d[i].x + site3d[i].y + 1.0;
	break;
    }
  //Bloch_phase_sin[i] = sin(2.0 * pi / 3.0 * summ);
  //Bloch_phase_cos[i] = cos(2.0 * pi / 3.0 * summ);
  Bloch_phase_sin[i] = sin(pi * summ);
  Bloch_phase_cos[i] = cos(pi * summ);
  }
}

void Kagome::neutron_initialization(vector<vector<double> > wavevectors){
  
  neutron_Q = wavevectors;
  int Num_Qpoints = wavevectors.size();
  double pi = 3.14159265358979323846;
  double sq = sqrt(3);

  for (int i = 0; i < Num_Qpoints; i++){
    neutron_cos.push_back(vector<double> (N));
    neutron_sin.push_back(vector<double> (N));
  }
  
  
  for(int i = 0; i < Num_Qpoints; i++){
    double Qx = -1.0 * wavevectors[i][0] - wavevectors[i][1];
    double Qy = sq * wavevectors[i][0] - wavevectors[i][1] / sq;
    neutron_Q[i][0] = Qx;
    neutron_Q[i][1] = Qy;
    
    for (int j = 0; j < N; j++){
      double QR = 2.0 * pi * Qx * coor[j].x + 2.0 * pi * Qy *coor[j].y;
      neutron_cos[i][j] = cos(QR);
      neutron_sin[i][j] = sin(QR);
    }
  }
}
