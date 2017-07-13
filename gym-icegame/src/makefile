#CXX = icpc
CXX = g++
CXXFLAGS = -Wall -O3 -fPIC -std=c++11

PYTHON_VERSION = 2.7
PYTHON_INCLUDE = /usr/include/python$(PYTHON_VERSION)

OBJECT_ALL = sample.o lattice.o observable.o measurement.o hamiltonian.o binning.o timer.o heattreatment.o rng.o analysis.o
OBJECT_MAIN = mc.o
OBJECT_ENG = cal_energy.o
LDFLAGS   = -lgsl -lgslcblas `libpng-config --ldflags`

BOOST_INC = /usr/include/boost
BOOST_LIB = /usr/local/lib

SQUAREICE=squareice

mcrun : $(OBJECT_ALL) $(OBJECT_MAIN)
	$(CXX) $(CXXFLAGS) $(OBJECT_ALL) $(OBJECT_MAIN)  -o mcrun

clean:
	rm -f $(OBJECT_ALL) $(OBJECT_MAIN) *.o *.so

cleandata:
	rm *.data *.txt info Summary_* SquareLattice

erun : $(OBJECT_ALL) $(OBJECT_ENG)
	$(CXX) $(CXXFLAGS) $(OBJECT_ALL) $(OBJECT_ENG) -o erun

sqice.o: sqice.cc
	$(CXX) $(CXXFLAGS) -I$(PYTHON_INCLUDE) -I$(BOOST_INC) -fPIC -c sqice.cc

icegame: sqice.o
	$(CXX) $(CXXFLAGS) -shared -Wl,-rpath,libicegame.so -o libicegame.so sqice.o $(OBJECT_ALL) -L$(BOOST_LIB) -lboost_python -L/usr/lib/python$(PYTHON_VERSION)/config -lpython$(PYTHON_VERSION)

cal_energy.o: lattice.o sample.hpp observable.hpp measurement.hpp hamiltonian.hpp binning.hpp timer.hpp heattreatment.hpp analysis.hpp cal_energy.cpp cal_energy.hpp
	$(CXX) $(CXXFLAGS) -c cal_energy.cpp -o cal_energy.o

mc.o: lattice.o sample.hpp observable.hpp measurement.hpp hamiltonian.hpp binning.hpp timer.hpp heattreatment.hpp analysis.hpp mc.cpp mc.hpp
	$(CXX) $(CXXFLAGS) -c mc.cpp -o mc.o

sample.o: sample.hpp sample.cpp rng.hpp 
	$(CXX) $(CXXFLAGS) -c sample.cpp -o sample.o

lattice.o: lattice.hpp lattice.cpp sample.hpp
	$(CXX) $(CXXFLAGS) -c lattice.cpp -o lattice.o

observable.o: observable.hpp observable.cpp sample.hpp lattice.hpp
	$(CXX) $(CXXFLAGS) -c observable.cpp -o observable.o

measurement.o: measurement.hpp measurement.cpp sample.hpp observable.hpp lattice.hpp
	$(CXX) $(CXXFLAGS) -c measurement.cpp -o measurement.o

hamiltonian.o: hamiltonian.hpp hamiltonian.cpp sample.hpp rng.hpp lattice.hpp
	$(CXX) $(CXXFLAGS) -c hamiltonian.cpp -o hamiltonian.o
 
binning.o: binning.hpp binning.cpp sample.hpp rng.hpp
	$(CXX) $(CXXFLAGS) -c binning.cpp -o binning.o

analysis.o: analysis.hpp analysis.cpp sample.hpp binning.hpp
	$(CXX) $(CXXFLAGS) -c analysis.cpp -o analysis.o

heattreatment.o: heattreatment.hpp heattreatment.cpp sample.hpp hamiltonian.hpp lattice.hpp
	$(CXX) $(CXXFLAGS) -c heattreatment.cpp -o heattreatment.o

rng.o: rng.hpp rng.cpp
	$(CXX) $(CXXFLAGS) -c rng.cpp -o rng.o 

timer.o: timer.hpp timer.cpp
	$(CXX) $(CXXFLAGS) -c timer.cpp -o timer.o
