#CXX = icpc
CXX = g++
CXXFLAGS = -Wall -O3 -fPIC -std=c++11

OS := $(shell uname)

PYTHON_VERSION = 2.7
PYTHON_INCLUDE = /usr/include/python$(PYTHON_VERSION)

OBJECT_ALL = sample.o lattice.o hamiltonian.o timer.o rng.o 

LDFLAGS   = -lgsl -lgslcblas `libpng-config --ldflags`

BOOST_INC = /usr/include/boost
BOOST_LIB = /usr/local/lib

icegame: sqice.o $(OBJECT_ALL)
	$(CXX) $(CXXFLAGS) -shared -Wl,-rpath,icegame.so -o icegame.so sqice.o $(OBJECT_ALL) -L$(BOOST_LIB) -lboost_python -L/usr/lib/python$(PYTHON_VERSION)/config -lpython$(PYTHON_VERSION)

clean:
	rm -f $(OBJECT_ALL) $(OBJECT_MAIN) *.o* *.so *.data *.txt info Summary_* SquareLattice

sqice.o: sqice.cc
	$(CXX) $(CXXFLAGS) -I$(PYTHON_INCLUDE) -I$(BOOST_INC) -fPIC -c sqice.cc
