#CXX = icpc
CXX = g++
CXXFLAGS = -Wall -O3 -fPIC -std=c++11

PYTHON_VERSION = 2.7
PYTHON_INCLUDE = /usr/include/python$(PYTHON_VERSION)

OBJECT_MAIN = saw.o
LDFLAGS   = -lgsl -lgslcblas `libpng-config --ldflags`

BOOST_INC = /usr/include/boost
BOOST_LIB = /usr/local/lib

export CPLUS_INCLUDE_PATH="$CPLUS_INCLUDE_PATH:/usr/include/python2.7/"

main:
	$(CXX) $(CXXFLAGS) -shared -o libsawmaze.so saw.cc -L$(BOOST_LIB) -lboost_python -L/usr/lib/python$(PYTHON_VERSION)/ -L/usr/lib/python$(PYTHON_VERSION)/config -lpython$(PYTHON_VERSION)

game: saw.o
	$(CXX) $(CXXFLAGS) -shared -Wl,-rpath,libsawmaze.so -o libsawmaze.so saw.o $(OBJECT_ALL) -L$(BOOST_LIB) -lboost_python -L/usr/lib/python$(PYTHON_VERSION)/config -lpython$(PYTHON_VERSION)

saw.o: saw.cc
	$(CXX) $(CXXFLAGS) -I$(PYTHON_INCLUDE) -I$(BOOST_INC) -fPIC -c saw.cc

clean:
	rm -f $(OBJECT_ALL) $(OBJECT_MAIN) *.o *.so
