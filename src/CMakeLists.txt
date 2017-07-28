SET(INCLUDE_DIR ../include)
INCLUDE_DIRECTORIES(${INCLUDE_DIR})

PYTHON_ADD_MODULE(icegame sqice.cc hamiltonian.cpp lattice.cpp sample.cpp timer.cpp rng.cpp observable.cpp)
FILE(COPY test_icegame.py DESTINATION .)
ADD_TEST(NAME src COMMAND ${PYTHON_EXECUTABLE} test_icegame.py)
