project (MyProject)
cmake_minimum_required(VERSION 2.8)

find_package(PythonLibs)
include_directories (${PYTHON_INCLUDE_DIRS})

find_package(Boost 1.45.0 COMPONENTS python)
include_directories (${Boost_INCLUDE_DIRS})


add_library (
    sawmaze SHARED  
    saw.cc
)

target_link_libraries (sawmaze
    boost_python
    ${PYTHON_LIBRARIES}
    ${Boost_LIBRARIES}
)