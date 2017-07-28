#!/bin/bash
set -x
#cd ${0%%$(basename $0)}

mkdir build
cd build

if [[ "$OSTYPE" == "linux-gnu" ]]; then
    cmake -DCMAKE_BUILD_TYPE=DEBUG .. && make 
elif [[ "$OSTYPE" == "darwin"* ]]; then
    PYTHON_VERSION=`python -c "import sys;t='{v[0]}.{v[1]}'.format(v=list(sys.version_info[:2]));sys.stdout.write(t)";`
    PYTHON_LIBRARY=/usr/local/Frameworks/Python.framework/Versions/$PYTHON_VERSION/lib/libpython$PYTHON_VERSION.dylib
    PYTHON_INCLUDE_DIR=/usr/local/Frameworks/Python.framework/Versions/$PYTHON_VERSION/Headers/
    PYTHON_NUMPY_INCLUDE_DIRS=/usr/local/Frameworks/Python.framework/Versions/$PYTHON_VERSION/lib/python$PYTHON_VERSION/site-packages/numpy/core/include/
    cmake -DPYTHON_LIBRARY=$PYTHON_LIBRARY -DPYTHON_INCLUDE_DIR=$PYTHON_INCLUDE_DIR -DPYTHON_NUMPY_INCLUDE_DIRS=$PYTHON_NUMPY_INCLUDE_DIRS -DCMAKE_BUILD_TYPE=DEBUG .. && make 
else
    : #Unkown
fi
