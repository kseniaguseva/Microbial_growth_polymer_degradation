CXXFLAGS=-Wall -O3 -g -std=c++17 -fPIC -lpython3.11  -lboost_python3 -lboost_numpy3 -I. -I/usr/include/python3.11 -I/usr/lib/python3.11/site-packages/numpy/core/include/ -shared -Wl,-soname,compounds.so

compounds.so: compounds.hh microorganism.hh enzymes.hh compounds.cc demangle.cc numpy_bind.hh
	g++ compounds.cc demangle.cc -o compounds.so ${CXXFLAGS}
