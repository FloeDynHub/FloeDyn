#!/bin/sh
NODEBUG="-DNDEBUG -DBOOST_UBLAS_NDEBUG"
# WARNING="-Wall -Wextra"  -fopenmp
WARNING=""
CPPLIST="../src/floe/collision/matlab/detector.cpp ../src/floe/collision/matlab/periodic_detector.cpp ../src/floe/dynamics/dynamics_manager.cpp ../src/floe/dynamics/periodic_dynamics_manager.cpp ../src/floe/generator/generator.cpp ../src/floe/lcp/solver/LCP_solver.cpp ../src/floe/lcp/solver/lemke_eigen.cpp ../src/floe/lcp/solver/lexicolemke.cpp"
g++-5 FLOE.cpp ${CPPLIST} -D PBC -std=c++11 ${WARNING} -I./../src/ -I/usr/local/include -I/usr/local/include/eigen3 -O3 ${NODEBUG} -pipe -fopenmp -mtune=native -lboost_system -lmatio -lhdf5 -lhdf5_cpp -lCGAL -lgmp -lmpfr -lboost_thread -o FLOE_PBC -v
# g++ -g -std=c++11 ${WARNING} -I./../src/ -I/usr/local/include -I/usr/local/include/eigen3 -O2 -pipe -fopenmp FLOE_PBC.cpp -lmatio -lhdf5 -lhdf5_cpp -lboost_timer -lboost_system -o FLOE_PBC

# Build FLOE :
# g++-5 FLOE.cpp -std=c++11 ${WARNING} -I./../src/ -I/usr/local/include -I/usr/local/include/eigen3 -O3 ${NODEBUG} -pipe -fopenmp -mtune=native -lboost_system -lmatio -lhdf5 -lhdf5_cpp -lCGAL -lgmp -lmpfr -lboost_thread -o FLOE


