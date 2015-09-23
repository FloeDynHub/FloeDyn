#!/bin/sh
NODEBUG="-DNDEBUG -DBOOST_UBLAS_NDEBUG"
# WARNING="-Wall -Wextra"  -fopenmp
WARNING=""
g++ -std=c++11 ${WARNING} -I./../../src/ -I/usr/local/include -I/usr/local/include/eigen3 -O3 ${NODEBUG} -pipe -fopenmp -mtune=native STEST_PERIODIC.cpp -lmatio -lhdf5 -lhdf5_cpp -lboost_timer -lboost_system -o PERIO
# g++ -g -std=c++11 ${WARNING} -I./../../src/ -I/usr/local/include -I/usr/local/include/eigen3 -O2 -pipe -fopenmp STEST_PERIODIC.cpp -lmatio -lhdf5 -lhdf5_cpp -lboost_timer -lboost_system -o PERIOG
