#!/bin/sh
#GCCOPT="-floop-interchange -floop-block -floop-parallelize-all"
# NODEBUG="-DNDEBUG -DBOOST_UBLAS_NDEBUG"
NODEBUG=""
# clang++ -g -std=c++11 -Wall -Wextra -I../../../ -I/usr/local/include/eigen3 -O3 ${NODEBUG} -pipe -march=native -mtune=native TEST_detector_from_mat.cpp -lmatio -lhdf5 -lboost_timer -lboost_system -o TEST_detector_from_mat


g++ -g -std=c++11 -Wall -Wextra -I../../../ -I/usr/local/include/eigen3 -O3 ${NODEBUG} -pipe -mtune=native -fopenmp TEST_detector_from_mat.cpp -lmatio -lhdf5 -lboost_timer -lboost_system -o TEST_detector_from_mat
