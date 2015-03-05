#!/bin/sh
#GCCOPT="-floop-interchange -floop-block -floop-parallelize-all"
NODEBUG="-DNDEBUG -DBOOST_UBLAS_NDEBUG"
clang++ -g -std=c++11 -Wall -Wextra -I../../../ -I/home/denis/Téléchargements/eigen3 -O3 ${NODEBUG} -pipe -march=native -mtune=native TEST_detector_from_mat.cpp -lmatio -lhdf5 -lboost_timer -lboost_system -o TEST_detector_from_mat
