#!/bin/sh
#GCCOPT="-floop-interchange -floop-block -floop-parallelize-all"
NODEBUG="-DNDEBUG -DBOOST_UBLAS_NDEBUG"
# clang++ -g -std=c++11 -Wall -Wextra -I../../../ -O3 ${NODEBUG} -pipe -march=native -mtune=native TEST_detector_from_mat.cpp -lmatio -lhdf5 -lboost_timer -lboost_system -o TEST_detector_from_mat
boost_root="/usr/local/boost_1_57_0/"
matio_root="/usr/local/matio-1.5.2/src/"
boost_build_path="/usr/local/boost_1_57_0/bin.v2/libs/(lib_name)/build/darwin-4.2.1/release/link-static/threading-multi/libboost_(lib_name).a"
boost_build_path_dy="/usr/local/boost_1_57_0/bin.v2/libs/(lib_name)/build/darwin-4.2.1/release/threading-multi/libboost_(lib_name).dylib"
clang++ -g -std=c++11 -Wall -Wextra -I../../../ -I${boost_root} -I${matio_root} -O3 ${NODEBUG} -pipe -march=native -mtune=native \
TEST_detector_from_mat.cpp -o TEST_detector_from_mat3 \
${boost_build_path//(lib_name)/timer} \
${boost_build_path_dy//(lib_name)/chrono} \
${boost_build_path_dy//(lib_name)/system} \
/usr/local/matio-1.5.2/src/.libs/libmatio.2.dylib




