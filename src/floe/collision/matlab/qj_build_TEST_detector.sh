#!/bin/sh
NODEBUG="-DNDEBUG -DBOOST_UBLAS_NDEBUG"
boost_root="/usr/local/boost_1_57_0/"
boost_builds="/usr/local/boost_1_57_0/bin.v2/libs/"
matio_root="/usr/local/matio-1.5.2/src/"

# clang++ -g -std=c++11 -Wall -Wextra -I../../../ -I${boost_root} -I${matio_root} -O3 ${NODEBUG} -pipe -march=native -mtune=native \
# TEST_detector.cpp -o TEST_detector \
# ${boost_builds}timer/build/darwin-4.2.1/release/link-static/threading-multi/libboost_timer.a \
# ${boost_builds}chrono/build/darwin-4.2.1/release/link-static/threading-multi/libboost_chrono.a \
# ${boost_builds}system/build/darwin-4.2.1/release/link-static/threading-multi/libboost_system.a

boost_build_path="/usr/local/boost_1_57_0/bin.v2/libs/(lib_name)/build/darwin-4.2.1/release/link-static/threading-multi/libboost_(lib_name).a"

clang++ -g -std=c++11 -Wall -Wextra -I../../../ -I${boost_root} -I${matio_root} -O3 ${NODEBUG} -pipe -march=native -mtune=native \
TEST_detector.cpp -o TEST_detector2 \
${boost_build_path//(lib_name)/timer} \
${boost_build_path//(lib_name)/chrono} \
${boost_build_path//(lib_name)/system}


# clang++ -g -std=c++11 -Wall -Wextra -I../../../ -I${boost_root} -I${matio_root} -O3 ${NODEBUG} -pipe -march=native -mtune=native \
# TEST_detector.cpp -o TEST_detector2 -lboost_timer -lboost_chrono -lboost_system