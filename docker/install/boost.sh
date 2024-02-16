cd deps
wget https://boostorg.jfrog.io/artifactory/main/release/1.72.0/source/boost_1_72_0.tar.bz2
tar --bzip2 -xf boost_1_72_0.tar.bz2
cd boost_1_72_0
./bootstrap.sh
./b2 install