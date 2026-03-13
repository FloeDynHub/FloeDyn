cd deps
wget https://sourceforge.net/projects/boost/files/boost/1.72.0/boost_1_72_0.tar.bz2/download -O boost_1_72_0.tar.bz2
tar --bzip2 -xf boost_1_72_0.tar.bz2
cd boost_1_72_0
./bootstrap.sh
./b2 install
