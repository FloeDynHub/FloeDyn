cd deps
wget https://sourceforge.net/projects/boost/files/boost/1.72.0/boost_1_72_0.tar.bz2/download -O boost_1_72_0.tar.bz2
tar --bzip2 -xf boost_1_72_0.tar.bz2
cd boost_1_72_0
# Fix Boost 1.72 / GCC 11+ incompatibility: PTHREAD_STACK_MIN may be a function call, not a constant
sed -i 's/#if PTHREAD_STACK_MIN > 0/#if defined(PTHREAD_STACK_MIN) \&\& PTHREAD_STACK_MIN > 0/' boost/thread/pthread/thread_data.hpp
CFLAGS="-march=x86-64" CXXFLAGS="-march=x86-64" ./bootstrap.sh
./b2 install cxxflags="-march=x86-64" -j2
