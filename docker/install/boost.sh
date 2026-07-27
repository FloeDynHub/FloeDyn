cd deps
wget https://sourceforge.net/projects/boost/files/boost/1.76.0/boost_1_76_0.tar.bz2/download -O boost_1_76_0.tar.bz2
tar --bzip2 -xf boost_1_76_0.tar.bz2
cd boost_1_76_0
# Boost < 1.77 / glibc 2.34+ incompatibility: PTHREAD_STACK_MIN may be a function call, not a constant
sed -i 's/#if PTHREAD_STACK_MIN > 0/#if defined(PTHREAD_STACK_MIN) \&\& PTHREAD_STACK_MIN > 0/' boost/thread/pthread/thread_data.hpp
CFLAGS="-march=x86-64" CXXFLAGS="-march=x86-64" ./bootstrap.sh
# wave and context fail on gcc:12 due to flex_string.hpp/-Wfree-nonheap-object (not needed by FloeDyn)
./b2 install cxxflags="-march=x86-64 -Wno-free-nonheap-object" -j2 || true
