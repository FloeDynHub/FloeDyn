#!/bin/bash
# Build the Boost FloeDyn needs (1.76, same as docker/install/boost.sh) into $HOME/install-gnu-4.7/boost_1_76,
# with the floe profile's gcc. Guix offers no Boost 1.76 (only >= 1.88, whose boost::geometry traits
# changed and break src/floe/geometry/core/access.hpp), so it is built by hand, once per account.
# Only system/program_options/thread are compiled (boost::geometry is header-only): a few minutes.
#
# Run it inside a job, not on the frontend (heavy processes get killed there), after guix_init.bash
# has built the floe profile once:
#     oarsub -S ./hpc/build_boost.sh      (from FloeDyn/; downloads the tarball first if needed)
#OAR -n floe_build_boost
#OAR --project pr-sasip
#OAR -l /core=8,walltime=1:00:00
#OAR --stdout build_boost.out
#OAR --stderr build_boost.err

set -e
source /applis/site/guix-start.sh
# gcc 14 of the floe profile (same ABI as the FloeDyn build) with its search paths (CPATH... for the kernel
# headers): a mere PATH is not enough, g++ then misses linux/types.h
GUIX_PROFILE=${GUIX_USER_PROFILE_DIR}/floe; source ${GUIX_PROFILE}/etc/profile

BOOST_ROOT="${BOOST_ROOT:-$HOME/install-gnu-4.7/boost_1_76}"
SRC_DIR="$(dirname "${BOOST_ROOT}")"
mkdir -p "${SRC_DIR}" && cd "${SRC_DIR}"
[ -f boost_1_76_0.tar.bz2 ] || wget https://sourceforge.net/projects/boost/files/boost/1.76.0/boost_1_76_0.tar.bz2/download -O boost_1_76_0.tar.bz2
echo "f0397ba6e982c4450f27bf32a2a83292aba035b827a5623a14636ea583318c41  boost_1_76_0.tar.bz2" | sha256sum -c
rm -rf boost_1_76_0 && tar --bzip2 -xf boost_1_76_0.tar.bz2 && cd boost_1_76_0
# Boost < 1.77 / glibc 2.34+ incompatibility: PTHREAD_STACK_MIN may be a function call, not a constant
sed -i 's/#if PTHREAD_STACK_MIN > 0/#if defined(PTHREAD_STACK_MIN) \&\& PTHREAD_STACK_MIN > 0/' boost/thread/pthread/thread_data.hpp
./bootstrap.sh --with-toolset=gcc --prefix="${BOOST_ROOT}" --with-libraries=system,program_options,thread
./b2 install -j8 toolset=gcc
echo "Boost 1.76 installed in ${BOOST_ROOT}"
