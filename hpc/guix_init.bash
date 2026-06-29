# Reproducible FloeDyn build environment on CIMENT/Gricad clusters (kraken, dahu, ...).
#
# Reproducibility relies on TWO pinned files, both versioned in this repo:
#   hpc/channels.scm          -> exact Guix commit  => exact package *versions*
#   hpc/manifest_floedyn.scm  -> the package list   (versions come from the pin)
# The point is to build the profile *through* the pin with `guix time-machine`.
# (The old `refresh_guix floe` resolved the manifest against whatever Guix state
#  each user had pulled -- which is why the build broke on another account.)
#
# Usage:  cd FloeDyn && source hpc/guix_init.bash

HPC_DIR="$(cd "$(dirname "${BASH_SOURCE[0]:-hpc/guix_init.bash}")" && pwd)"
source /applis/site/guix-start.sh

PROFILE="${GUIX_USER_PROFILE_DIR}/floe"

# Build/refresh the 'floe' profile from the PINNED channels + manifest.
# Same channels.scm + manifest => identical toolchain for everyone, now or later.
guix time-machine -C "${HPC_DIR}/channels.scm" -- \
     package -p "${PROFILE}" -m "${HPC_DIR}/manifest_floedyn.scm"

# Boost is currently NOT taken from Guix: FloeDyn wants Boost 1.72 (boost::geometry),
# so we drop the manifest's boost and use a hand-built 1.72 from $HOME instead.
# NOTE: this hand-built Boost is not in the repo -> not yet reproducible for others.
# (Planned fix: take Boost from the pinned channel, or commit a build script.)
guix remove -p "${PROFILE}" boost
BOOST_ROOT="${BOOST_ROOT:-$HOME/install-gnu-4.7/boost_1_72}"
export CFLAGS=-I${BOOST_ROOT}/include
export LDFLAGS=-L${BOOST_ROOT}/lib
export LD_LIBRARY_PATH=${BOOST_ROOT}/lib
# Put the floe profile first in PATH so the build picks the profile's mpicc/mpicxx (openmpi 4.x), whose
# --showme flags drive the MPI build (wscript), rather than any system MPI.
export PATH=${PROFILE}/bin:${PATH}

python3 ./waf configure --gcc --default-search-path "${PROFILE}"

# Example mpi run command :
# mpirun -np 10 build/FLOE_MPI io/inputs/in_2800f_75p_tpCrm.h5 -t 1000 -o 60 --obl 0 --fmodes 0 0
