# Lightweight RUNTIME environment for FloeDyn on CIMENT/Gricad clusters — for use INSIDE OAR jobs.
#
# Unlike hpc/guix_init.bash (which BUILDS the 'floe' profile: guix time-machine ... package + guix remove),
# this script only ACTIVATES the already-built profile read-only. It performs NO profile mutation, so many
# jobs can source it concurrently without hitting "profile is locked by another process".
#
# Build the profile ONCE (login node / dahu-workflow1):  source hpc/guix_init.bash
# Then in each OAR job script:                            source hpc/guix_env.bash ; build/FLOE ...
#
# Usage:  cd FloeDyn && source hpc/guix_env.bash

source /applis/site/guix-start.sh   # defines GUIX_USER_PROFILE_DIR (no lock, no mutation)
PROFILE="${GUIX_USER_PROFILE_DIR}/floe"

# FloeDyn binaries find the profile's libraries (hdf5, matio, ...) via their rpath to the /gnu/store paths,
# so no LD_LIBRARY_PATH is needed for those. Only the hand-built Boost 1.72 (outside Guix) must be found.
BOOST_ROOT="${BOOST_ROOT:-$HOME/install-gnu-4.7/boost_1_72}"
export LD_LIBRARY_PATH="${BOOST_ROOT}/lib${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"

# For MPI jobs: make the profile's mpirun/openmpi available.
export PATH="${PROFILE}/bin:${PATH}"
