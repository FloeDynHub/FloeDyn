module purge
source /applis/ciment/v2/env.bash
module use -a /home/rabatmat/modules
module load intel-devel/2015
# module load gcc/4.9.3_gcc-4.6.2 openmpi/1.6.5_gcc-4.7.2
module load python/2.7.6_gcc-4.6.2
module load boost-gcc4.7 matio-gcc4.7 eigen-gcc4.7 hdf5-1.8.15-gcc4.8.2
#module load gcc/4.8.2_gcc-4.4.6
module load cgal/4.6.3_gcc-4.6.2 mpfr/3.1.1_gcc-4.6.2 qt/4.8.3_gcc-4.6.2 cereal1.2

