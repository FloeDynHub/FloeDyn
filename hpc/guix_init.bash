source /applis/site/guix-start.sh
refresh_guix floe
export CFLAGS=-I/home/qjouet-ext/install-gnu-4.7/boost_1_72/include
export LDFLAGS=-L/home/qjouet-ext/install-gnu-4.7/boost_1_72/lib
export LD_LIBRARY_PATH=/home/qjouet-ext/install-gnu-4.7/boost_1_72/lib

python3 ./waf configure --gcc --default-search-path $GUIX_USER_PROFILE_DIR/floe
guix remove -p $GUIX_USER_PROFILE_DIR/floe boost

# Example mpi run command :
# mpirun -np 10 build/FLOE_MPI io/inputs/in_2800f_75p_tpCrm.h5 -t 1000 -o 60 --obl 0 --fmodes 0 0