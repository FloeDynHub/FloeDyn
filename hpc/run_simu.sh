#!/bin/bash
#OAR -l /nodes=1/cpu=1/core=4,walltime=00:01:00
#OAR -n floedyn_example
#OAR --project floedyn

input_file=$1
nb_time_steps=$2

# Absolute path to floedyn runner 
EXE=$HOME/Floe_Cpp/build/FLOE_MPI
# Absolute path to input files 
INPUTDIR=$HOME/Floe_Cpp/io/library

nbcores=`cat $OAR_NODE_FILE|wc -l`
rundir=/bettik/$USER/floedyn/F_${OAR_JOB_ID}
echo ${OAR_JOB_ID} ${HOSTNAME} $input_file $nb_time_steps ${OAR_JOB_NAME}>> jobs_params

mkdir -p $rundir/io/inputs
mkdir -p $rundir/io/outputs
cp $INPUTDIR/DataTopaz01.mat $rundir/io/library/
cd $rundir
source /applis/site/guix-start.sh
refresh_guix floedyn
export MPIPREFIX=$GUIX_USER_PROFILE_DIR/floedyn

mpirun  --machinefile $OAR_NODE_FILE --mca orte_rsh_agent "oarsh"  -mca btl_openib_allow_ib 1 -np $nbcores --prefix $MPIPREFIX $EXE -i $input_file -t $nb_time_steps
