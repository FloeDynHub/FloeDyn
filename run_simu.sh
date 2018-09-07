#!/bin/bash

source /applis/site/nix.sh

input_file=$1
nb_time_steps=$2

#OAR --project siconos
#OAR -l /nodes=1/cpu=1/core=2,walltime=02:01:00

# Choisir l'executable (chemin absolu!)
EXE=FLOE_MPI

# Liste des noeuds
export NODES=`awk -v ORS=, '{print}' $OAR_FILE_NODES|sed 's/,$//'`

# Number of cores
nbcores=`cat $OAR_NODE_FILE|wc -l`
# Number of nodes
nbnodes=`cat $OAR_NODE_FILE|sort|uniq|wc -l`
#Name of the first node
firstnode=`head -1 $OAR_NODE_FILE`
#Number of cores allocated on the first node (it is the same on all the nodes)
pernode=`grep "$firstnode\$" $OAR_NODE_FILE|wc -l`

rundir=/nfs_scratch/$USER/floedyn/F_${OAR_JOB_ID}

echo ${OAR_JOB_ID} ${HOSTNAME} $input_file $nb_time_steps ${OAR_JOB_NAME}>> jobs_params
mkdir -p rundir
#cd rundir

#INPUTS_PATH=/home/$USER/Softs/Floe_cpp/io/inputs
which mpirun
echo $input_file
echo $nb_time_steps
which $EXE
mpirun --machinefile $OAR_NODE_FILE -np $nbcores $EXE -i $input_file -t $nb_time_steps
