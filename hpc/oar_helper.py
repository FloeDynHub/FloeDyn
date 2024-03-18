#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Prints command for submitting a job on froggy/dahu or any cluster running OAR with OpenMPI.
This helper can calculate the number of cores needed for a given number of ocean parcels.
--nbside : number of ocean parcels per side = 5 for a 5x5 grid
"""

def print_oar_commands(**opts):
    N = opts.get("nbside")
    nb_proc = 4 * N * (N-1) + 2
    nb_core_per_node = 16
    nb_nodes = (nb_proc + nb_core_per_node -1) / nb_core_per_node
    oar_cmd = "oarsub -I --project=floedyn -l /nodes={0},walltime=00:30:00".format(nb_nodes)
    run_cmd = " ".join([
        "mpirun",
        "-np {}".format(nb_proc),
        "--rank-by socket:span",
        "-x LD_LIBRARY_PATH",
        "--prefix $openmpi_DIR",
        "--machinefile $OAR_NODE_FILE",
        '-mca plm_rsh_agent "oarsh"',
        "build/FLOE_MPI",
        "<args...>"
    ])
    print(oar_cmd)
    print(run_cmd)


####################
# Running module : #
####################


def run():
    import argparse
    parser = argparse.ArgumentParser(description='froggy commands')
    parser.add_argument('-n', '--nbside', type=int, default=1, help='number of ocean parcels per side')
    OPTIONS = parser.parse_args()

    print_oar_commands(**vars(OPTIONS))

if __name__ == "__main__":
    run()

