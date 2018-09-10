#!/bin/bash

# 
oarsub -S "./run_simu.sh $HOME/Softs/Floe_Cpp/io/inputs/in_11200f_80p_UKXcR.h5 20" --name=floedyn_2018_09 -l /nodes=2/cpu=1/core=6,walltime=02:01:00

