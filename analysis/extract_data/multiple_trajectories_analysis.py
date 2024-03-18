import numpy as np
import matplotlib.pyplot as plt
from Multiple_data_extractor import Multiple_data_extractor


list_simu=['out_BwULM.h5','out_MRJfE.h5','out_svUGq.h5','out_V0V02.h5']

## data extractor for a list of simulation
MySimu=Multiple_data_extractor(list_simu)
## MySimu.data is a list of Data_extractor object
MySimu0=MySimu.data[0];
MySimu0.print_shape();

## plot mass center position for all the simulation
MySimu.plot_mass_center_pos()

## plot the mean on all simulation of var_distance_from_origin and var_distance_from_origin for all simulation
## plot_var_distance_from_origin(nb_group,group_type => as Data_extractor, here all simulations are separate in 2 group by diameter
## col=[[color for the first group[mean,simu]] [color for the second group[mean,simu]]]
## here mean of var distance for the first group is in black and var distance for one imulation is in grey .. etc..
##lab=label,time_log=give time axis in log scale) 
MySimu.plot_var_distance_from_origin(2,'diameter',col=[['k','0.5'],['b','c']])

## plot the mean on all simulation of var_distance_from_origin and its polyfit like in Data_extractor
MySimu.plot_var_distance_from_origin_polyfit(tau=[2e4,4e5])

## same as plot_var_distance_from_origin but for variance of the dispersion of the pair
MySimu.plot_var_dispersion_pair(2,'diameter',col=[['k','0.5'],['b','c']])


plt.show(block = True)

