
import numpy as np
import matplotlib.pyplot as plt
from Data_extractor import Data_extractor



#print('position')
#print(mySimulation.floe_pos[:,0,0])
#print('trajectories')
#print(mySimulation.floe_trajectories()[:,0,0])
#print('speed position')
#print(mySimulation.floe_speed_pos()[:,0,0])
#print('speed')
#print(mySimulation.floe_speed[:,0,0])
#print('speed traj')
#print(mySimulation.floe_speed_traj()[:,0,0])


## with 1 floe with and without OBL ##
Simu=Data_extractor('./out_JmDeG.h5')
SimuOBL=Data_extractor('./out_uxM19.h5')


fig = plt.figure(figsize=(8,6))
plt.title('mass center position')
ax1 = fig.add_subplot(211,)
## plot mass center position(col=color of the line,lab=label) 
Simu.plot_mass_center_pos('b','without OBL')
plt.legend(loc='lower left')
ax2 = fig.add_subplot(212,sharey=ax1)
SimuOBL.plot_mass_center_pos('r','with OBL')
plt.legend(loc='lower left')

fig = plt.figure(figsize=(8,6))
plt.title('FFT on speed')
ax1 = fig.add_subplot(211,)
## plot_fft_speed(floe_num)indice of the floe,col=color of the line,lab=label,omega_log=give the frequency axix is log scale) :
Simu.plot_fft_speed(0,'b','without OBL',True)
plt.legend(loc='upper right')
ax2 = fig.add_subplot(212,sharey=ax1)
SimuOBL.plot_fft_speed(0,'r','with OBL',True)
plt.legend(loc='upper right')



### with many floe #####
Simu0=Data_extractor('./out_hehVr.h5')

## print different shape
print("data shapes")
Simu0.print_shape();

## print floe eare and diameter
print("floes area")
print(Simu0.floe_area)
print("floes diameter")
print(Simu0.floe_diameter)

## list of floe sorted by size 
# by area
print("indices of floes by area")
print(Simu0.ind_sorted_floe_area)
# by diameter
print("indices of floes by diameter")
print(Simu0.ind_sorted_floe_diam)

fig = plt.figure(figsize=(8,5))
plt.title('mass center position')
Simu0.plot_mass_center_pos('b',' ')


fig = plt.figure(figsize=(8,5))
plt.title('Trajectories of floes')
## plot_traj(floe_num=indice of the floe,col=color of the line,lab=label) :
Simu0.plot_traj(0,'k','floe num 0')
Simu0.plot_traj(10,'r','floe num 10')
Simu0.plot_traj(20,'b','floe num 20')
Simu0.plot_traj(-1,'g','floe num -1')
plt.legend(loc='lower right')

fig = plt.figure(figsize=(8,5))
plt.title('variance distance from origine')
## plot_var_distance_from_origin_polyfit(col=color of the line ,tau=separate time axis in 3 phase [0,tau[0]] [tau[0],tau[1]] [tau[1],end_time],
## signal_log=give log10(var_distance_from_origine),time_log=give time axis in log scale) 
Simu0.plot_var_distance_from_origin_polyfit(tau=[13000,200000]);
plt.legend(loc='lower right')


##plot_var_distance_from_origin(nb_group, group_type => the floe group is separate in nb_roup groups sorted by 'area' or 'diameter' or 'no' (no sorted)
## col=color of the line,lab=label,omega_log=give the frequency axix is log scale) :
fig = plt.figure(figsize=(8,10))
ax1 = fig.add_subplot(311,)
plt.title('variance distance from originefor group sorted by diameter')
Simu0.plot_var_distance_from_origin(4,'diameter',['k','r','b','g'])
plt.legend(loc='upper left')
ax2 = fig.add_subplot(312,sharey=ax1)
plt.title('dispersion pair for group sorted by area')
Simu0.plot_var_distance_from_origin(4,'area',['k','r','b','g'])
plt.legend(loc='upper left')
ax3 = fig.add_subplot(313,sharey=ax1)
plt.title('dispersion pair for group sorted by area')
Simu0.plot_var_distance_from_origin(1,'no',['k','r','b','g'])
plt.legend(loc='upper left')

## plot_var_dispersion_pair => same as plot_var_distance_from_origin but for variation of the dispersion of the 'pair'
fig = plt.figure(figsize=(8,5))
plt.title('dispersion pair')
Simu0.plot_var_dispersion_pair(1,'no');
plt.legend(loc='lower right')

fig = plt.figure(figsize=(8,6))
ax1 = fig.add_subplot(211,)
plt.title('variance distance from originefor group sorted by diameter')
Simu0.plot_var_dispersion_pair(4,'diameter',['k','r','b','g'])
ax1.set_xscale('log') 
plt.legend(loc='upper left')
ax2 = fig.add_subplot(212,sharey=ax1)
plt.title('dispersion pair for group sorted by area')
Simu0.plot_var_dispersion_pair(4,'area',['k','r','b','g'])
ax2.set_xscale('log') 
plt.legend(loc='upper left')

plt.show(block = True)

