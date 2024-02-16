
### test extraction hdf5

import h5py
from math import pi
import numpy as np
import matplotlib.pyplot as plt
from Data_extractor import Data_extractor



class Multiple_data_extractor(object):
	def __init__(self,filename_string_list):
		self.filename=filename_string_list
		self.nb_simu=len(self.filename)
		self.datah5=[]
		self.data=[]
		for i in range(0,self.nb_simu):
			self.datah5.append(h5py.File(self.filename[i], 'r'));
			self.data.append(Data_extractor(self.filename[i]));
		self.nb_time=self.data[0].nb_times; ### must be the same for all simulation !!!

		
		### calcul les variation de la distance à l'origine sur chaque sous groupe
	def var_dist_to_origin(self,nb_group=1,group_type='no'):
		var=np.zeros((self.nb_simu,nb_group,self.nb_time))
		##diff=np.zeros((self.nb_simu,self.nb_time-1))
		for i in range(0,self.nb_simu):
			var[i,:,:]=self.data[i].var_distance_from_origin(nb_group,group_type)
			##diff[i,:]=self.data[i].diffusivity(var[i,:]);
		return var
		
		
	def mean_var_dist(self,nb_group=1,group_type='no'):
		var=self.var_dist_to_origin(nb_group,group_type); ##((simu number,group number,time))
		mean_var=np.zeros((nb_group,self.nb_time));
		for i in range(0,self.nb_time):
			for j in range(0,nb_group):
				mean_var[j,i]=np.mean(var[:,j,i])
		return [var,mean_var]
		
		
	def var_dispersion_pair(self,nb_group=1,group_type='no'):
		var=np.zeros((self.nb_simu,nb_group,self.nb_time))
		for i in range(0,self.nb_simu):
			var[i,:,:]=self.data[i].var_dispersion_pair(nb_group,group_type)
		return var
		
	def mean_var_dispersion_pair(self,nb_group=1,group_type='no'):
		var=self.var_dispersion_pair(nb_group,group_type); ##((simu number,group number,time))
		mean_var=np.zeros((nb_group,self.nb_time));
		for i in range(0,self.nb_time):
			for j in range(0,nb_group):
				mean_var[j,i]=np.mean(var[:,j,i])
		return [var,mean_var]
		

	def plot_var_distance_from_origin(self,nb_group=1,group_type='no',col=[['k','0.5']],lab='var distance from origin',time_log=True) :
		var_mean=self.mean_var_dist(nb_group,group_type) ##(var[simu,group,time] mean_var[group,time])
		dt=self.data[0].dt;
		fig = plt.figure(figsize=(8,5))

		for j in range(0,nb_group) :
			for i in range(0,self.nb_simu):
				## plot variance for each simulation in the group
				plt.plot(self.data[i].time[1:],np.log10(var_mean[0][i,j,1:]),color=col[j%len(col)][1]);
			## plot mean variance for each group
			plt.plot(self.data[0].time[1:],np.log10(var_mean[1][j,1:]),color=col[j%len(col)][0],label=lab+'group number %d'%j)
		plt.ylabel(r'$\log_{10}(<X(t)-X(0)>)(t)$')
		plt.xlabel(r'$t$')
		if(time_log) :
			plt.xscale('log')
		plt.legend(loc='upper left')

	def plot_var_distance_from_origin_polyfit(self,tau=[7000,70000],col=['k','0.5','b','r','g'],time_log=True,signal_log=True) :
		var_mean=self.mean_var_dist(1,'no') ##(var[simu,group,time] mean_var[group,time])
		time=self.data[0].time;
		dt=self.data[0].dt;
		fig = plt.figure(figsize=(8,5))
		for i in range(0,self.nb_simu):
			plt.plot(self.data[i].time[1:],np.log10(var_mean[0][i,0,1:]),color=col[1],label='simu %d'%i);
		#var=np.zeros(time.shape[0])
		#for j in range(0,time.shape[0]):
		#	var=np.mean(var_mean[0][:,0,j]);
		var=np.mean(var_mean[0][:,0,:],axis=0);
		
		if len(tau)==2:
			phase1=time[1:int(tau[0]/dt)];
			coeff1=np.polyfit(np.log10(phase1),np.log10(var[1:phase1.shape[0]+1]),1);
			phase2=time[int(tau[0]/dt):int(tau[1]/dt)];
			coeff2=np.polyfit(np.log10(phase2),np.log10(var[(phase1.shape[0]+1):(phase2.shape[0]+phase1.shape[0]+1)]),1);
			phase3=time[int(tau[1]/dt):];
			coeff3=np.polyfit(np.log10(phase3),np.log10(var[(phase2.shape[0]+phase1.shape[0]+1):]),1);
			
		if (signal_log):
			if len(tau)==2:
				plt.plot(phase1,coeff1[1]+np.log10(phase1)*coeff1[0],color=col[2],label='$10^{%d}t^{%d}$'%(coeff1[1],coeff1[0]));
				plt.plot(phase2,coeff2[1]+np.log10(phase2)*coeff2[0],color=col[3],label='$10^{%d}t^{%d}$'%(coeff2[1],coeff2[0]));
				plt.plot(phase3,coeff3[1]+np.log10(phase3)*coeff3[0],color=col[4],label='$10^{%d}t^{%d}$'%(coeff3[1],coeff3[0]));
			plt.plot(time[1:],np.log10(var[1:]),color=col[0],label='mean var distance');
			plt.ylabel(r'$\log_{10}(<X(t)-X(0)>)(t)$')
		else:
			if len(tau)==2:
				plt.plot(phase1,np.power(10.0,coeff1[1])*np.power(phase1,coeff1[0]),color=col[2],label='$e^{%d}t^{%d}$'%(coeff1[1],coeff1[0]));
				plt.plot(phase2,np.power(10.0,coeff1[1])*np.power(phase2,coeff2[0]),color=col[3],label='$e^{%d}t^{%d}$'%(coeff2[1],coeff2[0]));
				plt.plot(phase3,np.power(10.0,coeff1[1])*np.power(phase3,coeff3[0]),color=col[4],label='$e^{%d}t^{%d}$'%(coeff3[1],coeff3[0]));
			plt.plot(time,var[:],color=col[0],label='mean var distance');
			plt.ylabel(r'$\log_{10}(<X(t)-X(0)>)(t)$')
		plt.xlabel(r'$t$')
		if (time_log) :
			plt.xscale('log')
		plt.legend(loc='lower right')



	def plot_mass_center_pos(self,color_list=['k','b','r','g','y','m'],label='mass center pos') :
		fig = plt.figure(figsize=(8,5))
		for i in range(0,self.nb_simu) :
			self.data[i].plot_mass_center_pos(col=color_list[i%len(color_list)],lab=label+" simu %d"%i)
		plt.legend(loc='upper left')



	def plot_var_dispersion_pair(self,nb_group=1,group_type='no',col=[['k','0.5']],lab='var dispersion pair',time_log=True) :
		var_mean=self.mean_var_dispersion_pair(nb_group,group_type) ##(var[simu,group,time] mean_var[group,time])
		dt=self.data[0].dt;
		time=self.data[0].time;
		fig = plt.figure(figsize=(8,5))
		for j in range(0,nb_group) :
			for i in range(0,self.nb_simu):
				## plot variance for each simulation in the group
				plt.plot(self.data[i].time[1:],np.log10(var_mean[0][i,j,1:]),color=col[j%len(col)][1]);
			## plot mean variance for each group
			plt.plot(time[1:],np.log10(var_mean[1][j,1:]),color=col[j%len(col)][0],label=lab+'group number %d'%j)
		plt.ylabel(r'$\log_{10}(<X(t)-X(0)>)(t)$')
		plt.xlabel(r'$t$')
		if (time_log) :
			plt.xscale('log')
		plt.legend(loc='upper left')

	def plot_var_dispersion_pair_asymptote(self,tau=[7000,70000],col=['k','0.5','b','r','g'],time_log=True,signal_log=True) :	
		var_mean=self.mean_var_dispersion_pair(1,'no') ##(var[simu,group,time] mean_var[group,time])
		time=self.data[0].time;
		dt=self.data[0].dt;
		fig = plt.figure(figsize=(8,5))
		for i in range(0,self.nb_simu):
			plt.plot(self.data[i].time[1:],np.log10(var_mean[0][i,0,1:]),color=col[1],label='simu %d'%i);
		var=np.mean(var_mean[0][:,0,:],axis=0);
		
		if len(tau)==2:
			phase1=time[1:int(tau[0]/dt)+1];
			coeff1=np.polyfit(np.log10(phase1),np.log10(var[1:phase1.shape[0]+1]),1);
			phase2=time[int(tau[0]/dt):int(tau[1]/dt)];
			coeff2=np.polyfit(np.log10(phase2),np.log10(var[(phase1.shape[0]+1):(phase2.shape[0]+phase1.shape[0]+1)]),1);
			phase3=time[int(tau[1]/dt):];
			coeff3=np.polyfit(np.log10(phase3),np.log10(var[(phase2.shape[0]+phase1.shape[0]):]),1);
			
		if (signal_log):
			if len(tau)==2:
				plt.plot(phase1,coeff1[1]+np.log10(phase1)*coeff1[0],color=col[2],label='$10^{%d}t^{%d}$'%(coeff1[1],coeff1[0]));
				plt.plot(phase2,coeff2[1]+np.log10(phase2)*coeff2[0],color=col[3],label='$10^{%d}t^{%d}$'%(coeff2[1],coeff2[0]));
				plt.plot(phase3,coeff3[1]+np.log10(phase3)*coeff3[0],color=col[4],label='$10^{%d}t^{%d}$'%(coeff3[1],coeff3[0]));
			plt.plot(time[1:],np.log10(var[1:]),color=col[0],label='mean var distance');
			plt.ylabel(r'$\log_{10}(<X(t)-X(0)>)(t)$')
		else:
			if len(tau)==2:
				plt.plot(phase1,np.power(10.0,coeff1[1])*np.power(phase1,coeff1[0]),color=col[2],label='$e^{%d}t^{%d}$'%(coeff1[1],coeff1[0]));
				plt.plot(phase2,np.power(10.0,coeff1[1])*np.power(phase2,coeff2[0]),color=col[3],label='$e^{%d}t^{%d}$'%(coeff2[1],coeff2[0]));
				plt.plot(phase3,np.power(10.0,coeff1[1])*np.power(phase3,coeff3[0]),color=col[4],label='$e^{%d}t^{%d}$'%(coeff3[1],coeff3[0]));
			plt.plot(time,var[:],color=col[0],label='mean var distance');			
			plt.ylabel(r'$(<X(t)-X(0)>)(t)$')
		plt.xlabel(r'$t$')
		if(time_log) :
			plt.xscale('log')
		plt.legend(loc='lower right')	


if __name__ == "__main__" :

	#myoutputfile='./out_Lh0OX.h5'
	myTraj=Trajectory_analysis(['./out_Lh0OX.h5'])
	myTraj.data[0].print_shape();
	
	
           
