
### test extraction hdf5

import h5py
from math import pi
import numpy as np
import matplotlib.pyplot as plt
from numpy.fft import fft, fftshift



class Data_extractor(object):
	def __init__(self,filename_string):
		self.filename=filename_string
		self.myData= h5py.File(self.filename, 'r')
		
		self.floe_states=np.array(self.myData['/'].get('floe_states') )
		## [times,floe_id,states]
		##state=[x,y,theta,ux,uy,u_theta,impulse,x-translation,y-translation]
		
		self.nb_floes=self.floe_states.shape[1]
		self.floe_pos=self.floe_states[:,:,:2]## [times,floe_id,pos(x,y)]
		self.floe_speed=self.floe_states[:,:,3:5]## [times,floe_id,speed(x,y)]
		
		self.kinematicE=np.array(self.myData['/'].get('Kinetic Energy') )
		##[time,floe_id]
		
		self.time=np.array(self.myData['/'].get('time'))
		self.nb_times=self.time.shape[0]
		self.dt=self.time[1]-self.time[0]
		
		self.OBL_speed=np.array(self.myData['/'].get('OBL_speed'))
		##[time,(OBL_x,OBL_y)]
		self.mass_center=np.array(self.myData['/'].get('mass_center'))
		##[time,pos(x,y)]
		
		### floe_shapes [nb_of_vertices,coordinate,)]
		self.floe_shapes=[];
		for i in range(0,self.nb_floes):
			self.floe_shapes.append(np.array(self.myData['/floe_shapes/'].get(str(i))));

		self.floe_area=np.zeros(self.nb_floes); ### area of each floe
		for i in range(0,self.nb_floes):
			self.floe_area[i]=0.5*np.abs(np.dot(self.floe_shapes[i][:,0],np.roll(self.floe_shapes[i][:,1],1))-np.dot(self.floe_shapes[i][:,1],np.roll(self.floe_shapes[i][:,0],1)))
			
		self.floe_diameter=np.zeros(self.nb_floes); ### diameter of each floe
		for i in range(0,self.nb_floes):
			diam_max=0.0
			for j in range(0,self.floe_shapes[i].shape[0]-1): ###for each vertices
				for k in range(j+1,self.floe_shapes[i].shape[0]):
					diam_max=max(diam_max,np.sqrt(np.power(self.floe_shapes[i][j,0]-self.floe_shapes[i][k,0],2)+np.power(self.floe_shapes[i][j,1]-self.floe_shapes[i][k,1],2)));
			self.floe_diameter[i]=diam_max; #### va y avoir un pb la, c'est pas par copie !?		
			
		### sort floe by diameter
		self.ind_sorted_floe_diam=np.argsort(self.floe_diameter);			
		
		### sort floe by area
		self.ind_sorted_floe_area=np.argsort(self.floe_area);
		
	def print_shape(self):
		print("output file readed")
		print(self.filename)
		print("nb of time")
		print(self.nb_times)
		print("nb of floe")
		print(self.nb_floes)
		print("floe state [times,floe_id,state],state=[x,y,theta,ux,uy,u_theta,impulse,x-translation,y-translation]")		
		print(self.floe_states.shape)
		print("floe position [times,floe_id,pos(x,y)]")
		print(self.floe_pos.shape)
		print("floe speed [times,floe_id,speed(x,y)]")
		print(self.floe_pos.shape)
		print("OBL data [time,(OBL_x,OBL_y)]")
		print(self.OBL_speed.shape)
		print("floe shape = list (array[verticices,(x,y)])")		
		print(len(self.floe_shapes),self.floe_shapes[0].shape)
		
		
		
	## OPERATION ON FLOE POSITION		
		
	## calculate floe speed from floe position	
	def floe_speed_pos(self):
		##[speed in x or y,time, floe id] ] 
		floe_speed=(self.floe_states[1:,:,:]-self.floe_states[:-1,:,:])/self.dt
		#floe_speed_y=(self.floe_states[1:,:,1]-self.floe_states[:-1,:,1])/self.dt
		return floe_speed #np.array([floe_speed_x,floe_speed_y])
	
	## calculate trajectories = position - mass center position	
	def floe_trajectories(self):
		trajectories=np.zeros(self.floe_pos.shape) ## [times,floe_id,pos(x,y)]
		for i in range(0,self.floe_pos.shape[1]): ### i is the floe id
			trajectories[:,i,:]=self.floe_pos[:,i,:]-self.mass_center
		return trajectories ## [times,floe_id,pos]
		
#	## calculate speed from trajectories	
#	def floe_speed_traj(self):
#		##[speed in x or y,time, floe id] ] 
#		traj=self.floe_trajectories()
#		floe_speed=(traj[1:,:,:]-traj[:-1,:,:])/self.dt
#		return floe_speed #[times,floe_id,speed]
		

	### separate floe in sorted group of area/diamter/no_sorted, return a list of nb_group of group_size indices of floe  
	def group_list(self,nb_group=1,group_type='no'):
		group_size=np.floor(self.nb_floes/nb_group);
		list_group=[];
		if group_type=='area':
			print("floe sorted by area")
			sorted_indice=self.ind_sorted_floe_area
		elif group_type=='diameter':
			print("floe sorted by diameter")
			sorted_indice=self.ind_sorted_floe_diam
		else :
			sorted_indice=np.arange(0,self.nb_floes,dtype=int);
		for i in range(0,nb_group-1):
			list_group.append(sorted_indice[int(i*group_size):int(group_size*(i+1))]);
		
		list_group.append(sorted_indice[int((nb_group-1)*group_size):]);
		return list_group
		
		
	## calculate distance from the origin(t=0) using trajectories
	def distance_from_origin(self):
		distance=np.zeros((self.nb_times,self.nb_floes))
		traj=self.floe_trajectories()
		for i in range(0,self.nb_floes):
			distance[:,i]=np.square(traj[:,i,0]-traj[0,i,0])+np.square(traj[:,i,1]-traj[0,i,1])
		return distance
	

	## variance of the distance from origine
	def var_distance_from_origin(self,nb_group=1,group_type='no'):
		## sorted_list is a list of nb_group list of indice, sorted by the criteria 'group_type'
		sorted_group_list=self.group_list(nb_group,group_type);
		distance=self.distance_from_origin()  ###[times,floe_id]
		var=np.zeros((nb_group,self.nb_times))  ###[floe group, time]
		for j in range(0,nb_group):
			for i in range(0,self.nb_times):
				var[j,i]=np.var(distance[i,sorted_group_list[j]])
		return var
		
		## variance of the distance from origine
	def var_dispersion_pair(self,nb_group=1,group_type='no'):
		distance=self.distance_from_origin() ###[times,floe_id]
		sorted_group_list=self.group_list(nb_group,group_type);
		var=np.zeros((nb_group,self.nb_times))
		for i in range(0,self.nb_times):
			for j in range(0,nb_group) :
				group=sorted_group_list[j];
				for k in range(0,group.shape[0]-1):
					var[j,i]=var[j,i]+np.sum(np.square(distance[i,group[k]]-distance[i,group[(k+1):]]));
				var[j,:]=0.5*var[j,:]/(group.shape[0]*(group.shape[0]-1))	
		return var
		
		
	## diffusivity a changer !
##	def diffusivity(self,variance=None):
##		if variance==None :
##			variance=self.var_distance_from_origin();
##		K=0.5*(variance[1:]-variance[:-1])/self.dt;
##		return K
	
	## calculate autocorrelation of trajectories for each floe
	def autocorr_floe_trajectories(self):
		traj=self.floe_speed_traj(); ### [times,floe_id,pos(x,y)]
		autocorr_floe=np.zeros(traj.shape)  
		for i in range(0,traj.shape[1]): ## for all floe
			autocorr_floe[:,i,0]=autocorr(traj[:,i,0]);
			autocorr_floe[:,i,1]=autocorr(traj[:,i,1]);
		return autocorr_floe

	## mean autocorrelation on the floe population 
	def autocorr_trajectories(self):
		autocorr_floe=self.autocorr_floe_trajectories();##  [times,floe_id,pos(x,y)]
		autocorr_mean=np.zeros((autocorr_floe.shape[0],autocorr_floe.shape[2])) ##  [times,pos(x,y)]
		for i in range(0,autocorr_floe.shape[1]): ## for all time
			autocorr_mean[i,0]=np.mean(autocorr_floe[i,:,0]);
			autocorr_mean[i,1]=np.mean(autocorr_floe[i,:,1]);
		return autocorr_mean
				
		
	##### FFT and FFT tools

	def frequencies_shifted(self,signal_shape=1):
		freqs = np.fft.fftfreq(signal_shape,d=self.time[1]-self.time[0])*pi*2
		#print(freqs.shape)
		return fftshift(freqs)	
		
	def fftshift_each_floe(self,value) :
		## value=(time,floe,value.x or value.y)
		fft_speed=np.zeros((value.shape[1],value.shape[0]))
		for i in range(0,value.shape[1]):
		### for each floe, calculate fourier transform on the postion see as complex value FT(x(t)+iy(t))
			fft_speed[i,:]=self.fftshift_1_floe(value[:,i,:]);
		return fft_speed
		
	def fftshift_1_floe(self,signal) :
		## value=(time,signal.x or signal.y)
		return np.abs(fftshift(fft(signal[:,0]*complex(1,0)+complex(0,1)*signal[:,1])));
		
	def fftshift_speed(self) :
		speed=self.floe_speed
		return self.fftshift_each_floe(speed)
	
	def fftshift_traj(self) :
		traj=self.floe_traj()
		return self.fftshift_each_floe(traj)
	
	def fftshift_pos(self) :
		pos=self.floe_pos
		return self.fftshift_each_floe(pos)
		
	def fftshift_autocorr_floe(self) :
		autocorr_floe=self.autocorr_floe_trajectories();
		return self.fftshift_each_floe(autocorr_floe)		
		
	def fftshift_autocorr_mean(self) :
		autocorr_mean=self.autocorr_trajectories();
		fft_autocorr=np.abs(fftshift(fft(autocorr_mean[:,0]*complex(1,0)+complex(0,1)*autocorr_mean[:,1])))
		return fft_autocorr
		
		
		
	#### PLOTTING
		
	def plot_pos(self,floe_num=0,col='k',lab='pos') :
		plt.plot(self.floe_pos[:,floe_num,0],self.floe_pos[:,floe_num,1],color=col,label=lab)
		plt.xlabel('x')
		plt.ylabel('y')  
	
	def plot_traj(self,floe_num=0,col='k',lab='traj') :
		traj=self.floe_trajectories()
		plt.plot(traj[:,floe_num,0],traj[:,floe_num,1],color=col,label=lab)
		plt.xlabel('x')
		plt.ylabel('y')  
		
	def plot_mass_center_pos(self,col='k',lab='mass center pos') :
		plt.plot(self.mass_center[:,0],self.mass_center[:,1],color=col,label=lab)
		plt.xlabel(r'$x$')
		plt.ylabel(r'$y$') 
	
	def plot_speed(self,floe_num=0,col='k',lab='speed') :
		plt.plot(self.floe_speed[:,floe_num,3],self.floe_speed[:,floe_num,4],color=col,label=lab)
		plt.xlabel('u_x')
		plt.ylabel('u_y')
		
	def plot_fft(self,signal,col='k',lab='fft',omega_log=False) :
		omega=self.frequencies_shifted(signal.shape[0])
		if(omega_log) :
			ind_0=int(np.floor(omega.shape[0]/2))
			plt.plot(omega[(ind_0+1):],signal[(ind_0+1):],color=col,label=lab+"$\omega>0$")
			plt.plot(-omega[0:ind_0],signal[0:ind_0],color=col,linestyle='--',label=lab+" $\omega<0$")
			plt.plot(omega[(ind_0+1)]/10,signal[ind_0],color=col,marker='*')
			plt.xscale('log')
		else :
			plt.plot(omega,signal,color=col,label=lab)
		#plt.yscale('log')
		plt.xlabel(r'$\omega$')

		
	def plot_fft_speed(self,floe_num=0,col='k',lab='fft speed',omega_log=False) :
		signal=self.fftshift_1_floe(self.floe_speed[:,floe_num,:])
		self.plot_fft(signal,col,lab,omega_log) 
		plt.ylabel('$|fft(u_x(t)+iu_y(t))|$') 
		
#	def plot_fft_speed_traj(self,floe_num=0,col='k',lab='fft speed',omega_log=False) :
#		signal=self.fftshift_1_floe(self.floe_speed[:,floe_num,:]);
#		self.plot_fft(signal,col,lab,omega_log) 
#		plt.ylabel('$|fft(u_x(t)+iu_y(t))|$')
		
	def plot_fft_autocorr_mean(self,col='k',lab='fft speed',omega_log=False) :
		signal=self.fftshift_autocorr_mean();
		self.plot_fft(signal,col,lab,omega_log) 
		plt.ylabel('autoccorelation speed trajectories')
		
	def plot_var_distance_from_origin(self,nb_group=1,group_type='no',col=['k','r','b','g'],signal_log=True,time_log=True) :
		var=self.var_distance_from_origin(nb_group,group_type)
		##[floe group, time]
		for i in range(0,nb_group):
			if (signal_log):
				plt.plot(self.time[1:],np.log10(var[i,1:]),color=col[i%len(col)],label='gp num %d'%i)
				plt.ylabel(r'$\log_{10}(<X(t)-X(0)>)(t)$')
			else:
				plt.plot(self.time,var[i,:],color=col[i%len(col)],label='gp num %d'%i)
				plt.ylabel(r'$<X(t)-X(0)>(t)$')
				plt.yscale('log')
		plt.xlabel(r'$t$')
		if(time_log) :
			plt.xscale('log')

	## calculate approxiation of var distance_from_origin in 3 differentes phase 
	## tau=separate time axis in 3 phase [0,tau[0]] [tau[0],tau[1]] [tau[1],end_time]
	def plot_var_distance_from_origin_polyfit(self,col=['k','r','b','g'],tau=[1e3,1e5],signal_log=True,time_log=True) :
		var=self.var_distance_from_origin(1,'no')
		##[floe group, 
		if len(tau)==2:
			phase1=self.time[1:int(tau[0]/self.dt)];
			coeff1=np.polyfit(np.log10(phase1),np.log10(var[0,1:phase1.shape[0]+1]),1);
			phase2=self.time[int(tau[0]/self.dt):int(tau[1]/self.dt)];
			coeff2=np.polyfit(np.log10(phase2),np.log10(var[0,(phase1.shape[0]+1):(phase2.shape[0]+phase1.shape[0]+1)]),1);
			phase3=self.time[int(tau[1]/self.dt):];
			coeff3=np.polyfit(np.log10(phase3),np.log10(var[0,(phase2.shape[0]+phase1.shape[0]+1):]),1);
			
		if (signal_log):
			plt.plot(self.time[1:],np.log10(var[0,1:]),color=col[0],label='var distance');
			if len(tau)==2:
				plt.plot(phase1,coeff1[1]+np.log10(phase1)*coeff1[0],color=col[1],label='$10^{%d}t^{%d}$'%(coeff1[1],coeff1[0]));
				plt.plot(phase2,coeff2[1]+np.log10(phase2)*coeff2[0],color=col[2],label='$10^{%d}t^{%d}$'%(coeff2[1],coeff2[0]));
				plt.plot(phase3,coeff3[1]+np.log10(phase3)*coeff3[0],color=col[3],label='$10^{%d}t^{%d}$'%(coeff3[1],coeff3[0]));
			plt.ylabel(r'$\log_{10}(<X(t)-X(0)>)(t)$')
		else:
			plt.plot(self.time,var[0,:],color=col[0],label='var distance');
			if len(tau)==2:
				plt.plot(phase1,np.power(10.0,coeff1[1])*np.power(phase1,coeff1[0]),color=col[1],label='$e^{%d}t^{%d}$'%(coeff1[1],coeff1[0]));
				plt.plot(phase2,np.power(10.0,coeff1[1])*np.power(phase2,coeff2[0]),color=col[2],label='$e^{%d}t^{%d}$'%(coeff2[1],coeff2[0]));
				plt.plot(phase3,np.power(10.0,coeff1[1])*np.power(phase3,coeff3[0]),color=col[3],label='$e^{%d}t^{%d}$'%(coeff3[1],coeff3[0]));
			plt.ylabel(r'$<X(t)-X(0)>(t)$')
			plt.yscale('log')
		plt.xlabel(r'$t$')
		if(time_log) :
			plt.xscale('log')


	def plot_var_dispersion_pair(self,nb_group=1,group_type='no',col=['k','r','b','g'],signal_log=True,time_log=True) :
		var=self.var_dispersion_pair(nb_group,group_type)
		##[floe group, time]
		for i in range(0,nb_group):
			if (signal_log):
				plt.plot(self.time[1:],np.log10(var[i,1:]),color=col[i%len(col)],label="gp num %d"%i)
				plt.ylabel(r'$\log_{10}(<X(t)-X(0)>)(t)$')
			else:
				plt.plot(self.time,var[i,:],color=col[i%len(col)],label='gp num %d'%i)
				plt.ylabel(r'$<X(t)-X(0)>(t)$')
				plt.yscale('log')
		##print(np.log10(var[i,1:]))
		plt.xlabel(r'$t$')
		if(time_log) :
			plt.xscale('log')
			


	def plot_norm_OBL(self,col='k',lab='OBL speed') :
		OBL=np.sqrt(np.square(self.OBL_speed[:,0])+np.square(self.OBL_speed[:,0]));
		plt.plot(self.time,OBL,color=col,label=lab)
		plt.xlabel(r'$t$')
		plt.ylabel(r'$|OBL|$') 
		

def autocorr(x):
	result = np.correlate(x, x, mode='full')
	return result[result.size//2:]
		


if __name__ == "__main__" :

	myoutput=Data_extractor('./out_uxM19.h5')
	myoutput.print_shape();
	
           
