# -*- coding: utf-8 -*-
"""
Created on Mon Apr 13 17:11:17 2020

@author: 31641
"""

import numpy as np
import MDAnalysis as mda
import matplotlib.pyplot as plt
import math
 

#Change index for molecules 


#Loading topology (.gro) and trajectory file (.xtc)
u=mda.Universe ('confout.gro','traj_comp.xtc')
print (u)
l=len(u.trajectory)
print ('Number of frames', l)
TDM=np.zeros((l,3))
#for ts in u.trajectory:
#  print("Frame: {0:5d}, Time: {1:8.3f} ps".format(ts.frame, u.trajectory.time))

for ts in u.trajectory[0:501:1]:
    print ('Frame number', ts)
#Selecting Mg atoms (BChl)    
    Mg=u.atoms.select_atoms ('name MG')
#Extracting positions of Mg atoms
    Mg_positions=Mg.positions
#print ('Positions of all Mg atoms', Mg_positions)
#Determining number of BChl molecules = sites
    number_of_sites=int(Mg_positions.size/3)
    # print ('Number of sites is', number_of_sites)

   

#Selecting N atoms 
    N=u.atoms.select_atoms ('name NA')
#Extracting positions of N atoms
    N_positions=N.positions 
    #print ('Positions of all N atoms',N_positions)

    #Positions of NA and NC atoms (BChl) 
    NA_positions=np.zeros((number_of_sites,3))
    NC_positions=np.zeros((number_of_sites,3))
    for i in range (number_of_sites):
        NA_positions[i,:]=N_positions[4*i]
        NC_positions[i,:]=N_positions[4*i+2]
#print ('Positions of NA atoms', NA_positions)
#print ('Positions of NC atoms', NC_positions)

#Convention:y axis of the molecules is defined as axis passing through the N atoms of rings A and C
#Qy transition dipole moment: magnitude and direction information
    NC_NA_distance = NC_positions - NA_positions
#print ('Distances between the NA and NC', NC_NA_distance)

    tdm_length=5.48 #Debye
    i=2558
    unit_transition_dipole=tdm_length*NC_NA_distance[i]/np.linalg.norm(NC_NA_distance[i])
    TDM [ts.frame,:] =unit_transition_dipole


orientation_correlation_function =np.zeros ((l,l))
angles=np.zeros((l,l))
for i in range (l):
    p=np.zeros(l) 
    theta=np.zeros(l)
    for j in range (0,l-i):
            theta[j]=np.arccos(np.dot(TDM[i],TDM[j])/(np.linalg.norm (TDM[i])*np.linalg.norm (TDM[j])))
            #since I observed appearance of the nun values in matrix of angles I tried to fix them in this way
            if np.isnan(theta[j])=='True':
                alpha=0
                p[j]=(1/(l+1-j))*(1/5)*(3*(np.cos(alpha)**2)-1)
            else:
                
                
                p[j]=(1/(l+1-j))*(1/5)*(3*(np.cos(theta[j])**2)-1)

    angles[i,:]=theta
    orientation_correlation_function [i,:]=p


#            #print (orientation_correlation_function[i])
print (angles)
orientation_correlation_function_all=np.sum (orientation_correlation_function,axis=0)
print (orientation_correlation_function_all)


        


x=np.linspace (0, 10, len (u.trajectory))          
plt.plot (x, orientation_correlation_function_all,'o')
# plt.xlim (0,2)
# plt.ylim (0.370)
plt.xlabel ('t(ps)')
plt.ylabel('Orientation correlation function')
plt.draw()
            

    