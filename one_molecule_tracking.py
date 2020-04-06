#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 21 16:22:16 2020

@author: vesnaeric
"""

import numpy as np
import matplotlib.pyplot as plt
import MDAnalysis as mda
from scipy.stats import norm

Molecule_6000=[]
#Loading topology (.gro) and trajectory file (.xtc)
u=mda.Universe ('mix_nvt.gro','mix_nvt.xtc')
print (u)

print ('Number of frames', len(u.trajectory))
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
    print ('Number of sites is', number_of_sites)

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
    unit_transition_dipole =np.zeros ((number_of_sites,3))
    for i in range (number_of_sites):
        unit_transition_dipole [i,:]=tdm_length*NC_NA_distance[i]/np.linalg.norm(NC_NA_distance[i])
#print ('Unit transition dipole momemnt of Qy transition', unit_transition_dipole)

#Building Hamiltonian for this system
    site_energy= 0 #cm-1
    Hamiltonian= np.eye(number_of_sites,number_of_sites)*site_energy
    for i in range (number_of_sites-1):
        for j in range (i+1, number_of_sites):
            tdi=unit_transition_dipole [i]
            tdj=unit_transition_dipole [j]
            r_ij=Mg_positions [i]-Mg_positions [j] #distance between sites
            r_ij_length=np.linalg.norm(r_ij)
            r_ij_length3=r_ij_length*r_ij_length*r_ij_length
            r_ij_length5=r_ij_length3*r_ij_length*r_ij_length
            J_ij=((np.dot(tdi,tdj))/(r_ij_length3)-3*(np.dot(tdi,r_ij)*np.dot(tdj,r_ij)/(r_ij_length5)))*5.04*1000 #5.04 unit conversion factor #1000 converting to nm
            Hamiltonian [i][j]= Hamiltonian [j][i] =J_ij
#print ('System Hamiltonian:', Hamiltonian)
    Coupling_strength_S=0
    for i in range (number_of_sites):
        if i==6000:
            Coupling_strength_S=Coupling_strength_S + np.sum (Hamiltonian [i])
            print ('Coupling strength', Coupling_strength_S )
            Molecule_6000.append (Coupling_strength_S)

ax = plt.subplot(111)
ax.plot(np.linspace (0, 10, len (u.trajectory)), Molecule_6000, 'r--', lw=2, label=r"$R_G$")
ax.set_xlabel("time (ps)")
ax.set_ylabel(r"coupling strength of molecule with index 6000 $R_G$ ($\AA$)")
ax.figure.savefig("tracking molecule 6000.pdf")
plt.draw()



