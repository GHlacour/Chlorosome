# -*- coding: utf-8 -*-
"""
Created on Fri Mar 27 13:46:01 2020

@author: 31641
"""


import numpy as np
import matplotlib.pyplot as plt
import MDAnalysis as mda
import random


#Loading topology (.gro) and trajectory file (.xtc)
u=mda.Universe ('confout.gro','traj_comp.xtc')
print (u)
print ('Number of frames', len(u.trajectory))
#for ts in u.trajectory:
#  print("Frame: {0:5d}, Time: {1:8.3f} ps".format(ts.frame, u.trajectory.time))

dynamics_of_H_bonding=[]

for ts in u.trajectory[0:501:1]:
    print ('Frame number', ts)

    #Selecting donor's O atom:
    Donor=u.atoms.select_atoms ('name OH')
    #Positions of donor's O atom:
    Donor_positions=Donor.positions
    #print ('Donor positions', Donor_positions)

    #Determining number of BChl molecules = sites
    number_of_sites=int(Donor_positions.size/3)
    #print ('Number of sites is', number_of_sites)

    #Determining number of interactions:
    number_of_interactions=int((number_of_sites*number_of_sites-number_of_sites)/2)

    
    #Selecting position of H atom:
    H=u.atoms.select_atoms ('name HO')
    #Position of H atom:
    H_positions=H.positions 
    #print ('H_positions', H_positions)

    #Selecting acceptor: 
    A=u.atoms.select_atoms ('name O_2')
    #Positions of O_2 groups:
    A_positions=A.positions
    Acceptor_positions=np.zeros((number_of_sites, 3))
    for i in range (number_of_sites):
        Acceptor_positions[i,:]=A_positions[i*2]
       # print ('Acceptor_positions', Acceptor_positions)


    #for i in range (number_of_sites-1):
    counter_counter =np.array (number_of_sites)
#Random sampling of molecules:
    samples=random.sample(range (number_of_sites), k=10)
    for i in samples:
        counter_donor=0    
        counter_acceptor=0
        counter=0
        print (i)    
            #Looking at the distance between the O and H on the chosen molecules
        r_d_H = Donor_positions [i] - H_positions [i]
        r_d_H_intensity  = np.linalg.norm (r_d_H)
        for j in range (0,number_of_sites-1):
                
                    #Distance between the O and H atom of the chosen molecules
            distance_d_H= Donor_positions[j]- H_positions[j]
            distance_d_H_intensity=np.linalg.norm (distance_d_H )
                    
                    #Molecule i as donor check
            r_d_a = Donor_positions[i]- Acceptor_positions[j]
            
            r_d_a_intensity  = np.linalg.norm (r_d_a)
                    
                    #Molecule i as acceptor check
            r_a_d =Donor_positions [j] - Acceptor_positions [i]
            r_a_d_intensity  = np.linalg.norm (r_a_d )          
            #print (r_a_d_intensity)
                    
                #Is distance between groups according to criterium (CHECKING RADIUS CRITERIUM)
            if r_d_a_intensity  <= 3.5:
                        #If distance is according to criterium, check the other condition: (CHECKING ANGLE CRITERIUM)
                alpha_rad=np.arccos(np.dot (r_d_a, r_d_H )/(r_d_a_intensity *r_d_H_intensity))
                alpha=alpha_rad *180/np.pi
                        
                if alpha <= 30:
                    counter_donor=counter_donor+1
                    print ('Counter donor', counter_donor, i)
                        #Checking is this molecule also an acceptor for hydrogen bonding:
            if r_a_d_intensity <= 3.5:
                beta_rad = np.arccos (np.dot (r_a_d, distance_d_H )/(r_a_d_intensity*distance_d_H_intensity))
                beta = beta_rad * 180/np.pi
                if beta <= 30:
                    counter_acceptor=counter_acceptor+1
                    print ('Counter acceptor', counter_acceptor, i)
    
         
          
