# -*- coding: utf-8 -*-
"""
Created on Thu Apr  2 14:38:52 2020

@author: 31641
"""


import numpy as np
import MDAnalysis as mda



#Loading topology (.gro) and trajectory file (.xtc)
u=mda.Universe ('confout.gro','traj_comp.xtc')
print (u)
print ('Number of frames', len(u.trajectory))
#for ts in u.trajectory:
#  print("Frame: {0:5d}, Time: {1:8.3f} ps".format(ts.frame, u.trajectory.time))

dynamics_of_H_bonding=[]

for ts in u.trajectory[0:1:1]:
    print ('Frame number', ts)
            #Selecting donor's O atom:
    Donor=u.atoms.select_atoms ('name OH')
    #Positions of donor's O atom:
    Donor_positions=Donor.positions
    print ('Donor positions', Donor_positions)
    
    #Determining number of BChl molecules = sites
    number_of_sites=int(Donor_positions.size/3)
    print ('Number of sites is', number_of_sites)
    
    #Determining number of interactions:
    number_of_interactions=int((number_of_sites*number_of_sites-number_of_sites)/2)
    
        
    #Selecting position of H atom:
    H=u.atoms.select_atoms ('name HO')
    #Position of H atom:
    H_positions=H.positions 
    print ('H_positions', H_positions)
    
    #Selecting acceptor: 
    A=u.atoms.select_atoms ('name O_2')
    #Positions of O_2 groups:
    A_positions=A.positions
    Acceptor_positions=np.zeros((number_of_sites, 3))
    for i in range (number_of_sites):
        Acceptor_positions[i,:]=A_positions[i*2]
    print ('Acceptor_positions', Acceptor_positions)
    
    acceptor_index=[]
    donor_index=[]
    
                      
    for i in range (number_of_sites-1):
            #Determining distance between the position of donor's oxygen and H atom through which the bond is created 
            r_d_H  = Donor_positions [i] - H_positions [i]
            r_d_H_intensity  = np.linalg.norm (r_d_H)
            for j in range (0,number_of_sites-1):
                #Probing molecule as a donor for hidrogen bonding:
                    r_d_a= Donor_positions[i]- Acceptor_positions[j]
                    r_d_a_intensity  = np.linalg.norm (r_d_a)
                #Is distance between groups according to criterium
                    if r_d_a_intensity <= 3.5:
                        numerator=np.dot (r_d_H, r_d_a)
                        denumerator=r_d_H_intensity*r_d_a_intensity
                        alpha_rad=np.arccos(numerator/denumerator)
                        alpha=alpha_rad*180/np.pi
                        if alpha <= 30:
                            donor_index.append (i)
                            acceptor_index.append (j)
                            #print ('angle criterium checked', i, j)
    acceptor=np.asarray (acceptor_index)
    donor=np.asarray (donor_index)
    
    acc=set (acceptor_index)
    don=set (donor_index)
    
    acc_and_don= acc.intersection (don)
    number_of_molecules_with_2_H_bonds=len(acc_and_don)
    
    print ('These molecule have 2 H bonds', acc_and_don, 'and their number is', number_of_molecules_with_2_H_bonds)
    print ('Number of acceptors', acceptor.size, 'Contribution acceptors',acceptor.size/number_of_sites)
    print ('Number of donors', donor.size, 'Contribution of donors', donor.size/number_of_sites)
    
    
    #np.savetxt ('akceptori.txt', acceptor)
    #np.savetxt ('donori.txt', donor)