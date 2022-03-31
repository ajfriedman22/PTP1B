#!/ usr / bin / env python

import mdtraj as md
import numpy as np
from sklearn.decomposition import PCA
from itertools import combinations
import argparse
from itertools import product

#Declare arguments
parser = argparse.ArgumentParser(description = 'Determination the number of h-bonds at every time point b/w section A and residues outside section B')
parser.add_argument('-t', required=True, help='File name for input trajectory')
parser.add_argument('-g', required=True, help= 'File name for input topology (gro format)')
parser.add_argument('-A', required=True, help= 'Residue range for section A (input format:# #)')
parser.add_argument('-B', required=True, help= 'Residue range for section B (input format:# #)')
parser.add_argument('-f', required=True, help= 'File output identifier')


#Import Arguments
args = parser.parse_args()
File_traj = args.t + '.xtc'
File_gro = args.g + '.gro'
File_name = args.f
res_sect_A = args.A
res_sect_B = args.B

#Load trajectories
traj = md.load(File_traj, top=File_gro)
top = traj.topology

#Process Trajectory
traj_ns = traj.remove_solvent() #Remove solvent from the trajectory leaving only protein (and ligand if applicable)
traj_prot = traj.atom_slice(top.select('protein')) #Select only atoms in the protein
print('Trajectory Loaded')

#Determine list of H-bonds present in the trajectory for over 60% of the frames
hbonds = md.baker_hubbard(traj_ns, freq=0.2, exclude_water=True, periodic=False)
label = lambda hbond : '%s -- %s' % (traj_ns.topology.atom(hbond[0]), traj_ns.topology.atom(hbond[2])) #Extract labels for h-bonds

#Seperate input into residues
res_sect_A_int = map(int, res_sect_A.split(' '))
res_sect_B_int = map(int, res_sect_B.split(' '))

#Set atom indices for the start and end of section A and B
if traj_prot.n_residues == 299:
    range_A = top.select(str(res_sect_A_int[0]) + ' <= resid and resid <= ' + str(res_sect_A_int[1])) #Seperate section A
    range_B = top.select(str(res_sect_B_int[0]) + ' <= resid and resid <= ' + str(res_sect_B_int[1])) #seperate section B
else:
    range_A = top.select(str(res_sect_A_int[0] - 1) + ' <= resid and resid <= ' + str(res_sect_A_int[1]) - 1) #Seperate section A
    range_B = top.select(str(res_sect_B_int[0] - 1) + ' <= resid and resid <= ' + str(res_sect_B_int[1]) - 1) #seperate section B

A_start = range_A[0]
A_end = range_A[-1]
B_start = range_B[0]
B_end = range_B[-1]

#Compute distances and angles for all bonds
da_distances = md.compute_distances(traj_ns, hbonds[:,[1,2]], periodic=False) #Compute distance between h-bond donor and acceptor
da_angles = md.compute_angles(traj_ns, hbonds[:,:], periodic=False) #Compute angle between h-bond donor and acceptor
[num_t, num_h] = np.shape(da_distances) #save values for number of frames(num_t) and number of bonds(num_b) to caculate

#Empty array for number of h-bonds fitting criteria over time
Hbond_num = np.zeros(num_t)
Hbond_name = []
file_name = open('Hbonds_bw_' + File_name + '_name.txt', 'w')

#Determine the number of bonds occuring between the atoms in the a3 helix and those outside it
for i in range(num_t): #Loop through all frames
    num = 0
    for j in range(num_h): #Loop through all bonds
        if da_distances[i,j] <= 0.25 and da_angles[i,j] >= 2.094: #If bond is present in this frame
            if (hbonds[j,1] <= A_end and hbonds[j,1] >= A_start and (hbonds[j,2] > B_end or hbonds[j,2] < B_start)) or (hbonds[j,2] <= A_end and hbonds[j,2] >= A_start and (hbonds[j,1] > B_end or hbonds[j,1] < B_start)): #If bond is between residues in section A and outside section B
                num += 1
                name = label(hbonds[j])
                if name not in Hbond_name:
                    Hbond_name.append(name)
                    file_name.write(str(name) + '\n')
    Hbond_num[i] = num

np.savetxt('Hbonds_bw_' + File_name + '.txt', Hbond_num)

