#!/ usr / bin / env python

import mdtraj as md
import numpy as np
from sklearn.decomposition import PCA
from itertools import combinations
import argparse
from itertools import product

#Declare arguments
parser = argparse.ArgumentParser(description = 'Determination the number of h-bonds at every time point b/w res 195-200 and residues outside the a3 helix')
parser.add_argument('-t', required=True, help='File name for input trajectory')
parser.add_argument('-g', required=True, help= 'File name for input topology (gro format)')

#Import Arguments
args = parser.parse_args()
File_traj = args.t + '.xtc'
File_gro = args.g + '.gro'

#Load trajectories
traj = md.load(File_traj, top=File_gro)
top = traj.topology

#Process Trajectory
traj_ns = traj.remove_solvent() #Remove solvent from the trajectory leaving only protein (and ligand if applicable)
print('Trajectory Loaded')

#Determine list of H-bonds present in the trajectory for over 60% of the frames
hbonds = md.baker_hubbard(traj_ns, freq=0.2, exclude_water=True, periodic=False)
label = lambda hbond : '%s -- %s' % (traj_ns.topology.atom(hbond[0]), traj_ns.topology.atom(hbond[2])) #Extract labels for h-bonds

#Set atom indices for the start and end of a3 helix and a3 helix section
if traj_ns.n_residues == 299:
    range_a3 = top.select('185 <= resid and resid <= 210') #a3 helix an immediate surrounding resid 186 to 206
    range_a3_sect = top.select('197 <= resid and resid <= 202') #resid 198 to 202
else:
    range_a3 = top.select('184 <= resid and resid <= 204')
    range_a3_sect = top.select('196 <= resid and resid <= 201')
a3_start = range_a3[0]
a3_end = range_a3[-1]
a3_sect_start = range_a3_sect[0]
a3_sect_end = range_a3_sect[-1]

#Compute distances and angles for all bonds
da_distances = md.compute_distances(traj_ns, hbonds[:,[1,2]], periodic=False) #Compute distance between h-bond donor and acceptor
da_angles = md.compute_angles(traj_ns, hbonds[:,:], periodic=False) #Compute angle between h-bond donor and acceptor
[num_t, num_h] = np.shape(da_distances) #save values for number of frames(num_t) and number of bonds(num_b) to caculate

#Empty array for number of h-bonds fitting criteria over time
Hbond_num = np.zeros(num_t)
Hbond_name = []
file_name = open('Hbonds_bw_a3_sect_and_not_a3_name.txt', 'w')

#Determine the number of bonds occuring between the atoms in the a3 helix and those outside it
for i in range(num_t): #Loop through all frames
    num = 0
    for j in range(num_h): #Loop through all bonds
        if da_distances[i,j] <= 0.25 and da_angles[i,j] >= 2.094: #If bond is present in this frame
            if (hbonds[j,1] <= a3_sect_end and hbonds[j,1] >= a3_sect_start and (hbonds[j,2] > a3_end or hbonds[j,2] < a3_start)) or (hbonds[j,2] <= a3_sect_end and hbonds[j,2] >= a3_sect_start and (hbonds[j,1] > a3_end or hbonds[j,1] < a3_start)):#If bond is between the alpha 3 helix residues 195-200 and residues outide the a3 helix
                num += 1
                name = label(hbonds[j])
                if name not in Hbond_name:
                    Hbond_name.append(name)
                    file_name.write(str(name) + '\n')
    Hbond_num[i] = num

np.savetxt('Hbonds_bw_a3_sect_and_not_a3.txt', Hbond_num)


