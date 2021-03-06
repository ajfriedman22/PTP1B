import mdtraj as md
import numpy as np
import argparse

parser = argparse.ArgumentParser(description = 'Determination of percent each h-bond within a list was formed for GROMACS Trajectory')
parser.add_argument('-t', required=True, help='File name for input trajectory')
parser.add_argument('-g', required=True, help= 'File name for input topology (gro format)')
parser.add_argument('-f', required=True, help= 'Distinction for group of hbonds of interest')
parser.add_argument('-p', required=True, help= 'File path for list of h-bonds of interest')

#Parse arguments
args = parser.parse_args()
File_traj = args.t + '.xtc'
File_gro = args.g + '.gro'
File_path = args.p
File_name = args.f

#Load trajectories
traj = md.load(File_traj, top=File_gro)
top = traj.topology

#Process Trajectory
traj_ns = traj.remove_solvent() #Remove solvent from the trajectory leaving only protein (and ligand if applicable)
print('Topology Loaded')

#Load atoms and names corresponding to h-bonds of note
bond_name = open(File_path + 'Hbonds_uncommon_' + File_name + '.txt', 'r').readlines()
bond_atom = open(File_path + 'Hbonds_uncommon_atom_' + File_name + '.txt', 'r').readlines()

#Convert atoms list to array
hbonds = np.zeros((len(bond_atom),3))
j = 0
for i in bond_atom:
    split_i = i.split(' ')
    for n in range(3):
        hbonds[j][n] = float(split_i[n])
    j += 1

#Determine the exact percentage of time that each h-bond of note is present
per = [] #Declare empty array for percentage of time h-bond is formed
da_distances = md.compute_distances(traj_ns, hbonds[:,[1,2]], periodic=False) #Compute distance between h-bond donor and acceptor
da_angles = md.compute_angles(traj_ns, hbonds[:,:], periodic=False) #Compute angle between h-bond donor and acceptor
[num_t, num_h] = np.shape(da_distances) #save values for number of frames(num_t) and number of bonds(num_b) to caculate
for j in range(num_h): #Loop through all h-bonds
    count = 0 #Initialize count
    for i in range(num_t): #Loop through all frames
        if da_distances[i,j] <= 0.25 and da_angles[i,j] >= 2.094: #If distance between donor and acceptor is less than 2.5A and the angle is greater than 120 degrees or ~ 2.094 radians
            count +=1
    per.append(100*count/num_t) #Percentage of time each h-bond is present in trajectory
np.savetxt('Hbonds_per_' + File_name + '_uncommon.txt',per)

print('Analysis Complete')
