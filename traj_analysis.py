#!/ usr / bin / env python

import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
from sklearn.decomposition import PCA
from itertools import combinations
import argparse

#Declare arguments
parser = argparse.ArgumentParser(description = 'Determination of DSSP, H-bonds, and PCA for GROMACS Trajectory of PTP1B')
parser.add_argument('-t', required=True, help='File name for input trajectory')
parser.add_argument('-g', required=True, help= 'File name for input topology (gro format)')
parser.add_argument('-f', required=True, help='base for all file names')
parser.add_argument('-a', required=False, default=False, type=bool, help='Is the alpha-7 helix present in the structure')

#Import Arguments
args = parser.parse_args()
File_traj = args.t + '.xtc'
File_gro = args.g + '.gro'
File_base = args.f
a7_present = args.a

#Load trajectories
traj = md.load(File_traj, top=File_gro)
top = traj.topology
#Remove solvent
traj_ns = traj.remove_solvent()
print('Topology Loaded')

#Only do DSSP if a7 is present
if a7_present==True:
    #Compute Secondary Structure
    dssp_list = md.compute_dssp(traj, simplified=False)
    file_dssp = open('DSSP_'+File_base +'.txt','w')
    frame_max,residue=dssp_list.shape
    imax = int(frame_max/10)
    for i in range(imax):
        for j in range(285,295):
            index = i*10
            file_dssp.write(dssp_list[index,j] + ' ')
        file_dssp.write('\n')
    file_dssp.close()
    print('DSSP File Written')
else:
    print('DSSP Skipped')
#H-bond determination
hbonds = md.baker_hubbard(traj, freq=0.6, exclude_water=True, periodic=False)
label = lambda hbond : '%s -- %s' % (traj.topology.atom(hbond[0]), traj.topology.atom(hbond[2]))
file_object = open('Hbonds_'+File_base+'.txt', 'w')
for hbond in hbonds:
    file_object.write(label(hbond))
    file_object.write("\n")
file_object.close()
per = [] #Declare empty array for percentage of time h-bond is formed
da_distances = md.compute_distances(traj, hbonds[:,[1,2]], periodic=False) #Compute distance between h-bond donor and acceptor
da_angles = md.compute_angles(traj, hbonds[:,:], periodic=False) #Compute angle between h-bond donor and acceptor
[num_t, num_h] = np.shape(da_distances)
for j in range(num_h):
    count = 0 #Initialize count
    for i in range(num_t):
        if da_distances[i,j] <= 0.25 and da_angles[i,j] >= 2.094: #If distance between donor and acceptor is less than 2.5A and the angle is greater than 120 degrees or ~ 2.094 radians
            count +=1
    per.append(100*count/num_t) #Percentage of time the h-bond is present
np.savetxt('Hbonds_per_' + File_base + '.txt', per)
print('Hbond Analysis Written')

#Principle Component Analysis
pca = PCA(n_components=10) #initialize
traj_ns.superpose(traj_ns, 0) #align trajectory
reduced_cartesian = pca.fit_transform(traj_ns.xyz.reshape(traj_ns.n_frames, traj_ns.n_atoms * 3))
per_var = np.round(pca.explained_variance_ratio_ * 100, decimals =2)
accum_per_var = [ i for i in [ np . sum ( per_var [: j ]) for j in range (1 ,11)]]
print('PCA Completed')

#Plot PCA
fig = plt.figure()
plt.scatter(reduced_cartesian[:, 0], reduced_cartesian[:,1], marker='x', c=traj.time)
plt.xlabel('PC1')
plt.ylabel('PC2')
plt.title('Cartesian coordinate PCA for'+File_base)
cbar = plt.colorbar()
cbar.set_label('Time [ps]')
fig.savefig('PCA_'+File_base+'.png')


fig.ax = plt.figure(figsize=(10,5.5))

#Scree Plot
labels = [ 'PC' + str(x) for x in range(1 , len(per_var) +1) ]
ax1 = plt.subplot (1 ,2 ,1)
ax2 = ax1.twinx ()
ax1.bar( x = range(1 , len(per_var) +1) , height = per_var , tick_label = labels , alpha = 0.85)
ax2.plot(range(1 , len(per_var) +1), accum_per_var, color = 'r' , marker = 'o')
ax2.grid( True )
xlocs , xlabs = plt.xticks()

ax1.set_ylabel ( 'Percentage of explained variance (%)', size = '12')
ax2.set_ylabel ( 'Accumulated explained variance (%)', size = '12')
ax1.set_xlabel ( 'Principal Components' , size = '12')
ax1.set_ylim ([0 , max(per_var)*1.1 ] )
ax2.set_ylim ([0 , max(accum_per_var)*1.1 ] )
plt.title ( 'Scree Plot' , size = '14')
plt.grid ( True )
plt.savefig('Scree_'+File_base+'.png')
