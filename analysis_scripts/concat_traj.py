#!/ usr / bin / env python
import mdtraj as md
import argparse
import numpy as np
import sys

#Declare arguments
parser = argparse.ArgumentParser(description = 'Concatenate uncorrelated frames from multiple trajectories')
parser.add_argument('-d', required=True, type=str, help='Directory Path for PTP1B repository')
parser.add_argument('-l', required=True, type=str, help='Ligand Name')

#Import Arguments
args = parser.parse_args()
directory = args.d
lig = args.l

#Import custom modules
sys.path.insert(1, directory + '/util/')
import mdfunc

#Load trajectory file names
traj_gro = open('traj_gro_' + lig + '_file.txt', 'r').readlines()[0].split(' ')
traj_xtc = open('traj_xtc_' + lig + '_file.txt', 'r').readlines()[0].split(' ')
uncorr = open('uncorr_file_' + lig + '.txt', 'r').readlines()[0].split(' ')

#Load and concatenate trajectories
for i in range(len(traj_gro)):
    #Load files
    File_gro = traj_gro[i].strip()
    File_traj = traj_xtc[i].strip()
    File_uncorr = uncorr[i].strip()

    #Load trajectory
    traj = md.load(File_traj, top=File_gro)
    top = traj.topology
    if i != 2:
        traj_ns = traj.atom_slice(top.select('protein or resname AD or resname BBR')) #Select only atoms in the protein or ligand
    else:
        traj_ns = traj.atom_slice(top.select('(protein or resname AD or resname BBR) and not resid 0')) #Select only atoms in the protein or ligand

    #Limit trajectory to uncorrelated frames
    uncorr_ind_string = open(File_uncorr, 'r').readlines()
    uncorr_ind = np.zeros(len(uncorr_ind_string), dtype=int)
    for j in range(len(uncorr_ind_string)):
        uncorr_ind[j] = int(j)
    traj_uncorr = traj_ns.slice(uncorr_ind)

    print('Trajectory ' + str(i) + ' loaded')
    
    #Export full trajectory
    traj_uncorr.save_xtc('Full_ligand_traj' + str(i) + '.xtc')

