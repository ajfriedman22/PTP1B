#!/ usr / bin / env python
import mdtraj as md
import numpy as np
import argparse

#Declare arguments
parser = argparse.ArgumentParser(description = 'Determination of Percentage of time each individual h-bond is formed and for h-bonds Formed simultaneously')
parser.add_argument('-t', required=True, help='File name for input trajectory')
parser.add_argument('-g', required=True, help= 'File name for input topology (gro format)')
parser.add_argument('-p', required=True, help= 'File path for H-bonds of interest')
parser.add_argument('-d', required=False, default='none', help= 'File path for other selected h-bonds')

#Import Arguments
args = parser.parse_args()
File_traj = args.t + '.xtc'
File_gro = args.g + '.gro'
File_path = args.p
directory = args.d

#Load trajectories
traj = md.load_xtc(File_traj, top=File_gro)
top = traj.topology

traj_ns = traj.remove_solvent() #Remove solvent from the trajectory leaving only protein (and ligand if applicable)
traj_prot = traj.atom_slice(top.select('protein')) #Select only atoms in the protein
print('Trajectory Loaded')

#Make array for all atom elements
if directory != 'none':
    atom_note = open(directory + '/Hbonds_note_atom.txt', 'r').readlines()
if traj_prot.n_residues == 299 or traj_prot.n_residues == 287: #With or witout the a7 helix
    atom_bonds = open(File_path + 'Atom_uncomm.txt', 'r').readlines()

else: # Atom numbers if the first residue is missing
    atom_bonds = open(File_path + 'Atom_uncomm_1sug.txt', 'r').readlines()

#Make array for bond names
name_bonds = open(File_path + 'Hbond_uncommon.txt', 'r').readlines()
if directory != 'none':
    name_note = open(directory + 'Hbonds_note_atom.txt', 'r').readlines()

#Output file for bond 
output_single = open('Hbond_per_single.txt', 'w') #percentage of time each bond occurs
output_bond = open('Hbond_per_comb.txt', 'w') #Percentage of time two bonds occur
output_pos_bond = open('Hbond_pos_corr.txt', 'w')
output_neg_bond = open('Hbond_neg_corr.txt', 'w')
if directory != 'none':
    output_note = open('Hbond_lit_per.txt', 'w')

#Declare array for bond correlations
num_bonds = len(atom_bonds)
bond_comb_frac = np.zeros((num_bonds, num_bonds))
bond_corr = bond_comb_frac
bond_single_frac = np.zeros((num_bonds))

#track bond indicies
num_i = 0
#Determine the percent of time each bond combination is present
for i in atom_bonds:
    #track bond indicies
    num_j = 0
    for j in atom_bonds:
        if i != j:
            #Seperate atom indices into float format
            bond1 = i.split()
            bond2 = j.split()
            bond_d = np.array([[float(bond1[1]), float(bond1[2])] , [float(bond2[1]), float(bond2[2])]])
            bond_a = np.array([[float(bond1[0]), float(bond1[1]), float(bond1[2])],[float(bond2[0]), float(bond2[1]), float(bond2[2])]])
            #Compute distances and angles over time for both bonds
            dist = md.compute_distances(traj_ns, bond_d, periodic = False)
            angle = md.compute_angles(traj_ns, bond_a , periodic = False)
            #Determine the percent of time both bonds are formed
            count=0
            count_s=0 #single bonds
            for k in range(len(dist)):
                if dist[k][0] <= 0.25 and dist[k][1] <= 0.25 and angle[k][0] >= 2.094 and angle[k][1] >= 2.094:
                    count += 1
                if dist[k][0] <= 0.25 and angle[k][0] >= 2.094:
                    count_s += 1
            bond_comb_frac[num_i][num_j] = count/len(dist)
            output_bond.write(name_bonds[num_i] + name_bonds[num_j] + str(100*count/len(dist)) + '\n')
            if num_j == 0 or (num_i == 0 and num_j ==1):
                output_single.write(name_bonds[num_i] + str(100*count_s/len(dist)) + '\n')
                bond_single_frac[num_i] = count_s/len(dist)
        num_j += 1
    num_i += 1

#Determine bond correlation
num_i = 0
for i in atom_bonds:
    num_j = 0
    for j in atom_bonds:
        if i != j:
            if bond_single_frac[num_i] != 0 and bond_single_frac[num_j] != 0:
                bond_corr[num_i][num_j] = (bond_comb_frac[num_i][num_j] - (bond_single_frac[num_i]*bond_single_frac[num_j]))/(bond_single_frac[num_i]*bond_single_frac[num_j])
            if bond_single_frac[num_i] == 0 or bond_single_frac[num_j] == 0:
                bond_corr[num_i][num_j] = 0
            if bond_corr[num_i][num_j] >= 0.1:
                output_pos_bond.write(name_bonds[num_i] + '+' + name_bonds[num_j] + ' Corr = ' + str(bond_corr[num_i][num_j]) + '\n')
            if bond_corr[num_i][num_j] <= -0.1:
                output_neg_bond.write(name_bonds[num_i] + '+' + name_bonds[num_j] + ' Corr = ' + str(bond_corr[num_i][num_j]) + '\n')
        num_j += 1
    num_i += 1

if directory != 'none':
    #Determine the percent of time each bond combination is present
    num_i = 0
    for i in atom_note:
        #track bond indicies
        num_j = 0
        for j in atom_note:
            if i != j:
                #Seperate atom indices into float format
                bond1 = i.split()
                bond2 = j.split()
                if traj_prot.n_residues == 299 or traj_prot.n_residues == 287:
                    bond_d = np.array([[float(bond1[1]), float(bond1[2])] , [float(bond2[1]), float(bond2[2])]])
                    bond_a = np.array([[float(bond1[0]), float(bond1[1]), float(bond1[2])],[float(bond2[0]), float(bond2[1]), float(bond2[2])]])
                else:
                    bond_d = np.array([[float(bond1[1])-17, float(bond1[2])-17] , [float(bond2[1])-17, float(bond2[2])-17]])
                    bond_a = np.array([[float(bond1[0])-17, float(bond1[1])-17, float(bond1[2])-17],[float(bond2[0])-17, float(bond2[1])-17, float(bond2[2])-17]])

                #Compute distances and angles over time for both bonds
                dist = md.compute_distances(traj_ns, bond_d, periodic = False)
                angle = md.compute_angles(traj_ns, bond_a , periodic = False)
                #Determine the percent of time both bonds are formed
                count=0
                count_s=0 #single bonds
                for k in range(len(dist)):
                    if dist[k][0] <= 0.25 and dist[k][1] <= 0.25 and angle[k][0] >= 2.094 and angle[k][1] >= 2.094:
                        count += 1
                    if dist[k][0] <= 0.25 and angle[k][0] >= 2.094:
                        count_s += 1
                bond_comb_frac[num_i][num_j] = count/len(dist)
                output_bond.write(name_bonds[num_i] + name_bonds[num_j] + str(100*count/len(dist)) + '\n')
                if num_j == 0 or (num_i == 0 and num_j ==1):
                    output_note.write(name_bonds[num_i] + str(100*count_s/len(dist)) + '\n')
                    bond_single_frac[num_i] = count_s/len(dist)
            num_j += 1
        num_i += 1

print('Hbond Percentages Calculated')
