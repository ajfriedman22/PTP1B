#!/ usr / bin / env python

import mdtraj as md
import numpy as np
import argparse

#Declare arguments
parser = argparse.ArgumentParser(description = 'Determination of Percentages for h-bonds Formed Together')
parser.add_argument('-t', required=True, help='File name for input trajectory')
parser.add_argument('-g', required=True, help= 'File name for input topology (gro format)')
parser.add_argument('-w', required=True, help='File name for WPD loop distances')
parser.add_argument('-s', required=False, type = bool, default = False, help='True if the base was 1sug or False if not')
parser.add_argument('-f', required=False, type = bool, default = False, help='True if the files are in a sub-folder')

#Import Arguments
args = parser.parse_args()
File_traj = args.t + '.xtc'
File_gro = args.g + '.gro'
File_WPD = args.w + '.xvg'
lig = 'AD_all'
check_s = args.s
check_f = args.f

#Load trajectories
traj = md.load_xtc(File_traj, top=File_gro)
print('Trajectory Loaded')

#Make array for all atom elements
if check_s == True and check_f == False:
    atom_bonds = open('../compare_' + lig + '/Hbond/Atom_uncomm_1sug.txt', 'r').readlines()
if check_s == False and check_f == False:
    atom_bonds = open('../compare_' + lig +'/Hbond/Atom_uncomm.txt', 'r').readlines()
if check_s == True and check_f == True:
    atom_bonds = open('../../compare_' + lig +'/Hbond/Atom_uncomm_1sug.txt', 'r').readlines()
if check_s == False and check_f == True:
    atom_bonds = open('../../compare_' + lig +'/Hbond/Atom_uncomm.txt', 'r').readlines()

#Make array for bond names
if check_f == False:
    name_bonds = open('../compare_' + lig +'/Hbond/Hbond_uncommon.txt', 'r').readlines()
else:
    name_bonds = open('../../compare_' + lig +'/Hbond/Hbond_uncommon.txt', 'r').readlines()

#Output file for bond 
output_single = open('Hbond_per_single.txt', 'w') #percentage of time each bond occurs
output_bond = open('Hbond_per_comb.txt', 'w') #Percentage of time two bonds occur
output_pos_bond = open('Hbond_pos_corr.txt', 'w')
output_neg_bond = open('Hbond_neg_corr.txt', 'w')

#Output files for correlation to WPD loop open or closed
output_pos_WPD = open('Hbond_WPD_open.txt', 'w')
output_neg_WPD = open('Hbond_WPD_close.txt', 'w')

#Get WPD distances
t, w = [], []
with open(File_WPD) as f:
    for _ in range(17):
        next(f)
    for line in f:
        cols = line.split()
        if len(cols) == 2:
            t.append(float(cols[0]))
            w.append(float(cols[1]))
count_WPD = 0
for i in range(len(t)):
    if w[i] >= 1:
        count_WPD += 1
frac_WPD = count_WPD/len(t)
print(frac_WPD)

#Declare array for bond correlations
num_bonds = len(atom_bonds)
bond_comb_frac = np.zeros((num_bonds, num_bonds))
bond_corr = bond_comb_frac
bond_single_frac = np.zeros((num_bonds))
WPD_hbond_frac = np.zeros((num_bonds))

#track bond indicies
num_i = 0
num_j = 0
#Determine the percent of time each bond combination is present
for i in atom_bonds:
    num_j = 0
    for j in atom_bonds:
        if i != j:
            #Seperate atom indices into float format
            bond1 = i.split()
            bond2 = j.split()
            bond_d = np.array([[float(bond1[1]), float(bond1[2])] , [float(bond2[1]), float(bond2[2])]])
            bond_a = np.array([[float(bond1[0]), float(bond1[1]), float(bond1[2])],[float(bond2[0]), float(bond2[1]), float(bond2[2])]])
            #Compute distances and angles over time for both bonds
            dist = md.compute_distances(traj, bond_d, periodic = False)
            angle = md.compute_angles(traj, bond_a , periodic = False)
            #Determine the percent of time both bonds are formed
            count=0
            count_s=0 #single bonds
            count_w = 0 #WPD open + h-bond present
            for k in range(len(dist)):
                if dist[k][0] <= 0.25 and dist[k][1] <= 0.25 and angle[k][0] >= 2.094 and angle[k][1] >= 2.094:
                    count += 1
                if dist[k][0] <= 0.25 and angle[k][0] >= 2.094:
                    count_s += 1
                if dist[k][0] <= 0.25 and angle[k][0] >= 2.094 and w[k] >= 1:
                    count_w += 1
            WPD_hbond_frac[num_i] = count_w/len(dist)
            bond_comb_frac[num_i][num_j] = count/len(dist)
            output_bond.write(name_bonds[num_i] + name_bonds[num_j] + str(100*count/len(dist)))
            output_bond.write('\n')
            if num_j == 0 or (num_i == 0 and num_j ==1):
                output_single.write(name_bonds[num_i] + str(100*count_s/len(dist)) + '\n')
                bond_single_frac[num_i] = count_s/len(dist)
        num_j += 1
    num_i += 1

#Determine bond correlation
num_i = 0
for i in atom_bonds:
    if bond_single_frac[num_i] == 0 or frac_WPD==0:
        WPD_hbond_corr = 0
    if bond_single_frac[num_i] != 0 and frac_WPD != 0:
        WPD_hbond_corr = (WPD_hbond_frac[num_i] - (frac_WPD*bond_single_frac[num_i]))/(bond_single_frac[num_i]*frac_WPD)
    if WPD_hbond_corr >= 0.1:
        output_pos_WPD.write(name_bonds[num_i] + ' Corr = ' + str(WPD_hbond_corr) + '\n')
    if WPD_hbond_corr <= -0.1:
        output_neg_WPD.write(name_bonds[num_i] + ' Corr = ' + str(WPD_hbond_corr) + '\n')
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
