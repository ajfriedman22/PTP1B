import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt
from itertools import product
import argparse

#Declare arguments
parser = argparse.ArgumentParser(description = 'Determination of ligand binding and interactions with Helices of PTP1B')
parser.add_argument('-t', required=True, help='File name for input trajectory')
parser.add_argument('-g', required=True, help= 'File name for input topology (gro format)')
parser.add_argument('-f', required=True, help='base for all file names')

#Import Arguments
args = parser.parse_args()
File_traj = args.t + '.xtc'
File_gro = args.g + '.gro'
File_base = args.f

#Load trajectories
traj = md.load(File_traj, top=File_gro)
traj_ns = traj.remove_solvent()

#Set Residue Pairs
if traj_ns.n_residues == 300 or traj_ns.n_residues == 288:
    group_l = [299]
    group_3 = [186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199, 200]
    group_4 = [221, 222, 223, 224, 225, 226, 227, 228, 229, 230, 231, 232, 234, 235, 236, 237, 238]
    group_5 = [245, 246, 247, 248, 249, 250, 251, 252]
    group_6 = [264, 265, 266, 267, 268, 269, 270, 271, 272, 273, 274, 275, 276, 277, 278, 279, 280]
    group_bend = [281, 282, 283, 284, 285, 286]
    group_7 = [287, 288, 289, 290, 291, 292, 293, 294, 295, 296, 297, 298]

if traj_ns.n_residues == 299 or traj_ns.n_residues == 287:
    group_l = [298]
    group_3 = [185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199]
    group_4 = [220, 221, 222, 223, 224, 225, 226, 227, 228, 229, 230, 231, 232, 234, 235, 236, 237]
    group_5 = [244, 245, 246, 247, 248, 249, 250, 251]
    group_6 = [263, 264, 265, 266, 267, 268, 269, 270, 271, 272, 273, 274, 275, 276, 277, 278, 279]
    group_bend = [281, 282, 283, 284, 285, 286]
    group_7 = [286, 287, 288, 289, 290, 291, 292, 293, 294, 295, 296, 297]
pair_a3 = list(product(group_l, group_3))
pair_a4 = list(product(group_l, group_4))
pair_a5 = list(product(group_l, group_5))
pair_a6 = list(product(group_l, group_6))
pair_bend = list(product(group_l, group_bend))
if traj_ns.n_residues == 300 or traj_ns.n_residues == 299:
    pair_a7 = list(product(group_l, group_7))

#compute distances
[dist_a3, pairs] = md.compute_contacts(traj_ns, contacts=pair_a3, scheme='closest', ignore_nonprotein = False, periodic=True, soft_min = False)
[dist_a4, pairs] = md.compute_contacts(traj_ns, contacts=pair_a4, scheme='closest', ignore_nonprotein = False, periodic=True, soft_min = False)
[dist_a5, pairs] = md.compute_contacts(traj_ns, contacts=pair_a5, scheme='closest', ignore_nonprotein = False, periodic=True, soft_min = False)
[dist_a6, pairs] = md.compute_contacts(traj_ns, contacts=pair_a6, scheme='closest', ignore_nonprotein = False, periodic=True, soft_min = False)
[dist_bend, pairs] = md.compute_contacts(traj_ns, contacts=pair_bend, scheme='closest', ignore_nonprotein = False, periodic=True, soft_min = False)
if traj_ns.n_residues == 300 or traj_ns.n_residues == 299:
    [dist_a7, pairs] = md.compute_contacts(traj_ns, contacts=pair_a7, scheme='closest', ignore_nonprotein = False, periodic=True, soft_min = False)

#Compute time points in which the distance to any residue is less than
time, num_pairs_a3 = np.shape(dist_a3)
time, num_pairs_a4 = np.shape(dist_a4)
time, num_pairs_a5 = np.shape(dist_a5)
time, num_pairs_a6 = np.shape(dist_a6)
time, num_pairs_bend = np.shape(dist_bend)
if traj_ns.n_residues == 300 or traj_ns.n_residues == 299:
    time, num_pairs_a7 = np.shape(dist_a7)

#Set array for the crystal structure binding location and two alternatives
contact_crys, contact_alt, contact_alt2, contact_alt3 = [],[],[],[]

#Open files for the number of interactions with each protein region at each point in time
file_a3 = open('a3_inter.txt', 'w')
file_a4 = open('a4_inter.txt', 'w')
file_a5 = open('a5_inter.txt', 'w')
file_a6 = open('a6_inter.txt', 'w')
if traj_ns.n_residues > 298:
    file_a7 = open('a7_inter.txt', 'w')

for i in range(time):
    check_a3 = 0
    check_a4 = 0
    check_a5 = 0
    check_a6 = 0
    check_bend = 0
    check_a7 = 0

    for j in range(num_pairs_a3): #Determine # of contacts with a3 helix
        if dist_a3[i][j] <= 0.5:
            check_a3 += 1
    for k in range(num_pairs_a4): #Determine # of contacts with a4 helix
        if dist_a4[i][k] <= 0.5:
            check_a4 += 1
    for l in range(num_pairs_a5):  #Determine # of contacts with a5 helix
        if dist_a5[i][l] <= 0.5:
            check_a5 += 1
    for m in range(num_pairs_a6):  #Determine # of contacts with a6 helix
        if dist_a6[i][m] <= 0.5:
            check_a6 += 1
    for o in range(num_pairs_bend): #Determine the # of contacts between a6 and a7 helices
        if dist_bend[i][o] <= 0.5:
            check_bend += 1
    
    #Output number of interactions to files for a3-a6
    file_a3.write(str(check_a3) + '\n')
    file_a4.write(str(check_a4) + '\n')
    file_a5.write(str(check_a5) + '\n')
    file_a6.write(str(check_a6) + '\n')
    
    #Save number of interactions
    total_contact_alt = check_a4 + check_a5 + check_a6 #Alternative binding location 1
    total_contact_alt2 = check_a3 + check_a4 + check_a6 #Alternative binding loaction 2
    if traj_ns.n_residues == 300 or traj_ns.n_residues == 299:
        for n in range(num_pairs_a7):#Determine # of contacts with the a7 helix pt1
            if dist_a7[i][n] <= 0.5:
                check_a7 += 1

        #Output number of interactions to file for a7
        file_a7.write(str(check_a7) + '\n')
        total_contact_crys = check_a3 + check_a6 + check_a7 #Crystal structure binding location
        total_contact_alt3 = check_a6 + check_bend + check_a7
        if check_a3 >= 1 and check_a6 >= 1 and check_a7 >= 1 and total_contact_crys >= 5 and check_a5 == 0 and check_bend == 0:
            contact_crys.append(1)
        else:
            contact_crys.append(0)
        if check_bend >= 1 and check_a6 >= 1 and check_a7 >= 1 and total_contact_alt3 >= 5 and check_a4 == 0 and check_a5 == 0:
            contact_alt3.append(1)
        else:
            contact_alt3.append(0)
    if traj_ns.n_residues == 288 or traj_ns.n_residues == 287:
        total_contact_crys = check_a3 + check_a6
        total_contact_alt3 = check_a6 + check_bend
        if check_a3 >= 1 and check_a6 >= 1 and total_contact_crys >= 5 and check_a5 == 0 and check_bend == 0 and check_a4 == 0:
            contact_crys.append(1)
        else:
            contact_crys.append(0)
        if check_bend >= 1 and check_a6 >= 1 and total_contact_alt3 >= 5 and check_a4 == 0 and check_a5 == 0:
            contact_alt3.append(1)
        else:
            contact_alt3.append(0)

    if check_a4 >= 1 and check_a6 >= 1 and total_contact_alt >= 5 and check_a3 == 0 and check_a7 == 0:
        contact_alt.append(1)
    else:
        contact_alt.append(0)
    if check_a3 >= 1 and check_a4 >= 1 and check_a6 >= 1 and total_contact_alt2 >= 5 and check_a5 == 0 and check_a7 == 0:
        contact_alt2.append(1)
    else:
        contact_alt2.append(0)
\
#Determine the fraction the ligand is in each location
AD_frac_crys = sum(contact_crys)/len(contact_crys)
AD_frac_alt = sum(contact_alt)/len(contact_alt)
AD_frac_alt2 = sum(contact_alt2)/len(contact_alt2)
AD_frac_alt3 = sum(contact_alt3)/len(contact_alt3)
AD_unbound = 1 - AD_frac_crys - AD_frac_alt - AD_frac_alt2 - AD_frac_alt3

#Print the percent time in each location
print('Crystal binding location: ' + str(100 * AD_frac_crys))
print('Alternative binding location: ' + str(100 * AD_frac_alt))
print('Alternative binding location 2: ' + str(100 * AD_frac_alt2))
print('Alternative binding location 3: ' + str(100 * AD_frac_alt3))

#Print pie chark of AD binding
labels = ['Crystal', 'Alt', 'Alt2', 'Alt3', 'Unbound']
sizes = [AD_frac_crys*100, AD_frac_alt*100, AD_frac_alt2*100, AD_frac_alt3*100, 100*AD_unbound]
fig1, ax1 = plt.subplots()
ax1.pie(sizes, labels=labels, autopct='%1.1f%%', startangle=90)
ax1.axis('equal')
fig1.savefig('AD_bind_' + File_base + '_pie.png')
plt.close(fig1)

