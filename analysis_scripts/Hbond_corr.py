#!/ usr / bin / env python
import mdtraj as md
import numpy as np
import argparse
import sys

#Import custom modules
sys.path.insert(1, '/ocean/projects/cts160011p/afriedma/code/PTP1B/util/')
import mdfunc

def split(word):
    return [char for char in word]

def convert(s): 
    # initialization of string to "" 
    str1 = "" 
  
    # using join function join the list s by 
    # separating words by str1 
    return(str1.join(s))
def sep_num(s):
    head = s.rstrip('0123456789')
    tail = s[len(head):]
    return tail

def deter_bond(top, res1, res2, name1, name2, i):
    bond = np.zeros(3)
    bond[0] = top.select('resid ' + str(res1[i]) + ' and name ' + str(name1[i]))
    bond[2] = top.select('resid ' + str(res2[i]) + ' and name ' + str(name2[i]))
    bond[1] = bond[0]+1
    return bond

#Declare arguments
parser = argparse.ArgumentParser(description = 'Determination of Percentage of time each individual h-bond is formed and for h-bonds Formed simultaneously')
parser.add_argument('-t', required=True, help='File name for input trajectory')
parser.add_argument('-g', required=True, help= 'File name for input topology (gro format)')
parser.add_argument('-p', required=True, help= 'File path for H-bonds of interest')

#Import Arguments
args = parser.parse_args()
File_traj = args.t + '.xtc'
File_gro = args.g + '.gro'
File_path = args.p

#Load trajectories
traj_bb, traj_prot, traj_ns, traj_a7, miss_first = mdfunc.mdtraj_load(File_traj, File_gro, [285, 297])
top = traj_ns.topology

if miss_first == True:
    offset = 2
else:
    offset = 1
print('Trajectory Loaded')


#Make array for bond names
name_bonds = open(File_path + 'Hbond_uncommon.txt', 'r').readlines()
name_note = open('/ocean/projects/cts160011p/afriedma/code/PTP1B/analysis_scripts/Hbonds_note.txt', 'r').readlines()

#Output file for bond 
output_single = open('Hbond_per_single.txt', 'w') #percentage of time each bond occurs
output_note = open('Hbond_lit_per.txt', 'w')

#Seperate residue names and atom names for each bond
res1, name1, res2, name2 = [],[],[],[]
for i in name_bonds:
    line = split(i.strip())
    #Determine indicise of dashed
    ind = []
    for j in range(len(line)):
        if line[j] == '-':
            ind.append(j)
    res1.append(str(int(sep_num(convert(line[0:ind[0]]))) - offset))
    name1.append(convert(line[ind[0]+1:ind[1]]).strip())
    res2.append(str(int(sep_num(convert(line[ind[2]+1:ind[3]]))) - offset))
    name2.append(convert(line[ind[3]+1:]).strip())

#Declare array for bond correlations
num_bonds = len(name_bonds)
bond_single_frac = np.zeros((num_bonds))
print(traj_prot.n_residues)

#track bond indicies
#Determine the percent of time each bond combination is present
for i in range(num_bonds):
    if int(res1[i]) >= 0 and int(res2[i]) >= 0 and int(res1[i]) < (traj_prot.n_residues-1) and int(res2[i]) < (traj_prot.n_residues-1):
        bond = deter_bond(top, res1, res2, name1, name2, i)
                
        #Seperate atom indices into float format
        bond_d = np.array([[float(bond[1]), float(bond[2])]])
        bond_a = np.array([[float(bond[0]), float(bond[1]), float(bond[2])]])

        #Compute distances and angles over time for both bonds
        dist = md.compute_distances(traj_ns, bond_d, periodic = False)
        angle = md.compute_angles(traj_ns, bond_a , periodic = False)
        
        #Determine the percent of time both bonds are formed
        count_s=0 #single bonds
        for k in range(len(dist)):
            if dist[k][0] <= 0.25 and angle[k][0] >= 2.094:
                count_s += 1
        output_single.write(name_bonds[i] + str(100*count_s/len(dist)) + '\n')
        bond_single_frac[i] = count_s/len(dist)
    else:
        output_single.write(name_bonds[i] + '0\n')

#Determine the percent of time each bond combination is present
#Seperate residue names and atom names for each bond
res1, name1, res2, name2 = [],[],[],[]
for i in name_note:
    line = split(i.strip())
    #Determine indicise of dashed
    ind = []
    for j in range(len(line)):
        if line[j] == '-':
            ind.append(j)
    res1.append(str(int(sep_num(convert(line[0:ind[0]]))) - offset))
    name1.append(convert(line[ind[0]+1:ind[1]]).strip())
    res2.append(str(int(sep_num(convert(line[ind[2]+1:ind[3]]))) - offset))
    name2.append(convert(line[ind[3]+1:]).strip())

num_bonds = len(name_note)
for i in range(len(name_note)):
    bond = deter_bond(top, res1, res2, name1, name2, i)

    #Seperate atom indices into float format
    bond_d = np.array([[float(bond[1]), float(bond[2])]])
    bond_a = np.array([[float(bond[0]), float(bond[1]), float(bond[2])]])
    
    #Compute distances and angles over time for both bonds
    dist = md.compute_distances(traj_ns, bond_d, periodic = False)
    angle = md.compute_angles(traj_ns, bond_a , periodic = False)
    
    #Determine the percent of time both bonds are formed
    count_s=0 #single bonds
    for k in range(len(dist)):
        if dist[k][0] <= 0.25 and angle[k][0] >= 2.094:
            count_s += 1
    output_note.write(name_bonds[i] + str(100*count_s/len(dist)) + '\n')

print('Hbond Percentages Calculated')
