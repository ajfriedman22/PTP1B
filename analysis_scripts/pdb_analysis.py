#!/ usr / bin / env python
import math
import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
from sklearn.decomposition import PCA
from itertools import combinations
import argparse
from itertools import product
from statistics import stdev
import sys
from mpl_toolkits import mplot3d

#Import custom modules
sys.path.insert(1,'/ocean/projects/cts160011p/afriedma/code/MD-Analysis/util')
import mdfunc
import uncorr
import plot

#Declare arguments
parser = argparse.ArgumentParser(description = 'Analysis of Crystal Structures')
parser.add_argument('-p', required=True, help='File name for input PDB')
parser.add_argument('-l', required=False, default='none', help= 'Name of Ligand')

#Import Arguments
args = parser.parse_args()
File_pdb = args.p + '.pdb'
lig = args.l

if lig != 'none':
    lig_check = True
else:
    lig_check = False

#Load PDB File
pdb = md.load_pdb(File_pdb)
pdb_top = pdb.topology

#Determine if full protein is present
pdb_prot = pdb.atom_slice(pdb_top.select('protein')) #Select only atoms in the protein
if pdb_prot.n_residues == 299:
    miss_first = False
else:
    miss_first = True

#Open Output file for all analysis
output = open('PDB_analysis.txt', 'w')

#Compute WPD loop distance
if miss_first == False:
    res1 = 180 #residue 181 but from 0 start
    res2 = 214 #residue 215 but from 0 start
else:
    res1 = 179 #residue 181 but from 0 start and missing residue 1
    res2 = 213 #residue 215 but from 0 start and missing residue 1

residues = np.zeros((1,2))

residues[0][0] = pdb_top.select('resid ' + str(res1) + ' and name CA') #Select protein from the trajectory select residue 181
residues[0][1] = pdb_top.select('resid ' + str(res2) + ' and name CA') #Select protein from the trajectory select residue 215

WPD_dist = md.compute_distances(pdb, residues, periodic=False) #Compute distance between h-bond donor and acceptor

output.write('WPD Loop Distance: ' + str(WPD_dist) + ' nm\n')

if lig_check == True:
    #Compute Radius of Gyration for the Ligand
    lig_pdb = pdb.atom_slice(pdb_top.select('resname ' + lig))#Seperate ligand from PDB
    
    lig_rg = md.compute_rg(lig_pdb) #Compute radius of gyration
    
    output.write('Ligand Radius of Gyration: ' + str(lig_rg) + ' nm\n')


