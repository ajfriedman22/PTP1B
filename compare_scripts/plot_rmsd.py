import mdtraj as md
import numpy as np
import argparse
from matplotlib import pyplot as plt

#Declare arguments
parser = argparse.ArgumentParser(description = 'Plot from Text file given Equilibrium Time')
parser.add_argument('-r', required=True, help='File name for input trajectory')
parser.add_argument('-e', required=True, type=float, help= 'Equilibrium Time')

#Import Arguments
args = parser.parse_args()
File = args.r
eq_time = args.e

#Read input file
rmsd = []
for i in open(File, 'r').readlines():
    rmsd.append(float(i)*10)

#Determine time vector
time = np.linspace(eq_time, 300, len(rmsd))

#Plot
fig = plt.figure()
plt.scatter(time, rmsd)
plt.ylabel('RMSD(A)')
plt.xlabel('Time(ns)')
fig.savefig('rmsd.png')
plt.close(fig)


