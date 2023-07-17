#Import Necessary Packages
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import sys

#Import custom modules
sys.path.insert(1,'/ocean/projects/cts160011p/afriedma/code/MD-Analysis/util')
import plot
import uncorr

def load_data(path):
    dist_193 = open('../../../' + path + '/dist_lig_193.txt', 'r').readlines()
    dist_276 = open('../../../' + path + '/dist_lig_276.txt', 'r').readlines()
    
    dist_193 = dist_193[:-1]
    dist_276 = dist_276[:-1]

    return dist_193, dist_276

#File paths for all input files
file_path_AD = ['1sug_dis_AD/analysis/config11',  '1sug_dis_AD/analysis/config_alt',  '1sug_dis_AD/analysis/config_alt2']

file_path_BBR = ['BBR_a7/analysis']

#open all files
dist_193_AD, dist_276_AD, dist_193_BBR, dist_276_BBR = [],[],[],[]

for i in range(len(file_path_AD)):
    dist_193, dist_276 = load_data(file_path_AD[i])
    for x in dist_193:
        dist_193_AD.append(float(x)*10) #convert to A from nm
    for y in dist_276:
        dist_276_AD.append(float(y)*10)

for i in range(len(file_path_BBR)):
    dist_193, dist_276 = load_data(file_path_BBR[i])
    for x in dist_193:
        dist_193_BBR.append(float(x)*10)
    for y in dist_276:
        dist_276_BBR.append(float(y)*10)

#Measure displacment of ligand relative to two references
displace_AD = np.zeros((len(dist_193_AD), 2))
displace_BBR = np.zeros((len(dist_193_BBR), 2))

min_193_AD = min(dist_193_AD)
min_276_AD = min(dist_276_AD)
min_193_BBR = min(dist_193_BBR)
min_276_BBR = min(dist_276_BBR)

for i in range(len(dist_193_AD)):
    displace_AD[i][0] = dist_193_AD[i] - min_193_AD
    displace_AD[i][1] = dist_276_AD[i] - min_276_AD
for i in range(len(dist_193_BBR)):
    displace_BBR[i][0] = dist_193_BBR[i] - min_193_BBR
    displace_BBR[i][1] = dist_276_BBR[i] - min_276_BBR

#Caluclate P value
st, p = stats.ttest_ind(displace_AD[:,0], displace_BBR[:,0], equal_var = False) #Welch's t-test b/w AD + BBR
st, p2 = stats.ttest_ind(displace_AD[:,1], displace_BBR[:,1], equal_var = False) #Welch's t-test b/w AD + BBR

plot.plot_gen_box(displace_AD[:,0], displace_BBR[:,0], 'AD', 'BBR', p, 10, ' ', r'Displacment ($\AA$)', r'Displacment Relative to the $\alpha$3 Helix', 'dist_lig_a3', 6, 6)
plot.plot_gen_box(displace_AD[:,1], displace_BBR[:,1], 'AD', 'BBR', p2, 10, ' ', r'Displacment ($\AA$)', r'Displacment Relative to the $\alpha$6 Helix', 'dist_lig_a6', 6, 6)


