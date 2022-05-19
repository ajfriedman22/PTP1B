#!/ usr / bin / env python

#Import packages
from matplotlib import pyplot as plt
import numpy as np
import sys
from scipy import stats

#Import custom modules
sys.path.insert(1,'/ocean/projects/cts160011p/afriedma/code/PTP1B/util')
import mdfunc

#Read in data from input files
WT_AD = open('../../WT/AD/analysis/DSSP_WT.txt', 'r').readlines()
WT_BBR = open('../../WT/BBR/analysis/DSSP_WT_BBR.txt', 'r').readlines()
F196A_AD = open('../../F196A/AD/analysis/DSSP_F196A.txt', 'r').readlines()
F196A_BBR = open('../../F196A/BBR/analysis/DSSP_F196A_BBR.txt', 'r').readlines()
L192A_AD = open('../../L192A/AD/analysis/DSSP_L192A.txt', 'r').readlines()
L192A_BBR = open('../../L192A/BBR/analysis/DSSP_L192A_BBR.txt', 'r').readlines()
L192F_AD = open('../../L192F/AD/analysis/DSSP_L192F.txt', 'r').readlines()
L192F_BBR = open('../../L192F/BBR/analysis/DSSP_L192F_BBR.txt', 'r').readlines()
L195A_AD = open('../../L195A/AD/analysis/DSSP_L195A.txt', 'r').readlines()
L195A_BBR = open('../../L195A/BBR/analysis/DSSP_L195A_BBR.txt', 'r').readlines()
L195F_AD = open('../../L195F/AD/analysis/DSSP_L195F.txt', 'r').readlines()
L195F_BBR = open('../../L195F/BBR/analysis/DSSP_L195F_BBR.txt', 'r').readlines()
L195N_AD = open('../../L195N/AD/analysis/DSSP_L195N.txt', 'r').readlines()
L195N_BBR = open('../../L195N/BBR/analysis/DSSP_L195N_BBR.txt', 'r').readlines()
S286A_AD = open('../../S286A/AD/analysis/DSSP_S286A_AD.txt', 'r').readlines()
S286A_BBR = open('../../S286A/BBR/analysis/DSSP_S286A_BBR.txt', 'r').readlines()
F280Y_AD = open('../../F280Y/AD/analysis/DSSP_F280Y.txt', 'r').readlines()
F280Y_BBR = open('../../F280Y/BBR/analysis/DSSP_F280Y_BBR.txt', 'r').readlines()
E276L_AD = open('../../E276L/AD/analysis/DSSP_E276L_AD.txt', 'r').readlines()
E276L_BBR = open('../../E276L/BBR/analysis/DSSP_E276L_BBR.txt', 'r').readlines()
E276F_AD = open('../../E276F/AD/analysis/DSSP_E276F_AD.txt', 'r').readlines()
E276F_BBR = open('../../E276F/BBR/analysis/DSSP_E276F_BBR.txt', 'r').readlines()
K279M_AD = open('../../K279M/AD/analysis/DSSP_K279M_AD.txt', 'r').readlines()
K279M_BBR = open('../../K279M/BBR/analysis/DSSP_K279M_BBR.txt', 'r').readlines()
K279W_AD = open('../../K279W/AD/analysis/DSSP_K279W_AD.txt', 'r').readlines()
K279W_BBR = open('../../K279W/BBR/analysis/DSSP_K279W_BBR.txt', 'r').readlines()
V287T_AD = open('../../V287T/AD/analysis/DSSP_V287T_AD.txt', 'r').readlines()
V287T_BBR = open('../../V287T/BBR/analysis/DSSP_V287T_BBR.txt', 'r').readlines()

#Declare array of zeros for num of times the residues are alpha helical and helix in general
per_struct = np.zeros([28, 8])
per_alpha_helix = np.zeros([28, 8])

#Seperate Characters in the string and record number that are in an alpha helix
data = [WT_AD, WT_BBR, F196A_AD, F196A_BBR, L192A_AD, L192A_BBR, L192F_AD, L192F_BBR, L195A_AD, L195A_BBR, L195F_AD, L195F_BBR, L195N_AD, L195N_BBR, S286A_AD, S286A_BBR, F280Y_AD, F280Y_BBR, E276L_AD, E276L_BBR, E276F_AD, E276F_BBR, K279M_AD, K279M_BBR, K279W_AD, K279W_BBR, V287T_AD, V287T_BBR]

for i in range(28):
    per_alpha_helix[i][:], per_struct[i][:], alpha_per_mean, struct_per_mean, alpha_per_sem, struct_per_sem = mdfunc.per_helix(data[i][:], False)

#Compare Ligand Bound Disordered Structures
num = [1, 2, 3, 4, 5, 6, 7, 8]
Method = ['Res 287', 'Res 288', 'Res 289', 'Res 290', 'Res 291', 'Res 292', 'Res 293', 'Res 294']
fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.set_title('Degree of Structure in AD Bound Disordered Conformations')
ax1.set_ylabel('% Residue was in alpha helix')
ax1.set_ylim(0,100)
ax1.plot(num,per_alpha_helix[0][:], color = 'red', label='WT AD')
plt.xticks(num, Method)
leg = ax1.legend()
fig.savefig('Helicity_cmpr_AD.png')
plt.close(fig)

#Caluclate mean and SD for percent helicity in the alpha 7 helix
mut = ['WT', 'F196A', 'L192A', 'L192F', 'L195A', 'L195F', 'L195N', 'S286A', 'F280Y', 'E276L', 'E276F', 'K279M', 'K279W', 'V287T']
mut_spef = ['WT', 'F196A', 'L192F', 'L195F', 'F280Y', 'E276F', 'V287T']
num = np.linspace(0, len(mut), num = len(mut))
num_spef = np.linspace(0, len(mut_spef), num = len(mut_spef))
mean_alpha_hel, err_alpha_hel = [],[]

for i in range(28):
    mean_alpha_hel.append(np.mean(per_alpha_helix[i][:]))
    err_alpha_hel.append(stats.sem(per_alpha_helix[i][:]))

#Set Colors
colors = ['blue']
alpha_AD = np.zeros(len(mut))
alpha_err_AD = np.zeros(len(mut))
i = 0
for n in [0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26]:
    #Calculated p-values
    st, p = stats.ttest_ind(per_alpha_helix[n][:], per_alpha_helix[0][:], equal_var = False) #Welch's t-test
    diff = mean_alpha_hel[n] - mean_alpha_hel[0]
    if diff > 0 and p < 0.05:
        colors.append('orange')
    elif diff > 0 and p < 0.05:
        colors.append('green')
    else:
        colors.append('gray')
    alpha_AD[i] = mean_alpha_hel[n]
    alpha_err_AD[i] = err_alpha_hel[n]
    i += 1

fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.set_title('Comparison of Helicity in the alpha-7 Helix')    
ax1.set_ylabel('Mean Precent Helicity')
ax1.bar(num, alpha_AD, color=colors, width=0.9)
plt.errorbar(num, alpha_AD, yerr= alpha_err_AD, fmt='o', color='black')
plt.xticks(num, mut, fontsize=8)
fig.savefig('Mean_helicity_AD_all.png')
plt.close(fig)

#Set Colors
colors = ['blue']
alpha_AD = np.zeros(len(mut_spef))
alpha_err_AD = np.zeros(len(mut_spef))
i = 0
for n in [0, 2, 6, 10, 16, 20, 26]:
    #Calculated p-values
    st, p = stats.ttest_ind(per_alpha_helix[n][:], per_alpha_helix[0][:], equal_var = False) #Welch's t-test
    diff = mean_alpha_hel[n] - mean_alpha_hel[0]
    if diff > 0 and p < 0.05:
        colors.append('orange')
    elif diff > 0 and p < 0.05:
        colors.append('green')
    else:
        colors.append('gray')
    alpha_AD[i] = mean_alpha_hel[n]
    alpha_err_AD[i] = err_alpha_hel[n]
    i += 1

fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.set_title('Comparison of Helicity in the alpha-7 Helix')    
ax1.set_ylabel('Mean Precent Helicity')
ax1.bar(num_spef, alpha_AD, color=colors, width=0.9)
plt.errorbar(num_spef, alpha_AD, yerr= alpha_err_AD, fmt='o', color='black')
plt.xticks(num_spef, mut_spef, fontsize=8)
fig.savefig('Mean_helicity_AD.png')
plt.close(fig)

#Set Colors
colors = ['blue']
alpha_BBR = np.zeros(len(mut))
alpha_err_BBR = np.zeros(len(mut))
i = 0
for n in [1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27]:
    #Calculated p-value
    st, p = stats.ttest_ind(per_alpha_helix[n][:], per_alpha_helix[1][:], equal_var = False) #Welch's t-test
    diff = mean_alpha_hel[n] - mean_alpha_hel[1]
    if diff > 0 and p < 0.05:
        colors.append('orange')
    elif diff > 0 and p < 0.05:
        colors.append('green')
    else:
        colors.append('gray')
    alpha_BBR[i] = mean_alpha_hel[n]
    alpha_err_BBR[i] = err_alpha_hel[n]
    i += 1

fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.set_title('Comparison of Helicity in the alpha-7 Helix')    
ax1.set_ylabel('Mean Precent Helicity')
ax1.bar(num, alpha_BBR, color=colors, width=0.9)
plt.errorbar(num, alpha_BBR, yerr= alpha_err_BBR, fmt='o', color='black')
plt.xticks(num, mut, fontsize=8)
fig.savefig('Mean_helicity_BBR_all.png')
plt.close(fig)

#Set Colors
colors = ['blue']
alpha_BBR = np.zeros(len(mut_spef))
alpha_err_BBR = np.zeros(len(mut_spef))
i = 0
for n in [1, 3, 7, 11, 17, 21, 27]:
    #Calculated p-values
    st, p = stats.ttest_ind(per_alpha_helix[n][:], per_alpha_helix[0][:], equal_var = False) #Welch's t-test
    diff = mean_alpha_hel[n] - mean_alpha_hel[0]
    if diff > 0 and p < 0.05:
        colors.append('orange')
    elif diff > 0 and p < 0.05:
        colors.append('green')
    else:
        colors.append('gray')
    alpha_BBR[i] = mean_alpha_hel[n]
    alpha_err_BBR[i] = err_alpha_hel[n]
    i += 1

fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.set_title('Comparison of Helicity in the alpha-7 Helix')    
ax1.set_ylabel('Mean Precent Helicity')
ax1.bar(num_spef, alpha_BBR, color=colors, width=0.9)
plt.errorbar(num_spef, alpha_BBR, yerr= alpha_err_BBR, fmt='o', color='black')
plt.xticks(num_spef, mut_spef, fontsize=8)
fig.savefig('Mean_helicity_BBR.png')
plt.close(fig)

