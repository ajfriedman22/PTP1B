#!/ usr / bin / env python

#Import packages
from matplotlib import pyplot as plt
import numpy as np
from scipy import stats
from scipy.stats import norm
from scipy.optimize import curve_fit
import sys

#Import custom modules
sys.path.insert(1,'/ocean/projects/cts160011p/afriedma/code/MD-Analysis/util')
import mdfunc

#Read in data from input files
WT = open('../../../Apo_dis/analysis/DSSP_alt_Apo.txt', 'r').readlines()
WT_close = open('../../../Apo_1SUG/analysis/1sug/DSSP_1sug.txt', 'r').readlines()
F196A = open('../../F196A/Apo/analysis/DSSP_F196A_Apo.txt', 'r').readlines()
L192F = open('../../L192F/Apo/analysis/DSSP_L192F_Apo.txt', 'r').readlines()
L195F = open('../../L195F/Apo/analysis/DSSP_L195F_Apo.txt', 'r').readlines()
F280Y = open('../../F280Y/Apo/analysis/DSSP_F280Y_Apo.txt', 'r').readlines()
E276F = open('../../E276F/Apo/analysis/DSSP_E276F_Apo.txt', 'r').readlines()
V287T = open('../../V287T/Apo/analysis/DSSP_V287T_Apo.txt', 'r').readlines()

#Declare array of zeros for num of times the residues are alpha helical and helix in general
num_struct = np.zeros([8, 10])
num_alpha_helix = np.zeros([8, 10])
per_struct = np.zeros([8, 8])
per_alpha_helix = np.zeros([8, 8])
per_alpha_sect = np.zeros([8, 4])

#Seperate Characters in the string and record number that are in an alpha helix
num_alpha_helix[0][:], num_struct[0][:], per_struct[0][:], per_alpha_helix[0][:], per_alpha_sect[0][:] = mdfunc.determine_helix(WT)
num_alpha_helix[1][:], num_struct[1][:], per_struct[1][:], per_alpha_helix[1][:], per_alpha_sect[1][:]= mdfunc.determine_helix(WT_close)
num_alpha_helix[2][:], num_struct[2][:], per_struct[2][:], per_alpha_helix[2][:], per_alpha_sect[2][:] = mdfunc.determine_helix(F196A)
num_alpha_helix[3][:], num_struct[3][:], per_struct[3][:], per_alpha_helix[3][:], per_alpha_sect[3][:] = mdfunc.determine_helix(L192F)
num_alpha_helix[4][:], num_struct[4][:], per_struct[4][:], per_alpha_helix[4][:], per_alpha_sect[4][:] = mdfunc.determine_helix(L195F)
num_alpha_helix[5][:], num_struct[5][:], per_struct[5][:], per_alpha_helix[5][:], per_alpha_sect[5][:] = mdfunc.determine_helix(F280Y)
num_alpha_helix[6][:], num_struct[6][:], per_struct[6][:], per_alpha_helix[6][:], per_alpha_sect[6][:] = mdfunc.determine_helix(E276F)
num_alpha_helix[7][:], num_struct[7][:], per_struct[7][:], per_alpha_helix[7][:], per_alpha_sect[7][:] = mdfunc.determine_helix(V287T)

#Compare Ligand Bound Disordered Structures
num = [1, 2, 3, 4, 5, 6, 7, 8]
Method = ['Res 287', 'Res 288', 'Res 289', 'Res 290', 'Res 291', 'Res 292', 'Res 293', 'Res 294']
fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.set_title('Degree of Structure in AD Bound Disordered Conformations')
ax1.set_ylabel('% Residue was in alpha helix')
ax1.set_ylim(0,100)
ax1.plot(num,per_alpha_helix[0][:], label='WT Open')
ax1.plot(num,per_alpha_helix[1][:], label='WT Close')
ax1.plot(num,per_alpha_helix[2][:], label='F196A')
ax1.plot(num,per_alpha_helix[3][:], label='L192F')
ax1.plot(num,per_alpha_helix[4][:], label='L195F')
ax1.plot(num,per_alpha_helix[5][:], label='F280Y')
ax1.plot(num,per_alpha_helix[6][:], label='E276F')
ax1.plot(num,per_alpha_helix[7][:], label='V287T')
plt.xticks(num, Method)
leg = ax1.legend()
fig.savefig('Helicity_cmpr_Apo.png')
plt.close(fig)

#Caluclate mean and SD for percent helicity in the alpha 7 helix
mut = ['WT Open', 'WT Close', 'F196A', 'L192F', 'L195F', 'F280Y', 'E276F', 'V287T']
num = [0, 1, 2, 3, 4, 5, 6, 7]
mean_alpha_hel, err_alpha_hel = [],[]
mean_alpha_sect, err_alpha_sect = [],[]
mean_mut_alpha, err_mut_alpha = [],[]
mean_struct, err_struct = [],[]
mean_mut_struct, err_mut_struct = [],[]
mean_mut_struct_lim, err_mut_struct_lim = [],[]
output = open('DSSP_alpha_apo.txt', 'w')#Output file
output2 = open('DSSP_struct_apo.txt', 'w')#Output file
for i in range(len(num)):
    mean_alpha_hel.append(np.mean(per_alpha_helix[i][:]))
    err_alpha_hel.append(stats.sem(per_alpha_helix[i][:]))
    mean_alpha_sect.append(np.mean(per_alpha_sect[i][:]))
    err_alpha_sect.append(stats.sem(per_alpha_sect[i][:]))
    mean_struct.append(np.mean(per_struct[i][:]))
    err_struct.append(stats.sem(per_struct[i][:]))
    output.write(mut[i] + ': ' + str(mean_alpha_hel[i]) + ' +/- ' + str(err_alpha_hel[i]) + '\n')
    output2.write(mut[i] + ': ' + str(mean_struct[i]) + ' +/- ' + str(err_struct[i]) + '\n')
    if i != 1:
        mean_mut_alpha.append(np.mean(per_alpha_helix[i][:]))
        err_mut_alpha.append(stats.sem(per_alpha_helix[i][:]))
        mean_mut_struct.append(np.mean(per_struct[i][:]))
        err_mut_struct.append(stats.sem(per_struct[i][:]))
        if i != 6:
            mean_mut_struct_lim.append(np.mean(per_struct[i][:]))
            err_mut_struct_lim.append(stats.sem(per_struct[i][:]))


#Set Colors
colors = ['blue', 'red']
for i in range(6):
    #Calculated p-values
    n = i + 2
    st, p = stats.ttest_ind(per_alpha_helix[n][:], per_alpha_helix[0][:], equal_var = False) #Welch's t-test
    diff = mean_alpha_hel[n] - mean_alpha_hel[0]
    if diff > 0 and p < 0.05:
        colors.append('orange')
    elif diff > 0 and p < 0.05:
        colors.append('green')
    else:
        colors.append('gray')

fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.set_title('Comparison of Helicity in the alpha-7 Helix')    
ax1.set_ylabel('Mean Precent Helicity')
ax1.bar(num, mean_alpha_hel, color=colors, width=0.9)
plt.errorbar(num, mean_alpha_hel, yerr= err_alpha_hel, fmt='o', color='black')
plt.xticks(num, mut, fontsize=8)
fig.savefig('Mean_helicity_Apo.png')
plt.close(fig)

#Set Colors
colors = ['blue', 'red']
for i in range(6):
    #Calculated p-values
    n = i + 2
    st, p = stats.ttest_ind(per_alpha_sect[n][:], per_alpha_sect[0][:], equal_var = False) #Welch's t-test
    diff = mean_alpha_sect[n] - mean_alpha_sect[0]
    if p < 0.05:
        colors.append('orange')
    else:
        colors.append('gray')

fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.set_title('Comparison of Helicity in the alpha-7 Helix Middle Section')    
ax1.set_ylabel('Mean Precent Helicity')
ax1.bar(num, mean_alpha_sect, color=colors, width=0.9)
plt.errorbar(num, mean_alpha_sect, yerr= err_alpha_sect, fmt='o', color='black')
plt.xticks(num, mut, fontsize=8)
fig.savefig('Mean_helicity_sect_Apo.png')
plt.close(fig)

#Plot DSSP against relative inhibition
F_BBR = [0.00, -0.259789014, -0.303991355, -0.488133372, 0.164257, -0.447825789, -0.243071805]
F_BBR_se = [0.00, 0.001039749, 0.001238509, 0.001227868, 0.000777901, 0.001397659, 0.001188003]
F_AD = [0.00, -0.204972546, -0.035556769, -0.2663284, 0.131296446, -0.286222165, -0.331721332]
F_AD_se = [0.00, 0.000911951, 0.001351002, 0.000930425, 0.000928346, 0.001033858, 0.001059996]

#Make plot
fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.set_title(r'Comparison of Helicity in the $\alpha$-7 Helix')   
ax1.set_xlabel(r'Mean Precent $\alpha$ Helicity')
ax1.set_ylabel('Relative Inhibition (F)')
ax1.scatter(mean_mut_alpha, F_AD, label = 'AD', color = 'blue')
plt.errorbar(mean_mut_alpha, F_AD, xerr= err_mut_alpha, yerr = F_AD_se, fmt='o', color='blue')
fig.savefig('helicity_v_F.png')
plt.close(fig)

# Define the Gaussian function
def Gauss(x, a, x0, sigma, b):
    return a*np.exp(-(x-x0)**2/(2*sigma**2)) + b
parameters, covariance = curve_fit(Gauss, mean_mut_struct, F_AD, p0 = [0.2, 35, 10, -0.25])
  
[fit_A, fit_B, fit_C, fit_D] = parameters

x = np.linspace(0, 75, num=1000)
Label = ['WT', 'F196A', 'L192F', 'L195F', 'F280Y', 'E276F', 'V287T']

fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.set_title(r'Comparison of Helicity in the $\alpha$-7 Helix')    
ax1.set_xlabel('Mean Precent Helicity')
ax1.set_ylabel('Relative Inhibition (F)')
ax1.scatter(mean_mut_struct, F_AD, label = 'AD', color = 'blue')
plt.errorbar(mean_mut_struct, F_AD, xerr= err_mut_struct, yerr = F_AD_se, fmt='o', color='blue')
# Loop for annotation of all points
for i in range(len(Label)):
    plt.annotate(Label[i], (mean_mut_struct[i], F_AD[i] + 0.05))
ax1.fill_between([0, 70], -0.1, 0.1, color = 'gray', alpha=0.5)
ax1.fill_between([0, 70], -0.5, -0.1, color = 'red', alpha=0.5)
ax1.fill_between([0, 70], 0.1, 0.3, color = 'green', alpha=0.5)
#plt.plot(x, Gauss(x, fit_A, fit_B, fit_C, fit_D), label = 'Gaussian Fit')
fig.savefig('helicity_v_F_struct.png')
plt.close(fig)

#Plot helicity vs dG for AD
dG = [0.00, 0.437, -0.375, -0.044, 0.634, -1.526, 0.056]
dG_MBAR = [0.00, 0.381, -0.297, -0.019, 0.236, -0.373, -0.026]
dG_MBAR_err = [0.00, 0.109, 0.095, 0.074, 0.072, 0.694, 0.106]

fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.set_title(r'Comparison of Helicity in the $\alpha$-7 Helix')
ax1.set_xlabel(r'$\delta\delta$G from WT')
ax1.set_ylabel(r'Mean Percent $\alpha$ Helicity')
ax1.scatter(mean_mut_alpha, dG_MBAR, label = 'AD', color = 'blue')
plt.errorbar(mean_mut_alpha, dG_MBAR, xerr= err_mut_alpha, yerr = dG_MBAR_err, fmt='o', color='blue')
fig.savefig('helicity_v_dG.png')
plt.close(fig)

#Plot ddG vs F
fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.set_title(r'Free Energy Difference of Mutation vs Relative Inhibition') 
ax1.set_ylabel(r'$\Delta\Delta$G from WT(kcal/mol)')
ax1.set_xlabel('Relative Inhibition(F)')
ax1.scatter(F_AD, dG_MBAR, label = 'AD', color = 'blue')
plt.errorbar(F_AD, dG_MBAR, xerr= F_AD_se, yerr = dG_MBAR_err, fmt='o', color='blue')
ax1.fill_between([-0.4, 0.15], -0.1, 0.1, color = 'gray', alpha=0.5)
ax1.fill_between([-0.4, 0.15], -1.0, -0.1, color = 'green', alpha=0.5)
ax1.fill_between([-0.4, 0.15], 0.1, 0.5, color = 'red', alpha=0.5)
fig.savefig('F_v_dG.png')
plt.close(fig)

dG = [0.00, 0.437, -0.375, -0.044, 0.634, 0.056]
dG_MBAR = [0.00, 0.381, -0.297, -0.019, 0.236, -0.026]
dG_MBAR_err = [0.00, 0.109, 0.095, 0.074, 0.072, 0.106]
Label_lim = ['WT', 'F196A', 'L192F', 'L195F', 'F280Y', 'V287T']

fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.set_title(r'Comparison of Helicity in the $\alpha$-7 Helix')    
ax1.set_ylabel(r'$\Delta\Delta$G from WT (kcal/mol)')
ax1.set_xlabel('Mean Percent Helicity')
ax1.scatter(mean_mut_struct_lim, dG_MBAR, label = 'AD', color = 'blue')
plt.errorbar(mean_mut_struct_lim, dG_MBAR, xerr= err_mut_struct_lim, yerr = dG_MBAR_err, fmt='o', color='blue')
# Loop for annotation of all points
for i in range(len(Label_lim)):
    plt.annotate(Label_lim[i], (mean_mut_struct_lim[i], dG_MBAR[i] + 0.05))
ax1.fill_between([0, 70], -0.1, 0.1, color = 'gray', alpha=0.5)
ax1.fill_between([0, 70], -0.5, -0.1, color = 'green', alpha=0.5)
ax1.fill_between([0, 70], 0.1, 0.5, color = 'red', alpha=0.5)
fig.savefig('helicity_v_dG_struct.png')
plt.close(fig)

