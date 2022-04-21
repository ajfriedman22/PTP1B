#!/ usr / bin / env python

#Import packages
from matplotlib import pyplot as plt
import numpy as np
from scipy import stats

#Function for determining helix orientation
def determine_helix(list_dssp):
    num_struct = np.zeros(10)
    num_alpha_helix = np.zeros(10)
    per_struct = np.zeros(8)
    per = np.zeros(8)
    per_sect = np.zeros(4)
    list_len = 0
    for char in list_dssp:
        if char[2]=='H':
            num_alpha_helix[0] +=1
        if char[2] == 'T' or char[0] == 'G' or char[0] == 'I':
            num_struct[0] += 1
        if char[4]=='H':
            num_alpha_helix[1] +=1
        if char[4] == 'T' or char[2] == 'G' or char[2] == 'I':
            num_struct[1] += 1
        if char[6]=='H':
            num_alpha_helix[2] += 1
        if char[6] == 'T' or char[4] == 'G' or char[4] == 'I':
            num_struct[2] += 1
        if char[8]=='H':
            num_alpha_helix[3] +=1
        if char[8] == 'T' or char[6] == 'G' or char[6] == 'I':
            num_struct[3] += 1
        if char[10]=='H':
            num_alpha_helix[4] += 1
        if char[10] == 'T' or char[8] == 'G' or char[8] == 'I':
            num_struct[4] += 1
        if char[12]=='H':
            num_alpha_helix[5] += 1
        if char[12] == 'T' or char[10] == 'G' or char[10] == 'I':
            num_struct[5] += 1
        if char[14]=='H':
            num_alpha_helix[6] += 1
        if char[14] == 'T' or char[12] == 'G' or char[12] == 'I':
            num_struct[6] += 1
        if char[16]=='H':
            num_alpha_helix[7] +=1
        if char[16] == 'T' or char[14] == 'G' or char[14] == 'I':
            num_struct[7] += 1
        if char[18]=='H':
            num_alpha_helix[8] +=1
        if char[18] == 'T' or char[16] == 'G' or char[16] == 'I':
            num_struct[8] += 1
        if char[20]=='H':
            num_alpha_helix[9] +=1
        if char[20] == 'T' or char[18] == 'G' or char[18] == 'I':
            num_struct[9] += 1
        list_len += 1
    print(list_len)

    per[0] = round(100*num_alpha_helix[2]/list_len)
    per[1] = round(100*num_alpha_helix[3]/list_len)
    per[2] = round(100*num_alpha_helix[4]/list_len) 
    per[3] = round(100*num_alpha_helix[5]/list_len) 
    per[4] = round(100*num_alpha_helix[6]/list_len) 
    per[5] = round(100*num_alpha_helix[7]/list_len) 
    per[6] = round(100*num_alpha_helix[8]/list_len) 
    per[7] = round(100*num_alpha_helix[9]/list_len)

    per_sect[0] = per[3]
    per_sect[1] = per[4]
    per_sect[2] = per[5]
    per_sect[3] = per[6]

    per_struct[0] = round(100*num_struct[2]/list_len)
    per_struct[1] = round(100*num_struct[3]/list_len)
    per_struct[2] = round(100*num_struct[4]/list_len) 
    per_struct[3] = round(100*num_struct[5]/list_len) 
    per_struct[4] = round(100*num_struct[6]/list_len) 
    per_struct[5] = round(100*num_struct[7]/list_len) 
    per_struct[6] = round(100*num_struct[8]/list_len) 
    per_struct[7] = round(100*num_struct[9]/list_len)
    return num_alpha_helix, num_struct, per_struct, per, per_sect

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
num_alpha_helix[0][:], num_struct[0][:], per_struct[0][:], per_alpha_helix[0][:], per_alpha_sect[0][:] = determine_helix(WT)
num_alpha_helix[1][:], num_struct[1][:], per_struct[1][:], per_alpha_helix[1][:], per_alpha_sect[1][:]= determine_helix(WT_close)
num_alpha_helix[2][:], num_struct[2][:], per_struct[2][:], per_alpha_helix[2][:], per_alpha_sect[2][:] = determine_helix(F196A)
num_alpha_helix[3][:], num_struct[3][:], per_struct[3][:], per_alpha_helix[3][:], per_alpha_sect[3][:] = determine_helix(L192F)
num_alpha_helix[4][:], num_struct[4][:], per_struct[4][:], per_alpha_helix[4][:], per_alpha_sect[4][:] = determine_helix(L195F)
num_alpha_helix[5][:], num_struct[5][:], per_struct[5][:], per_alpha_helix[5][:], per_alpha_sect[5][:] = determine_helix(F280Y)
num_alpha_helix[6][:], num_struct[6][:], per_struct[6][:], per_alpha_helix[6][:], per_alpha_sect[6][:] = determine_helix(E276F)
num_alpha_helix[7][:], num_struct[7][:], per_struct[7][:], per_alpha_helix[7][:], per_alpha_sect[7][:] = determine_helix(V287T)

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
mean_mut_alpha_sect, err_mut_alpha_sect = [],[]
for i in range(len(num)):
    mean_alpha_hel.append(np.mean(per_alpha_helix[i][:]))
    err_alpha_hel.append(stats.sem(per_alpha_helix[i][:]))
    mean_alpha_sect.append(np.mean(per_alpha_sect[i][:]))
    err_alpha_sect.append(stats.sem(per_alpha_sect[i][:]))
    if i > 1:
        mean_mut_alpha_sect.append(np.mean(per_alpha_sect[i][:]))
        err_mut_alpha_sect.append(stats.sem(per_alpha_sect[i][:]))

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
F_BBR = [-0.259789014, -0.303991355, -0.488133372, -0.447825789, -0.243071805, 0.164257023]
F_BBR_se = [0.001039749, 0.001238509, 0.001227868, 0.001397659, 0.001188003, 0.000777901]
F_AD = [-0.204972546, -0.035556769, -0.2663284, -0.286222165, -0.331721332, 0.131296446]
F_AD_se = [0.000911951, 0.001351002, 0.000930425, 0.001033858, 0.001059996, 0.000928346]

#Gaussian Fit
mean, std = norm.fit(F_AD)
mean2, std2 = norm.fit(F_BBR)
x = np.linspace(0, 100, num=1000)

#Make plot
fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.set_title(r'Comparison of Helicity in the $\alpha$-7 Helix Middle Section')    
ax1.set_xlabel('Mean Precent Helicity')
ax1.set_ylabel('Relative Inhibition (F)')
ax1.scatter(mean_mut_alpha_sect, F_AD, label = 'AD', color = 'blue')
plt.errorbar(mean_mut_alpha_sect, F_AD, xerr= err_mut_alpha_sect, yerr = F_AD, fmt='o', color='blue')
plt.plot(x, norm.pdf(x, mean, std))
ax1.scatter(mean_mut_alpha_sect, F_BBR, label = 'BBR', color = 'purple')
plt.errorbar(mean_mut_alpha_sect, F_BBR, xerr= err_mut_alpha_sect, yerr = F_BBR, fmt='o', color='purple')
plt.plot(x, norm.pdf(x, mean2, std2))
fig.savefig('helicity_v_F.png')
plt.close(fig)


