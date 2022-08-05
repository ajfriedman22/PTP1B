#Import Necessary Packages
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import sys

#Import custom modules
sys.path.insert(1,'/ocean/projects/cts160011p/afriedma/code/PTP1B/util')
import plot
import uncorr

def load_rms(path, file_name_i, combo):
    t, r, raw = [],[],[]
    if combo == False:
        file_name = '../../../' + path + '/rmsd_lig_' + file_name_i + '.xvg'
    else:
        file_name = '../../../' + path + '/rmsd_' + file_name_i + '.xvg'
    with open(file_name) as f:
        for _ in range(17):
            next(f)
        for line in f:
            cols = line.split()
            if len(cols) == 2:
                t.append(float(cols[0]))
                r.append(float(cols[1])*10)
                raw.append(float(cols[1]))
    #Remove correlated data
    t_uncorr = uncorr.ind(raw)
    rmsd_uncorr = uncorr.sort(r, t_uncorr)

    return t_uncorr, rmsd_uncorr

#File paths for all input files
file_path_AD = ['AD_rebuild_a7/analysis', 'AD/analysis', '1sug_AD/analysis', '1sug_no_a7_AD/analysis', 'AD_dis/analysis/config7', 'AD_dis/analysis/config9', 'AD_dis/analysis/config11', '1sug_dis_AD/analysis/config7', '1sug_dis_AD/analysis/config9',  '1sug_dis_AD/analysis/config11',  '1sug_dis_AD/analysis/config11_2',  '1sug_dis_AD/analysis/config_alt',  '1sug_dis_AD/analysis/config_alt2']
file_name_AD = ['a7_AD', 'AD', '1sug_AD', '1sug_na7_AD', 'AD_dis7', 'AD_dis9', 'AD_dis11', '1sug_dis_AD7', '1sug_dis_AD9', '1sug_dis_AD11', '1sug_dis_AD11_2', '1sug_dis_AD_alt', '1sug_dis_AD_alt2']

file_path_BBR = ['BBR_a7/analysis', 'BBR_1sug/analysis', 'BBR_dis/analysis/config7', 'BBR_dis/analysis/config9', 'BBR_dis/analysis/config11', 'BBR_1sug_dis/analysis/config7', 'BBR_1sug_dis/analysis/config11']
file_name_BBR = ['BBR_a7', 'BBR_1sug', 'BBR_dis7', 'BBR_dis9', 'BBR_dis11', 'BBR_dis7', 'BBR_dis11']

#open all files
RMSD_mean_AD = np.zeros(len(file_path_AD))
RMSD_mean_BBR = np.zeros(len(file_path_BBR))
RMSD_err_AD = np.zeros(len(file_path_AD))
RMSD_err_BBR = np.zeros(len(file_path_BBR))

#Indices for AD and BBR bound in close to crystal loaction
AD_ind = [9, 11, 12]
BBR_ind = [0, 3]
RMSD_AD, RMSD_BBR= [],[]
for i in range(len(file_path_AD)):
    #Load Data
    a, rmsd = load_rms(file_path_AD[i], file_name_AD[i], False)

    #Mean and SEM for each trajectory
    RMSD_mean_AD[i] = np.mean(rmsd)
    RMSD_err_AD[i] = stats.sem(rmsd)

    #If lig is bound to close to crystal save rmsd values
    if i in AD_ind:
        for j in rmsd:
            RMSD_AD.append(j)

for i in range(len(file_path_BBR)):
    #Load Data
    a, rmsd = load_rms(file_path_BBR[i], file_name_BBR[i], False)

    #Mean and SEM for each trajectory
    RMSD_mean_BBR[i] = np.mean(rmsd)
    RMSD_err_BBR[i] = stats.sem(rmsd)

    #If lig is bound to close to crystal save rmsd values
    if i in BBR_ind:
        for j in rmsd:
            RMSD_BBR.append(j)

#Name Labels
Label_AD = ['Open Ordered', 'Open Absent', 'Closed Ordered', 'Closed Absent', 'Open Dis', 'Open Dis', 'Open Dis', 'Closed Dis7', 'Closed Dis9', 'Closed Dis11', 'Closed Dis11', 'Closed Dis alt', 'Closed Dis Alt']
Label_BBR = ['Open Ordered', 'Closed Ordered', 'Open Dis7', 'Open Dis9', 'Open Dis11', 'Closed Dis7', 'Closed Dis11']

num = np.linspace(1, len(Label_AD)+1, num = len(Label_AD))
fig = plt.figure(figsize=(18,10))
ax1 = fig.add_subplot(111)
ax1.set_title("Comparison of AD RMSD")
ax1.set_ylabel(r'RMSD($\AA$)')
ax1.bar(num, RMSD_mean_AD)
plt.xticks(num, Label_AD, fontsize=14)
plt.errorbar(num, RMSD_mean_AD, yerr= RMSD_err_AD, fmt='o', color='black')
leg = ax1.legend(loc='upper right')
fig.savefig('RMSD_AD_compare.png') 
plt.close()

num = np.linspace(1, len(Label_BBR)+1, num = len(Label_BBR))
fig = plt.figure(figsize=(18,10))
ax1 = fig.add_subplot(111)
ax1.set_title("Comparison of BBR RMSD")
ax1.set_ylabel(r'RMSD($\AA$)')
ax1.bar(num, RMSD_mean_BBR)
plt.xticks(num, Label_BBR, fontsize=14)
plt.errorbar(num, RMSD_mean_BBR, yerr= RMSD_err_BBR, fmt='o', color='black')
leg = ax1.legend(loc='upper right')
fig.savefig('RMSD_BBR_compare.png') 
plt.close()

#Average AD and BBR comparison
st, p = stats.ttest_ind(RMSD_AD, RMSD_BBR, equal_var = False) #Welch's t-test b/w AD + BBR

plot.plot_gen_box(RMSD_AD, RMSD_BBR, 'AD', 'BBR', p, -1, '', r'RMSD($\AA$)', 'Comparison of RMSD b/w Ligands', 'AD_BBR_rmsd_cmpr', 6, 6)

#AD and BBR Combo simulations
#File paths for all input files
file_path_combo = ['AD_BBR/analysis']
file_name_combo = ['AD_combo', 'BBR_combo']

#Indices for AD and BBR bound in close to crystal loaction
RMSD_AD_combo, RMSD_BBR_combo= [],[]
for i in range(len(file_name_combo)):
    #Load Data
    a, rmsd = load_rms(file_path_combo[0], file_name_combo[i], True)

    #If lig is bound to close to crystal save rmsd values
    if i == 0:
        for j in rmsd:
            RMSD_AD_combo.append(j)
    if i == 1:
        for j in rmsd:
            RMSD_BBR_combo.append(j)

#Average AD and BBR comparison
st, p = stats.ttest_ind(RMSD_AD_combo, RMSD_BBR_combo, equal_var = False) #Welch's t-test b/w AD + BBR

plot.plot_gen_box(RMSD_AD_combo, RMSD_BBR_combo, 'AD', 'BBR', p, -1, '', r'RMSD($\AA$)', 'Comparison of RMSD b/w Ligands in Combo Simulation', 'AD_BBR_rmsd_combo_cmpr', 6, 6)

#Compare AD and BBR between solo simulations and combo simultaions
st, p = stats.ttest_ind(RMSD_AD, RMSD_AD_combo, equal_var = False) #Welch's t-test b/w AD + BBR
st, p2 = stats.ttest_ind(RMSD_BBR, RMSD_BBR_combo, equal_var = False) #Welch's t-test b/w AD + BBR

plot.plot_gen_box(RMSD_AD, RMSD_AD_combo, 'Solo', 'Combo', p, -1, '', r'RMSD($\AA$)', 'Comparison of RMSD b/w AD in Combo and Solo Simulation', 'AD_rmsd_combo_solo_cmpr', 6, 6)
plot.plot_gen_box(RMSD_BBR, RMSD_BBR_combo, 'Solo', 'Combo', p2, -1, '', r'RMSD($\AA$)', 'Comparison of RMSD b/w AD in Combo and Solo Simulation', 'AD_rmsd_combo_solo_cmpr', 6, 6)

#COM RMSD for ligand
file_path_AD = ['1sug_dis_AD/analysis/config11', '1sug_dis_AD/analysis/config_alt', '1sug_dis_AD/analysis/config_alt2', 'mutate/WT/AD/analysis']
file_path_BBR = ['BBR_a7/analysis', 'mutate/WT/BBR/analysis']
file_path_combo = ['AD_BBR/analysis']

AD_rmsd = np.zeros(len(file_path_AD))
BBR_rmsd = np.zeros(len(file_path_BBR))

for i in range(len(file_path_AD)):
    rmsd = open('../../../' + file_path_AD[i] + '/lig_rmsd.txt').readlines()
    AD_rmsd[i] = float(rmsd[0])*10
for i in range(len(file_path_BBR)):
    rmsd = open('../../../' + file_path_BBR[i] + '/lig_rmsd.txt').readlines()
    BBR_rmsd[i] = float(rmsd[0])*10

#Plot AD RMSD vs BBR
lig_mean = [np.mean(AD_rmsd), np.mean(BBR_rmsd)]
lig_sem = [stats.sem(AD_rmsd), stats.sem(BBR_rmsd)]
lig = ['AD', 'BBR']
num = [1,2]
st, p = stats.ttest_ind(AD_rmsd, BBR_rmsd, equal_var = False) #Welch's t-test b/w AD + BBR

fig = plt.figure(figsize=(7,5))
ax1 = fig.add_subplot(111)
ax1.set_title("Comparison of Ligand COM RMSD")
ax1.set_ylabel(r'RMSD($\AA$)')
ax1.bar(num, lig_mean, color=['blue', 'purple'])
plt.xticks(num, lig, fontsize=14)
plot.error_bar(1, 2, lig_mean[0], lig_mean[1], p, 0.5, 'k')
plt.errorbar(num, lig_mean, yerr= lig_sem, fmt='o', color='black')
leg = ax1.legend(loc='best')
fig.savefig('RMSD_com_lig_compare.png') 
plt.close()

AD_combo = open('../../../AD_BBR/analysis/lig_rmsd.txt').readlines()
AD_combo = float(AD_combo[0])*10
BBR_combo = open('../../../AD_BBR/analysis/lig2_rmsd.txt').readlines()
BBR_combo = float(BBR_combo[0])*10

lig_mean = [np.mean(AD_rmsd), np.mean(BBR_rmsd)]
lig_mean_combo = [AD_combo, BBR_combo]
lig_sem = [stats.sem(AD_rmsd), stats.sem(BBR_rmsd)]
num = [1,2,3,4]
lig = ['AD', 'Ternary AD', 'BBR', 'Ternary BBR']
num1 = [1,3]
num2 = [2,4]

st, p = stats.ttest_ind(AD_rmsd, AD_combo, equal_var = False) #Welch's t-test b/w AD solo and ternary
st, p1 = stats.ttest_ind(BBR_rmsd, BBR_combo, equal_var = False) #Welch's t-test b/w BBR solo and ternary

fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.set_title("Comparison of Ligand COM RMSD", fontsize = 16)
ax1.set_ylabel(r'RMSD($\AA$)', fontsize = 14)
ax1.bar(num1, lig_mean, alpha = 0.7, color = ['blue', 'purple'])
ax1.bar(num2, lig_mean_combo, alpha = 0.7, color = ['blue', 'purple'], hatch='///')
plt.xticks(num, lig, fontsize=12)
plt.yticks(fontsize=12)
plt.errorbar(num1, lig_mean, yerr= lig_sem, fmt='o', color='black')
plot.error_bar(1, 2, lig_mean[0], lig_mean_combo[0], p, 1, 'k')
plot.error_bar(3, 4, lig_mean[1], lig_mean_combo[1], p1, 1, 'k')
leg = ax1.legend(loc='best')
fig.savefig('RMSD_com_combo_lig_compare.png') 
plt.close()

#Load bootstrap RMSD COM
AD_solo_list = open('../concat/RMSD_bootstrap_AD.txt', 'r').readlines()
BBR_solo_list = open('../concat/RMSD_bootstrap_BBR.txt', 'r').readlines()
AD_combo_list = open('../concat/RMSD_bootstrap_AD_combo.txt', 'r').readlines()
BBR_combo_list = open('../concat/RMSD_bootstrap_BBR_combo.txt', 'r').readlines()

#Convert to numpy array
AD_solo = np.zeros(len(AD_solo_list))
BBR_solo = np.zeros(len(AD_solo_list))
AD_combo = np.zeros(len(AD_solo_list))
BBR_combo = np.zeros(len(AD_solo_list))

#Convert to float and numpy array
for i in range(len(AD_solo)):
    AD_solo[i] = float(AD_solo_list[i])*10
    BBR_solo[i] = float(BBR_solo_list[i])*10
    AD_combo[i] = float(AD_combo_list[i])*10
    BBR_combo[i] = float(BBR_combo_list[i])*10

st, p = stats.ttest_ind(AD_solo, AD_combo, equal_var = False) #Welch's t-test b/w AD solo and ternary
st, p1 = stats.ttest_ind(BBR_solo, BBR_combo, equal_var = False) #Welch's t-test b/w BBR solo and ternary
st, p1 = stats.ttest_ind(AD_solo, BBR_solo, equal_var = False) #Welch's t-test b/w solo AD vs BBR

#Calculate mean and sem
lig_mean = [np.mean(AD_solo), np.mean(BBR_solo)]
lig_sem = [stats.sem(AD_solo), stats.sem(BBR_solo)]
lig_mean_combo = [np.mean(AD_combo), np.mean(BBR_combo)]
lig_sem_combo = [stats.sem(AD_combo), stats.sem(BBR_combo)]

fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.set_title("Comparison of Ligand COM RMSD", fontsize = 18)
ax1.set_ylabel(r'RMSD($\AA$)', fontsize = 16)
ax1.bar(num1, lig_mean, alpha = 0.7, color = ['blue', 'purple'])
ax1.bar(num2, lig_mean_combo, alpha = 0.7, color = ['blue', 'purple'], hatch='///')
plt.xticks(num, lig, fontsize=14)
plt.yticks(fontsize=14)
plt.errorbar(num1, lig_mean, yerr= lig_sem, fmt='o', color='black')
plt.errorbar(num2, lig_mean_combo, yerr= lig_sem_combo, fmt='o', color='black')
plot.error_bar(1, 2, lig_mean[0], lig_mean_combo[0], p, 0.5, 'k')
plot.error_bar(3, 4, lig_mean[1], lig_mean_combo[1], p1, 1, 'k')
plot.error_bar(1, 3, lig_mean[0], lig_mean[1], p2, 1, 'k')
leg = ax1.legend(loc='best')
fig.savefig('RMSD_com_boot_combo_lig_compare.png') 
plt.close()

#Print % difference
output = open('COM_RMSD_per_diff.txt', 'w')
output.write('AD vs BBR: ' + str(100*(lig_mean[0] - lig_mean[1])/lig_mean[1]) + '%\n')
output.write('AD solo vs combo: ' + str(100*(lig_mean[0] - lig_mean_combo[0])/lig_mean_combo[0]) + '%\n')
output.write('BBR solo vs combo: ' + str(100*(lig_mean[1] - lig_mean_combo[1])/lig_mean_combo[1]) + '%\n')
output.close()

