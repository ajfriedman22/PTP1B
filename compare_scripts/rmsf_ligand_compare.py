#Import Necessary Packages
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import sys

#Import custom modules
sys.path.insert(1,'/ocean/projects/cts160011p/afriedma/code/MD-Analysis/util')
import plot
import uncorr

def load_rms(path, file_name_i, combo):
    t, r = [],[]
    if combo == False:
        file_name = '../../../' + path + '/rmsf_lig_' + file_name_i + '.xvg'
    else:
        file_name = '../../../' + path + '/rmsf_' + file_name_i + '.xvg'

    with open(file_name) as f:
        for _ in range(17):
            next(f)
        for line in f:
            cols = line.split()
            if len(cols) == 2:
                t.append(float(cols[0]))
                r.append(float(cols[1])*10)
    return t, r

#File paths for all input files
file_path_AD = ['AD_rebuild_a7/analysis', 'AD/analysis', '1sug_AD/analysis', '1sug_no_a7_AD/analysis', 'AD_dis/analysis/config7', 'AD_dis/analysis/config9', 'AD_dis/analysis/config11', '1sug_dis_AD/analysis/config7', '1sug_dis_AD/analysis/config9',  '1sug_dis_AD/analysis/config11',  '1sug_dis_AD/analysis/config11_2',  '1sug_dis_AD/analysis/config_alt',  '1sug_dis_AD/analysis/config_alt2']
file_name_AD = ['a7_AD', 'AD', '1sug_AD', '1sug_na7_AD', 'AD_dis7', 'AD_dis9', 'AD_dis11', '1sug_dis_AD7', '1sug_dis_AD9', '1sug_dis_AD11', '1sug_dis_AD11_2', '1sug_dis_AD_alt', '1sug_dis_AD_alt2']

file_path_BBR = ['BBR_a7/analysis', 'BBR_1sug/analysis', 'BBR_dis/analysis/config7', 'BBR_dis/analysis/config9', 'BBR_dis/analysis/config11', 'BBR_1sug_dis/analysis/config7', 'BBR_1sug_dis/analysis/config11']
file_name_BBR = ['BBR_a7', 'BBR_1sug', 'BBR_dis7', 'BBR_dis9', 'BBR_dis11', 'BBR_dis7', 'BBR_dis11']

#open all files
RMSF_mean_AD = np.zeros(len(file_path_AD))
RMSF_mean_BBR = np.zeros(len(file_path_BBR))
RMSF_err_AD = np.zeros(len(file_path_AD))
RMSF_err_BBR = np.zeros(len(file_path_BBR))

#Indices for AD and BBR bound in close to crystal loaction
AD_ind = [9, 11, 12]
BBR_ind = [0, 3]
RMSF_AD, RMSF_BBR= [],[]
for i in range(len(file_path_AD)):
    #Load Data
    a, rmsf = load_rms(file_path_AD[i], file_name_AD[i], False)

    #Mean and SEM for each trajectory
    RMSF_mean_AD[i] = np.mean(rmsf)
    RMSF_err_AD[i] = stats.sem(rmsf)

    #If lig is bound to close to crystal save rmsd values
    if i in AD_ind:
        for j in rmsf:
            RMSF_AD.append(j)

for i in range(len(file_path_BBR)):
    #Load Data
    a, rmsf = load_rms(file_path_BBR[i], file_name_BBR[i], False)

    #Mean and SEM for each trajectory
    RMSF_mean_BBR[i] = np.mean(rmsf)
    RMSF_err_BBR[i] = stats.sem(rmsf)

    #If lig is bound to close to crystal save rmsd values
    if i in BBR_ind:
        for j in rmsf:
            RMSF_BBR.append(j)

#Name Labels
Label_AD = ['Open Ordered', 'Open Absent', 'Closed Ordered', 'Closed Absent', 'Open Dis', 'Open Dis', 'Open Dis', 'Closed Dis7', 'Closed Dis9', 'Closed Dis11', 'Closed Dis11', 'Closed Dis alt', 'Closed Dis Alt']
Label_BBR = ['Open Ordered', 'Closed Ordered', 'Open Dis7', 'Open Dis9', 'Open Dis11', 'Closed Dis7', 'Closed Dis11']

num = np.linspace(1, len(Label_AD)+1, num = len(Label_AD))
fig = plt.figure(figsize=(18,10))
ax1 = fig.add_subplot(111)
ax1.set_title("Comparison of AD RMSD")
ax1.set_ylabel(r'RMSF($\AA$)')
ax1.bar(num, RMSF_mean_AD)
plt.xticks(num, Label_AD, fontsize=14)
plt.errorbar(num, RMSF_mean_AD, yerr= RMSF_err_AD, fmt='o', color='black')
leg = ax1.legend(loc='upper right')
fig.savefig('RMSF_AD_compare.png') 
plt.close()

num = np.linspace(1, len(Label_BBR)+1, num = len(Label_BBR))
fig = plt.figure(figsize=(18,10))
ax1 = fig.add_subplot(111)
ax1.set_title("Comparison of BBR RMSF")
ax1.set_ylabel(r'RMSF($\AA$)')
ax1.bar(num, RMSF_mean_BBR)
plt.xticks(num, Label_BBR, fontsize=14)
plt.errorbar(num, RMSF_mean_BBR, yerr= RMSF_err_BBR, fmt='o', color='black')
leg = ax1.legend(loc='upper right')
fig.savefig('RMSF_BBR_compare.png') 
plt.close()

#Average AD and BBR comparison
st, p = stats.ttest_ind(RMSF_AD, RMSF_BBR, equal_var = False) #Welch's t-test b/w AD + BBR

plot.plot_gen_box(RMSF_AD, RMSF_BBR, 'AD', 'BBR', p, -1, '', r'RMSF($\AA$)', 'Comparison of RMSF b/w Ligands', 'AD_BBR_rmsf_cmpr', 6, 6)

#AD and BBR Combo simulations
#File paths for all input files
file_path_combo = ['AD_BBR/analysis']
file_name_combo = ['AD_combo', 'BBR_combo']

#Indices for AD and BBR bound in close to crystal loaction
RMSF_AD_combo, RMSF_BBR_combo= [],[]
for i in range(len(file_name_combo)):
    #Load Data
    a, rmsf = load_rms(file_path_combo[0], file_name_combo[i], True)

    #If lig is bound to close to crystal save rmsd values
    if i == 0:
        for j in rmsf:
            RMSF_AD_combo.append(j)
    if i == 1:
        for j in rmsf:
            RMSF_BBR_combo.append(j)

#Compare AD and BBR between solo simulations and combo simultaions
st, p = stats.ttest_ind(RMSF_AD, RMSF_AD_combo, equal_var = False) #Welch's t-test b/w AD + BBR
st, p2 = stats.ttest_ind(RMSF_BBR, RMSF_BBR_combo, equal_var = False) #Welch's t-test b/w AD + BBR

plot.plot_gen_box(RMSF_AD, RMSF_AD_combo, 'Solo', 'Combo', p, -1, '', r'RMSF($\AA$)', 'Comparison of RMSF b/w AD in Combo and Solo Simulation', 'AD_rmsf_combo_solo_cmpr', 6, 6)
plot.plot_gen_box(RMSF_BBR, RMSF_BBR_combo, 'Solo', 'Combo', p2, -1, '', r'RMSF($\AA$)', 'Comparison of RMSF b/w AD in Combo and Solo Simulation', 'AD_rmsf_combo_solo_cmpr', 6, 6)

