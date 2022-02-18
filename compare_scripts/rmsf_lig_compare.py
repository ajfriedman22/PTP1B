#Import Necessary Packages
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

def load_rmsf(path, file_name):
    a, r = [],[]
    with open('../../../' + path + '/rmsf_lig_' + file_name + '.xvg') as f:
        for _ in range(17):
            next(f)
        for line in f:
            cols = line.split()
            if len(cols) == 2:
                a.append(float(cols[0]))
                r.append(float(cols[1]))
    return a, r

#File paths for all input files
file_path_AD = ['AD_rebuild_a7/analysis', 'AD/analysis', '1sug_AD/analysis', '1sug_no_a7_AD/analysis', 'AD_dis/analysis/config7', 'AD_dis/analysis/config9', 'AD_dis/analysis/config11', '1sug_dis_AD/analysis/config7', '1sug_dis_AD/analysis/config9',  '1sug_dis_AD/analysis/config11',  '1sug_dis_AD/analysis/config11_2',  '1sug_dis_AD/analysis/config_alt',  '1sug_dis_AD/analysis/config_alt2']
file_name_AD = ['a7_AD', 'AD', '1sug_AD', '1sug_na7_AD', 'AD_dis7', 'AD_dis9', 'AD_dis11', '1sug_AD_dis7', '1sug_AD_dis9', '1sug_AD_dis11', '1sug_AD_dis11_2', '1sug_AD_dis_alt', '1sug_AD_dis_alt2']

file_path_BBR = ['BBR_a7/analysis', 'BBR_1sug/analysis', 'BBR_dis/analysis/config7', 'BBR_dis/analysis/config9', 'BBR_dis/analysis/config11', 'BBR_1sug_dis/analysis/config7', 'BBR_1sug_dis/analysis/config11']
file_name_BBR = ['BBR_a7', 'BBR_1sug', 'BBR_dis7', 'BBR_dis9', 'BBR_dis11', 'BBR_1sug_dis7', 'BBR_1sug_dis11']

#open all files
RMSF_AD = np.zeros((len(file_path_AD), 15))
RMSF_BBR = np.zeros((len(file_path_BBR), 41))
RMSF_mean_AD = np.zeros(len(file_path_AD))
RMSF_mean_BBR = np.zeros(len(file_path_BBR))
RMSF_err_AD = np.zeros(len(file_path_AD))
RMSF_err_BBR = np.zeros(len(file_path_BBR))


for i in range(len(file_path_AD)):
    #Load Data
    a, RMSF_AD[i][:] = load_rmsf(file_path_AD[i], file_name_AD[i])
    #Mean and SEM for each trajectory
    RMSF_mean_AD[i] = np.mean(RMSF_AD[i][:])
    RMSF_err_AD[i] = stats.sem(RMSF_AD[i][:])

for i in range(len(file_path_BBR)):
    #Load Data
    a, RMSF_BBR[i][:] = load_rmsf(file_path_BBR[i], file_name_BBR[i])
    #Mean and SEM for each trajectory
    RMSF_mean_BBR[i] = np.mean(RMSF_BBR[i][:])
    RMSF_err_BBR[i] = stats.sem(RMSF_BBR[i][:])

#Name Labels
Label_AD = ['Open Ordered a7', 'Open Absent a7', 'Closed Ordered a7', 'Closed Absent a7', 'Open Dis', 'Open Dis', 'Open Dis', 'Closed Dis7', 'Closed Dis9', 'Closed Dis11', 'Closed Dis11', 'Closed Dis alt', 'Closed Dis Alt']
Label_BBR = ['Open Ordered', 'Closed Ordered', 'Open Dis7', 'Open Dis9', 'Open Dis11', 'Closed Dis7', 'Closed Dis11']

num = np.linspace(1, len(Label_AD)+1, num = len(Label_AD))
fig = plt.figure(figsize=(18,10))
ax1 = fig.add_subplot(111)
ax1.set_title("Comparison of AD RMSF")
ax1.set_ylabel('RMSF(nm)')
ax1.bar(num, RMSF_mean_AD)
plt.xticks(num, Label_AD, fontsize=14)
plt.errorbar(num, RMSF_mean_AD, yerr= RMSF_err_AD, fmt='o', color='black')
leg = ax1.legend(loc='upper right')
fig.savefig('RMSF_AD_compare.png') 

num = np.linspace(1, len(Label_BBR)+1, num = len(Label_BBR))
fig = plt.figure(figsize=(18,10))
ax1 = fig.add_subplot(111)
ax1.set_title("Comparison of BBR RMSF")
ax1.set_ylabel('RMSF(nm)')
ax1.bar(num, RMSF_mean_BBR)
plt.xticks(num, Label_BBR, fontsize=14)
plt.errorbar(num, RMSF_mean_BBR, yerr= RMSF_err_BBR, fmt='o', color='black')
leg = ax1.legend(loc='upper right')
fig.savefig('RMSF_BBR_compare.png') 

