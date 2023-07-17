#Import Necessary Packages
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import sys

#Import custom modules
sys.path.insert(1,'/ocean/projects/cts160011p/afriedma/code/PTP1B/util')
import plot
import uncorr


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
st, p2 = stats.ttest_ind(AD_solo, BBR_solo, equal_var = False) #Welch's t-test b/w solo AD vs BBR

#Calculate mean and sem
lig_mean = [np.mean(AD_solo), np.mean(BBR_solo)]
lig_sem = [stats.sem(AD_solo), stats.sem(BBR_solo)]
lig_mean_combo = [np.mean(AD_combo), np.mean(BBR_combo)]
lig_sem_combo = [stats.sem(AD_combo), stats.sem(BBR_combo)]

num = [1, 2, 3, 4]
lig = ['AD', 'Ternary AD', 'BBR', 'Ternary BBR']
num1 = [1,3]
num2 = [2,4]

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

