#!/ usr / bin / env python

#Import packages
from matplotlib import pyplot as plt
import numpy as np
from scipy import stats
import sys

#Import custom modules
sys.path.insert(1, '/ocean/projects/cts160011p/afriedma/code/PTP1B/util/')
import mdfunc
import plot

def load_data(dir_path):
    data = open('../../../' + dir_path + 'Hbond_allo_per.txt', 'r').readlines()
    
    bonds, per = [],[]
    for i in range(len(data)):
        if i%2 == 0:
            bonds.append(data[i])
        else:
            per.append(float(data[i].strip(' :\n')))
    return bonds, per

#Import data
dir_path = ['rebuild_a7/analysis/', 'AD_rebuild_a7/analysis/', 'Apo_1SUG/analysis/1sug/', 'Apo_1SUG/analysis/1sug2/', 'Apo_1SUG/analysis/1sug3/', '1sug_AD/analysis/', 
        '1sug_dis_AD/analysis/config7/', '1sug_dis_AD/analysis/config9/', '1sug_dis_AD/analysis/config11/', '1sug_dis_AD/analysis/config11_2/', '1sug_dis_AD/analysis/config_alt/', 
        '1sug_dis_AD/analysis/config_alt2/', 'AD_dis/analysis/config7/', 'AD_dis/analysis/config9/', 'AD_dis/analysis/config11/', '1sug_dis/analysis/config7/', '1sug_dis/analysis/config9/',
        '1sug_dis/analysis/config11/', 'rebuild_a7_high/config7/analysis/', 'rebuild_a7_high/config9/analysis/', 'rebuild_a7_high/config11/analysis/', 'BBR_a7/analysis/', 
        'BBR_1sug/analysis/', 'AD_BBR/analysis/', 'mutate/WT/AD/analysis/', 'mutate/WT/BBR/analysis/', 'Apo_dis/analysis/']

bonds, per_a7 = load_data(dir_path[0])

per_all = np.zeros((len(dir_path), len(per_a7)))

for i in range(len(dir_path)):
    bonds, per_all[i,:] = load_data(dir_path[i])

output = open('Hbond_network.txt', 'w')
#Plot data for all bonds
for i in range(len(bonds)):
    apo_open = np.array([per_all[18, i], per_all[19, i], per_all[20, i], per_all[26, i]])
    apo_closed = np.array([per_all[2, i], per_all[3, i]])
    AD = np.array([per_all[8, i], per_all[10, i], per_all[11, i], per_all[24, i]])
    BBR = np.array([per_all[21, i], per_all[25, i]])

    mean = [np.mean(apo_open), np.mean(apo_closed), np.mean(AD), np.mean(BBR)]
    error = [stats.sem(apo_open), stats.sem(apo_closed), stats.sem(AD), stats.sem(BBR)]

    #Detemine p-value
    st, p = stats.ttest_ind(apo_open, apo_closed, equal_var = False) #Welch's t-test b/w open apo and closed apo
    st, p2 = stats.ttest_ind(apo_open, AD, equal_var = False) #Welch's t-test b/w open apo and AD
    st, p3 = stats.ttest_ind(apo_open, BBR, equal_var = False) #Welch's t-test b/w open apo and BBR

    #Write the hbonds to file
    output.write(str(bonds[i]) + ': ' + str(p) + ' ' + str(p2) + ' ' + str(p3) + '\n')
    num = [1, 2, 3, 4]
    Method = ['Apo\n Open', 'Apo\n Closed', 'AD', 'BBR']

    #Plot Bar graph comparing averages
    fig = plt.figure()
    plt.title('Time Hbond ' + str(bonds[i]) + ' is formed')    
    plt.ylabel('% Time Hbond Formed')
    plt.bar(num, mean, color = ['gray', 'red', 'blue', 'purple'], width=0.9)
    plt.errorbar(num, mean, yerr=error, fmt='o', color='black')
    plt.xticks(num, Method, fontsize=8)
    plot.error_bar(1, 2, mean[0], mean[1], p, 5, 'black')
    plot.error_bar(1, 3, mean[0], mean[2], p2, 5, 'black')
    plot.error_bar(1, 4, mean[0], mean[3], p3, 5, 'black')
    fig.savefig(str(bonds[i]) + '_group_percent.png')
    plt.close(fig)

