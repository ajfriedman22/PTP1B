#Import Necessary Packages
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from itertools import product
import seaborn as sns
import pandas as pd
import random
import sys

#Import custom modules
sys.path.insert(1,'/ocean/projects/cts160011p/afriedma/code/PTP1B/util')
import plot

def plot_mean(inter, err, hel1, hel2, p, p1, p2):
    num = [5, 10, 15, 20]
    Method = ['Apo Open', 'Apo Closed', 'AD', 'BBR']
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.set_title('Comparison of Interactions b/w ' + hel1 + ' ' + hel2, fontsize = 20) 
    ax1.set_ylabel('Mean Number of Interactions', fontsize = 18)
    ax1.bar(num, inter, color = ['black', 'gray', 'blue', 'red'], width=4.5)
    plt.errorbar(num, inter, yerr= err, fmt='o', color='black')
    plt.xticks(num, Method, fontsize=14)

    plot.error_bar(5, 15, inter[0], inter[2], p, 1, 'k')
    plot.error_bar(5, 20, inter[0], inter[3], p1, 1, 'k')
    plot.error_bar(5, 10, inter[0], inter[1], p2, 1, 'k')
    fig.savefig(hel1 + '_' + hel2 + '_inter.png')
    plt.close(fig)

def plot_box(d_Apo_open, d_Apo_close, d_AD, d_BBR, inter1, inter2, p, p1, ylim):
    d_Apo_open_df = pd.DataFrame({'Apo Open': d_Apo_open})
    d_Apo_close_df = pd.DataFrame({'Apo Closed': d_Apo_close})
    d_AD_df = pd.DataFrame({'AD': d_AD})
    d_BBR_df = pd.DataFrame({'BBR': d_BBR})
    mean = np.array([np.mean(d_Apo_open), np.mean(d_Apo_close), np.mean(d_AD), np.mean(d_BBR)])

    df = pd.concat([d_Apo_open_df, d_Apo_close_df, d_AD_df, d_BBR_df])

    ax = sns.stripplot(data = df, dodge=True, alpha=0.05, zorder=1, palette='bright')
    ax = sns.pointplot(data = df, join=False, scale=0.75, palette='dark')
    
    plot.error_bar(0, 2, mean[0], mean[2], p, 1, 'k')
    plot.error_bar(0, 3, mean[0], mean[3], p1, 1, 'k')

    plt.ylabel('Mean # of Interactions', fontsize = 14)
    plt.xticks(fontsize = 12)
    plt.yticks(fontsize = 12)
    ax.set_ylim(0,ylim)
    if inter2 == 'L11':
        plt.title(r'Helical Interactions b/w $\alpha$-' + inter1 + ' and ' + inter2, fontsize = 18)
        plt.savefig('Hel_inter_a' + inter1 + '_' + inter2 + '_box.png')
    else:
        plt.title(r'Helical Interactions b/w $\alpha$-' + inter1 + r' and $\alpha$-' + inter2, fontsize = 18)
        plt.savefig('Hel_inter_a' + inter1 + '_a' + inter2 + '_box.png')
    plt.close()

#Make open arrays for time and atomic interances
da3_a6_AD, da3_a6_BBR = [],[]
da3_a6_Apo_open, da3_a6_Apo_close = [],[]

da3_a7_AD, da3_a7_BBR = [],[]
da3_a7_Apo_open, da3_a7_Apo_close = [],[]

da6_a7_AD, da6_a7_BBR = [],[]
da6_a7_Apo_open, da6_a7_Apo_close = [],[]

dL11_a7_AD, dL11_a7_BBR = [],[]
dL11_a7_Apo_open, dL11_a7_Apo_close = [],[]

#List of all directory paths for each group
dir_path_Apo_open = ['rebuild_a7_high/config9/analysis', 'rebuild_a7_high/config11/analysis', 'Apo_dis/analysis']
dir_path_Apo_close = ['Apo_1SUG/analysis/1sug', 'Apo_1SUG/analysis/1sug2']
dir_path_AD = ['mutate/WT/AD/analysis', '1sug_dis_AD/analysis/config11', '1sug_dis_AD/analysis/config_alt', '1sug_dis_AD/analysis/config_alt2']
dir_path_BBR = ['mutate/WT/BBR/analysis', 'BBR_a7/analysis']


#List interactions of interest
inters = ['a3_a6', 'a7_a3', 'a7_a6', 'a7_L11']

#Load all data
for i in range(len(dir_path_Apo_open)):
    for j in range(len(inters)):
        data = open('../../../' + dir_path_Apo_open[i] + '/' + inters[j] + '_inter.txt').readlines()
        if j == 0:
            for k in data:
                da3_a6_Apo_open.append(float(k))
        if j == 1:
            for k in data:
                da3_a7_Apo_open.append(float(k))
        if j == 2:
            for k in data:
                da6_a7_Apo_open.append(float(k))
        if j == 3:
            for k in data:
                dL11_a7_Apo_open.append(float(k))

for i in range(len(dir_path_Apo_close)):
    for j in range(len(inters)):
        data = open('../../../' + dir_path_Apo_close[i] + '/' + inters[j] + '_inter.txt').readlines()
        if j == 0:
            for k in data:
                da3_a6_Apo_close.append(float(k))
        if j == 1:
            for k in data:
                da3_a7_Apo_close.append(float(k))
        if j == 2:
            for k in data:
                da6_a7_Apo_close.append(float(k))
        if j == 3:
            for k in data:
                dL11_a7_Apo_close.append(float(k))

for i in range(len(dir_path_AD)):
    for j in range(len(inters)):
        data = open('../../../' + dir_path_AD[i] + '/' + inters[j] + '_inter.txt').readlines()
        if j == 0:
            for k in data:
                da3_a6_AD.append(float(k))
        if j == 1:
            for k in data:
                da3_a7_AD.append(float(k))
        if j == 2:
            for k in data:
                da6_a7_AD.append(float(k))
        if j == 3:
            for k in data:
                dL11_a7_AD.append(float(k))

for i in range(len(dir_path_BBR)):
    for j in range(len(inters)):
        data = open('../../../' + dir_path_BBR[i] + '/' + inters[j] + '_inter.txt').readlines()
        if j == 0:
            for k in data:
                da3_a6_BBR.append(float(k))
        if j == 1:
            for k in data:
                da3_a7_BBR.append(float(k))
        if j == 2:
            for k in data:
                da6_a7_BBR.append(float(k))
        if j == 3:
            for k in data:
                dL11_a7_BBR.append(float(k))

#Calculate mean and sem for interactions
da3_a6 = np.array([np.mean(da3_a6_Apo_open), np.mean(da3_a6_Apo_close), np.mean(da3_a6_AD), np.mean(da3_a6_BBR)])
da3_a7 = np.array([np.mean(da3_a7_Apo_open), np.mean(da3_a7_Apo_close), np.mean(da3_a7_AD), np.mean(da3_a7_BBR)])
da6_a7 = np.array([np.mean(da6_a7_Apo_open), np.mean(da6_a7_Apo_close), np.mean(da6_a7_AD), np.mean(da6_a7_BBR)])
dL11_a7 = np.array([np.mean(dL11_a7_Apo_open), np.mean(dL11_a7_Apo_close), np.mean(dL11_a7_AD), np.mean(dL11_a7_BBR)])

da3_a6_err = np.array([stats.sem(da3_a6_Apo_open), stats.sem(da3_a6_Apo_close), stats.sem(da3_a6_AD), stats.sem(da3_a6_BBR)])
da3_a7_err = np.array([stats.sem(da3_a7_Apo_open), stats.sem(da3_a7_Apo_close), stats.sem(da3_a7_AD), stats.sem(da3_a7_BBR)])
da6_a7_err = np.array([stats.sem(da6_a7_Apo_open), stats.sem(da6_a7_Apo_close), stats.sem(da6_a7_AD), stats.sem(da6_a7_BBR)])
dL11_a7_err = np.array([stats.sem(dL11_a7_Apo_open), stats.sem(dL11_a7_Apo_close), stats.sem(dL11_a7_AD), stats.sem(dL11_a7_BBR)])

#Run t-test between groups
st, p1 = stats.ttest_ind(da3_a6_Apo_open, da3_a6_AD, equal_var = False) #Welch's t-test b/w Apo Open + AD
st, p2 = stats.ttest_ind(da3_a6_Apo_open, da3_a6_BBR, equal_var = False) #Welch's t-test b/w Apo Open + BBR
st, p3 = stats.ttest_ind(da3_a6_Apo_open, da3_a6_Apo_close, equal_var = False) #Welch's t-test b/w Apo Open + closed

st, p4 = stats.ttest_ind(da3_a7_Apo_open, da3_a7_AD, equal_var = False) #Welch's t-test b/w Apo Open + AD
st, p5 = stats.ttest_ind(da3_a7_Apo_open, da3_a7_BBR, equal_var = False) #Welch's t-test b/w Apo Open + BBR
st, p6 = stats.ttest_ind(da3_a7_Apo_open, da3_a7_Apo_close, equal_var = False) #Welch's t-test b/w Apo Open + closed

st, p7 = stats.ttest_ind(da6_a7_Apo_open, da6_a7_AD, equal_var = False) #Welch's t-test b/w Apo Open + AD
st, p8 = stats.ttest_ind(da6_a7_Apo_open, da6_a7_BBR, equal_var = False) #Welch's t-test b/w Apo Open + BBR
st, p9 = stats.ttest_ind(da6_a7_Apo_open, da6_a7_Apo_close, equal_var = False) #Welch's t-test b/w Apo Open + closed

st, p10 = stats.ttest_ind(dL11_a7_Apo_open, dL11_a7_AD, equal_var = False) #Welch's t-test b/w Apo Open + AD
st, p11 = stats.ttest_ind(dL11_a7_Apo_open, dL11_a7_BBR, equal_var = False) #Welch's t-test b/w Apo Open + BBR
st, p12 = stats.ttest_ind(dL11_a7_Apo_open, dL11_a7_Apo_close, equal_var = False) #Welch's t-test b/w Apo Open + closed

plot_mean(da3_a6, da3_a6_err, 'a3', 'a6', p1, p2, p3)
plot_mean(da3_a7, da3_a7_err, 'a3', 'a7', p4, p5, p6)
plot_mean(da6_a7, da6_a7_err, 'a6', 'a7', p7, p8, p9)
plot_mean(dL11_a7, dL11_a7_err, 'L11', 'a7', p10, p11, p12)

#Create box plots for the data frames
plot_box(da3_a6_Apo_open, da3_a6_Apo_close, da3_a6_AD, da3_a6_BBR, '3', '6', p1, p2, 35)
plot_box(da3_a7_Apo_open, da3_a7_Apo_close, da3_a7_AD, da3_a7_BBR, '3', '7', p4, p5, 30)
plot_box(da6_a7_Apo_open, da6_a7_Apo_close, da6_a7_AD, da6_a7_BBR, '6', '7', p7, p8, 20)
plot_box(dL11_a7_Apo_open, dL11_a7_Apo_close, dL11_a7_AD, dL11_a7_BBR, '7', 'L11', p10, p11, 20)

output = open('a3_a6.txt', 'w')
output.write('Apo Open:' + str(da3_a6[0]) + '+/-' + str(da3_a6_err[0]) + '\n')
output.write('Apo Closed:' + str(da3_a6[1]) + '+/-' + str(da3_a6_err[1]) + '\n')
output.write('AD:' + str(da3_a6[2]) + '+/-'  + str(da3_a6_err[2]) + '\n')
output.write('BBR:' + str(da3_a6[3]) + '+/-'  + str(da3_a6_err[3]) + '\n')

output2 = open('a3_a7.txt', 'w')
output2.write('Apo Open:' + str(da3_a7[0]) + '+/-' + str(da3_a7_err[0]) + '\n')
output2.write('Apo Closed:' + str(da3_a7[1]) + '+/-' + str(da3_a7_err[1]) + '\n')
output2.write('AD:' + str(da3_a7[2]) + '+/-' + str(da3_a7_err[2]) + '\n')
output2.write('BBR:' + str(da3_a7[3]) + '+/-' + str(da3_a7_err[3]) + '\n')

output3 = open('a6_a7.txt', 'w')
output3.write('Apo Open:' + str(da6_a7[0]) + '+/-' + str(da6_a7_err[0]) + '\n')
output3.write('Apo Closed:' + str(da6_a7[1]) + '+/-' + str(da6_a7_err[1]) + '\n')
output3.write('AD:' + str(da6_a7[2]) + '+/-' + str(da6_a7_err[2]) + '\n')
output3.write('BBR:' + str(da6_a7[3]) + '+/-' + str(da6_a7_err[3]) + '\n')

output4 = open('a3_WPD.txt', 'w')
output4.write('Apo Open:' + str(dL11_a7[0]) + '+/-' + str(dL11_a7_err[0]) + '\n')
output4.write('Apo Closed:' + str(dL11_a7[1]) + '+/-' + str(dL11_a7_err[1]) + '\n')
output4.write('AD:' + str(dL11_a7[2]) + '+/-' + str(dL11_a7_err[2]) + '\n')
output4.write('BBR:' + str(dL11_a7[3]) + '+/-' + str(dL11_a7_err[3]) + '\n')

#Compute percentage difference from mean number of interactions for all three helices relative to the Apo closed State
output = open('Per_diff_open.txt', 'w')
label = ['Apo Open', 'AD', 'BBR']
label_open = ['Apo Closed', 'AD', 'BBR']
helices = ['a3 and a6', 'a3 and a7', 'a6 and a7', 'L-11 and a7']
per_diff_all = np.zeros((len(helices), len(label)))
index = [0, 2, 3]
index_open = [1, 2, 3]
for i in range(len(index)):
    n = index[i]
    per_diff_all[0][i] = ((da3_a6[n] - da3_a6[1]) / ((da3_a6[n] + da3_a6[1])/2)) * 100
    per_diff_all[1][i] = ((da3_a7[n] - da3_a7[1]) / ((da3_a7[n] + da3_a7[1])/2)) * 100
    per_diff_all[2][i] = ((da6_a7[n] - da6_a7[1]) / ((da6_a7[n] + da6_a7[1])/2)) * 100
    per_diff_all[3][i] = ((dL11_a7[n] - dL11_a7[1]) / ((dL11_a7[n] + dL11_a7[1])/2)) * 100
    m = index_open[i]
    output.write(label_open[i] + '\n')
    per_diff_open = ((da3_a6[m] - da3_a6[0]) / ((da3_a6[m] + da3_a6[0])/2)) * 100
    output.write('a3 and a6: ' + str(per_diff_open) + '\n')
    per_diff_open = ((da3_a7[m] - da3_a7[0]) / ((da3_a7[m] + da3_a7[0])/2)) * 100
    output.write('a3 and a7: ' + str(per_diff_open) + '\n')
    per_diff_open = ((da6_a7[m] - da6_a7[0]) / ((da6_a7[m] + da6_a7[0])/2)) * 100
    output.write('a6 and a7: ' + str(per_diff_open) + '\n')
    per_diff_open = ((dL11_a7[m] - dL11_a7[0]) / ((dL11_a7[m] + dL11_a7[0])/2)) * 100
    output.write('a7 and L11: ' + str(per_diff_open) + '\n')
print(per_diff_all)

#Make heatmap to show the changes in helical interactions relative to the Apo Closed State
ax = plt.figure(figsize=(8, 10), frameon=True) # no visible frame
ax = sns.heatmap(per_diff_all, annot=False, cmap = 'PuBu_r', cbar = True, vmin = -100, vmax = 0, cbar_kws={'label': 'Percentage Difference from Apo Closed'}, xticklabels = label, yticklabels = helices)
#ax.add_artist(lines.Line2D([0, 20], [7, 7], color = 'black', linestyle= '--', linewidth = 4))
plt.title('Helical Interactions Compared to Apo Closed WPD Loop')
plt.savefig('Hel_inter_cmpr.png')
plt.close()


