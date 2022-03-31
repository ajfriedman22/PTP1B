#Import Necessary Packages
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from itertools import product
def plot_func(dist, err, label, p, p1, p2):
    num = [5, 10, 15, 20]
    Method = ['Apo Open', 'Apo Closed', 'AD', 'BBR']
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.set_title('Comparison of Distance b/w ' + label, fontsize = 20) 
    ax1.set_ylabel('Distance b/w Residues', fontsize = 18)
    ax1.bar(num, dist, color = ['black', 'black', 'darkblue', 'darkred', 'blue', 'red'], width=4.5)
    plt.errorbar(num, dist, yerr= err, fmt='o', color='black')
    plt.xticks(num, Method, fontsize=14)
    if p < 0.05 and p > 0.01:
        x1, x2 = 5, 15 #Columns for Apo and AD
        y, h, col = (1.1*dist[[0, 2]].max()), 1, 'b'
        plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
        plt.text((x1+x2)*0.5, y+h, "*" , ha='center', va='bottom', color=col)
    if p < 0.01 and p > 0.001:
        x1, x2 = 5, 15 #Columns for Apo and AD
        y, h, col = (1.1*dist[[0, 2]].max()), 1, 'b'
        plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
        plt.text((x1+x2)*0.5, y+h, "**" , ha='center', va='bottom', color=col)
    if p < 0.001:
        x1, x2 = 5, 15 #Columns for Apo and AD
        y, h, col = (1.1*dist[[0, 2]].max()), 1, 'b'
        plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
        plt.text((x1+x2)*0.5, y+h, "***" , ha='center', va='bottom', color=col)
    if p1 < 0.05 and p1 > 0.01:
        x1, x2 = 5, 20 #Columns for Apo and BBR
        y, h, col = (1.1*dist[[0, 3]].max()), 1, 'r'
        plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
        plt.text((x1+x2)*0.5, y+h, "*" , ha='center', va='bottom', color=col)
    if p1 < 0.01 and p1 > 0.001:
        x1, x2 = 5, 20 #Columns for Apo and BBR
        y, h, col = (1.1*dist[[0, 3]].max()), 1, 'r'
        plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
        plt.text((x1+x2)*0.5, y+h, "**" , ha='center', va='bottom', color=col)
    if p1 < 0.001:
        x1, x2 = 5, 20 #Columns for Apo and BBR
        y, h, col = (1.1*dist[[0, 3]].max()), 1, 'r'
        plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
        plt.text((x1+x2)*0.5, y+h, "***" , ha='center', va='bottom', color=col)
    if p2 < 0.05 and p2 > 0.01:
        x1, x2 = 5, 10 #Columns for Apo and AD
        y, h, col = (1.1*dist[[0, 1]].max()), 1, 'k'
        plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
        plt.text((x1+x2)*0.5, y+h, "*" , ha='center', va='bottom', color=col)
    if p2 < 0.01 and p2 > 0.001:
        x1, x2 = 5, 10 #Columns for Apo and AD
        y, h, col = (1.1*dist[[0, 1]].max()), 1, 'k'
        plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
        plt.text((x1+x2)*0.5, y+h, "**" , ha='center', va='bottom', color=col)
    if p2 < 0.001:
        x1, x2 = 5, 10 #Columns for Apo and AD
        y, h, col = (1.1*dist[[0, 1]].max()), 1, 'k'
        plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
        plt.text((x1+x2)*0.5, y+h, "***" , ha='center', va='bottom', color=col)

    fig.savefig(label + '_inter.png')
    plt.close(fig)

#Make open arrays for time and atomic distances
pairs = ['151_191_', '152_297_', '178_150_', '179_191_', '185_191_', '189_295_', '200_282_', '200_287_', '264_185_', '276_292_', '280_287_']
label = ['L11_a3', 'L11_a7', 'WPD_L11', 'WPD_a3_top', 'WPD_a3', 'a3_a7_top', 'a3_a6', 'a3_a7', 'a3_a6_top', 'a6_a7_top', 'a6_a7']
eq_time = [5, 5, 50, 5, 75, 5, 5, 5, 70, 5, 60, 30]
tot_time = [200, 200, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300]

for i in range(len(pairs)):
    d_AD, d_BBR = [],[]
    d_Apo_open, d_Apo_close = [],[]
    d_WT_AD, d_WT_BBR, d_a7, d_dis9, d_dis11, d_1sug, d_1sug2, d_1sug3 = [],[],[],[],[],[],[],[]
    d_BBR_a7, d_AD_dis11, d_AD_alt, d_AD_alt2 = [],[],[],[]

    #Input Data for a3 and a6 residue dist
    for j in open('../../../mutate/WT/AD/analysis/' + pairs[i] + 'dist.txt').readlines():
        d_WT_AD.append(float(j))
    for j in open('../../../mutate/WT/BBR/analysis/' + pairs[i] + 'dist.txt').readlines():
        d_WT_BBR.append(float(j))
    for j in open('../../../rebuild_a7/analysis/' + pairs[i] + 'dist.txt').readlines():
        d_a7.append(float(j))
    for j in open('../../../rebuild_a7_high/config9/analysis/' + pairs[i] + 'dist.txt').readlines():
        d_dis9.append(float(j))
    for j in open('../../../rebuild_a7_high/config11/analysis/' + pairs[i] + 'dist.txt').readlines():
        d_dis11.append(float(j))
    for j in open('../../../Apo_1SUG/analysis/1sug/' + pairs[i] + 'dist.txt').readlines():
        d_1sug.append(float(j))
    for j in open('../../../Apo_1SUG/analysis/1sug2/' + pairs[i] + 'dist.txt').readlines():
        d_1sug2.append(float(j))
    for j in open('../../../Apo_1SUG/analysis/1sug3/' + pairs[i] + 'dist.txt').readlines():
        d_1sug3.append(float(j))
    for j in open('../../../BBR_a7/analysis/' + pairs[i] + 'dist.txt').readlines():
        d_BBR_a7.append(float(j))
    for j in open('../../../1sug_dis_AD/analysis/config11/' + pairs[i] + 'dist.txt').readlines():
        d_AD_dis11.append(float(j))
    for j in open('../../../1sug_dis_AD/analysis/config_alt/' + pairs[i] + 'dist.txt').readlines():
        d_AD_alt.append(float(j))
    for j in open('../../../1sug_dis_AD/analysis/config_alt2/' + pairs[i] + 'dist.txt').readlines():
        d_AD_alt2.append(float(j))
    
    d_Apo_open = np.concatenate((d_dis9, d_dis11))
    d_Apo_close = np.concatenate((d_1sug, d_1sug2, d_1sug3))
    d_AD = np.concatenate((d_AD_dis11, d_AD_alt, d_AD_alt2))
    d_BBR = d_BBR_a7

    #Calculate mean and sem for interactions
    d = np.array([np.mean(d_Apo_open), np.mean(d_Apo_close), np.mean(d_AD), np.mean(d_BBR)])

    d_err = np.array([stats.sem(d_Apo_open), stats.sem(d_Apo_close), stats.sem(d_AD), stats.sem(d_BBR)])

    #Run t-test between groups
    st, p1 = stats.ttest_ind(d_Apo_open, d_AD, equal_var = False) #Welch's t-test b/w Apo Open + AD
    st, p2 = stats.ttest_ind(d_Apo_open, d_BBR, equal_var = False) #Welch's t-test b/w Apo Open + BBR
    st, p3 = stats.ttest_ind(d_Apo_open, d_Apo_close, equal_var = False) #Welch's t-test b/w Apo Open + closed

    plot_func(d, d_err, label[i], p1, p2, p3)

    output = open(label[i] + '.txt', 'w')
    output.write('Apo Open:' + str(d[0]) + '+/-' + str(d_err[0]) + '\n')
    output.write('Apo Closed:' + str(d[1]) + '+/-' + str(d_err[1]) + '\n')
    output.write('AD:' + str(d[2]) + '+/-'  + str(d_err[2]) + '\n')
    output.write('BBR:' + str(d[3]) + '+/-'  + str(d_err[3]) + '\n')
    output.close()

