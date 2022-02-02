#Import Necessary Packages
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from itertools import product
import seaborn as sns

def plot_func(inter, err, hel1, hel2, p, p1, p2):
    num = [5, 10, 15, 20]
    Method = ['Apo Open', 'Apo Closed', 'AD', 'BBR']
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.set_title('Comparison of Interactions b/w ' + hel1 + ' ' + hel2, fontsize = 20) 
    ax1.set_ylabel('Distance b/w Residues', fontsize = 18)
    ax1.bar(num, inter, color = ['black', 'gray', 'blue', 'red'], width=4.5)
    plt.errorbar(num, inter, yerr= err, fmt='o', color='black')
    plt.xticks(num, Method, fontsize=14)
    if p < 0.05 and p > 0.01:
        x1, x2 = 5, 15 #Columns for Apo and AD
        y, h, col = (1.1*inter[[0, 2]].max()), 1, 'b'
        plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
        plt.text((x1+x2)*0.5, y+h, "*" , ha='center', va='bottom', color=col)
    if p < 0.01 and p > 0.001:
        x1, x2 = 5, 15 #Columns for Apo and AD
        y, h, col = (1.1*inter[[0, 2]].max()), 1, 'b'
        plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
        plt.text((x1+x2)*0.5, y+h, "**" , ha='center', va='bottom', color=col)
    if p < 0.001:
        x1, x2 = 5, 15 #Columns for Apo and AD
        y, h, col = (1.1*inter[[0, 2]].max()), 1, 'b'
        plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
        plt.text((x1+x2)*0.5, y+h, "***" , ha='center', va='bottom', color=col)
    if p1 < 0.05 and p1 > 0.01:
        x1, x2 = 5, 20 #Columns for Apo and BBR
        y, h, col = (1.1*inter[[0, 3]].max()), 1, 'r'
        plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
        plt.text((x1+x2)*0.5, y+h, "*" , ha='center', va='bottom', color=col)
    if p1 < 0.01 and p1 > 0.001:
        x1, x2 = 5, 20 #Columns for Apo and BBR
        y, h, col = (1.1*inter[[0, 3]].max()), 1, 'r'
        plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
        plt.text((x1+x2)*0.5, y+h, "**" , ha='center', va='bottom', color=col)
    if p1 < 0.001:
        x1, x2 = 5, 20 #Columns for Apo and BBR
        y, h, col = (1.1*inter[[0, 3]].max()), 1, 'r'
        plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
        plt.text((x1+x2)*0.5, y+h, "***" , ha='center', va='bottom', color=col)
    if p2 < 0.05 and p2 > 0.01:
        x1, x2 = 5, 10 #Columns for Apo and AD
        y, h, col = (1.1*inter[[0, 1]].max()), 1, 'k'
        plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
        plt.text((x1+x2)*0.5, y+h, "*" , ha='center', va='bottom', color=col)
    if p2 < 0.01 and p2 > 0.001:
        x1, x2 = 5, 10 #Columns for Apo and AD
        y, h, col = (1.1*inter[[0, 1]].max()), 1, 'k'
        plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
        plt.text((x1+x2)*0.5, y+h, "**" , ha='center', va='bottom', color=col)
    if p2 < 0.001:
        x1, x2 = 5, 10 #Columns for Apo and AD
        y, h, col = (1.1*inter[[0, 1]].max()), 1, 'k'
        plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
        plt.text((x1+x2)*0.5, y+h, "***" , ha='center', va='bottom', color=col)
    fig.savefig(hel1 + '_' + hel2 + '_inter.png')
    plt.close(fig)

def deter_corr_time(inter, eq_time, tot_time):
    corr_time = []
    #WT AD
    ref = inter[0]
    t_ref = 0
    time = np.linspace(eq_time, tot_time, num=len(inter))
    for i in range(len(inter)):
        if abs(inter[i] - ref) > 2:
            corr_time.append(time[i] - t_ref)
            ref = inter[i]
            t_ref = time[i]
    return np.mean(corr_time)
def corr_array(inter, corr_mean, eq_time, tot_time, new):
    #Determine how many frames correspond to the correlation time
    frame = np.round((corr_mean / (tot_time - eq_time)) * len(inter))

    #Load frames to overall list array
    for i in range(len(inter)):
        if i%frame == 0:
            new.append(inter[i])
    return new

#Make open arrays for time and atomic interances
da3_a6_AD, da3_a6_BBR = [],[]
da3_a6_Apo_open, da3_a6_Apo_close = [],[]

da3_a7_AD, da3_a7_BBR = [],[]
da3_a7_Apo_open, da3_a7_Apo_close = [],[]

da6_a7_AD, da6_a7_BBR = [],[]
da6_a7_Apo_open, da6_a7_Apo_close = [],[]

dL11_a7_AD, dL11_a7_BBR = [],[]
dL11_a7_Apo_open, dL11_a7_Apo_close = [],[]

da3_a6_WT_AD, da3_a6_WT_BBR, da3_a6_a7, da3_a6_dis9, da3_a6_dis11, da3_a6_1sug, da3_a6_1sug2, da3_a6_1sug3 = [],[],[],[],[],[],[],[]
da3_a6_BBR_a7, da3_a6_AD_dis11, da3_a6_AD_alt, da3_a6_AD_alt2 = [],[],[],[]

da3_a7_WT_AD, da3_a7_WT_BBR, da3_a7_a7, da3_a7_dis9, da3_a7_dis11, da3_a7_1sug, da3_a7_1sug2, da3_a7_1sug3 = [],[],[],[],[],[],[],[]
da3_a7_BBR_a7, da3_a7_AD_dis11, da3_a7_AD_alt, da3_a7_AD_alt2 = [],[],[],[]

da6_a7_WT_AD, da6_a7_WT_BBR, da6_a7_a7, da6_a7_dis9, da6_a7_dis11, da6_a7_1sug, da6_a7_1sug2, da6_a7_1sug3 = [],[],[],[],[],[],[],[]
da6_a7_BBR_a7, da6_a7_AD_dis11, da6_a7_AD_alt, da6_a7_AD_alt2 = [],[],[],[]

dL11_a7_WT_AD, dL11_a7_WT_BBR, dL11_a7_a7, dL11_a7_dis9, dL11_a7_dis11, dL11_a7_1sug, dL11_a7_1sug2, dL11_a7_1sug3 = [],[],[],[],[],[],[],[]
dL11_a7_BBR_a7, dL11_a7_AD_dis11, dL11_a7_AD_alt, dL11_a7_AD_alt2 = [],[],[],[]

#Input Data for a3 and a6 residue inter
for i in open("../../../mutate/WT/AD/analysis/a3_a6_inter.txt").readlines():
    da3_a6_WT_AD.append(float(i))
for i in open("../../../mutate/WT/BBR/analysis/a3_a6_inter.txt").readlines():
    da3_a6_WT_BBR.append(float(i))
for i in open("../../../rebuild_a7/analysis/a3_a6_inter.txt").readlines():
    da3_a6_a7.append(float(i))
for i in open("../../../rebuild_a7_high/config9/analysis/a3_a6_inter.txt").readlines():
    da3_a6_dis9.append(float(i))
for i in open("../../../rebuild_a7_high/config11/analysis/a3_a6_inter.txt").readlines():
    da3_a6_dis11.append(float(i))
for i in open("../../../Apo_1SUG/analysis/1sug/a3_a6_inter.txt").readlines():
    da3_a6_1sug.append(float(i))
for i in open("../../../Apo_1SUG/analysis/1sug2/a3_a6_inter.txt").readlines():
    da3_a6_1sug2.append(float(i))
for i in open("../../../Apo_1SUG/analysis/1sug3/a3_a6_inter.txt").readlines():
    da3_a6_1sug3.append(float(i))
for i in open("../../../BBR_a7/analysis/a3_a6_inter.txt").readlines():
    da3_a6_BBR_a7.append(float(i))
for i in open("../../../1sug_dis_AD/analysis/config11/a3_a6_inter.txt").readlines():
    da3_a6_AD_dis11.append(float(i))
for i in open("../../../1sug_dis_AD/analysis/config_alt/a3_a6_inter.txt").readlines():
    da3_a6_AD_alt.append(float(i))
for i in open("../../../1sug_dis_AD/analysis/config_alt2/a3_a6_inter.txt").readlines():
    da3_a6_AD_alt2.append(float(i))

#Determine approximate time scale for correlated motions
corr_time = []
corr_time.append(deter_corr_time(da3_a6_WT_AD, 5, 200))
corr_time.append(deter_corr_time(da3_a6_WT_BBR, 5, 200))
corr_time.append(deter_corr_time(da3_a6_a7, 50, 300))
corr_time.append(deter_corr_time(da3_a6_dis9, 5, 300))
corr_time.append(deter_corr_time(da3_a6_dis11, 75, 300))
corr_time.append(deter_corr_time(da3_a6_1sug, 5, 300))
corr_time.append(deter_corr_time(da3_a6_1sug2, 5, 300))
corr_time.append(deter_corr_time(da3_a6_1sug3, 5, 300))
corr_time.append(deter_corr_time(da3_a6_BBR_a7, 70, 300))
corr_time.append(deter_corr_time(da3_a6_AD_dis11, 5, 300))
corr_time.append(deter_corr_time(da3_a6_AD_alt, 60, 300))
corr_time.append(deter_corr_time(da3_a6_AD_alt2, 30, 300))

corr_mean = np.mean(corr_time)

#Limit arrays to only samples within the correlation time and convert to Angstrom
corr_array(da3_a6_WT_AD, corr_mean, 5, 200, da3_a6_AD)
corr_array(da3_a6_WT_BBR, corr_mean, 5, 200, da3_a6_BBR)
corr_array(da3_a6_a7, corr_mean, 50, 300, da3_a6_Apo_open)
corr_array(da3_a6_dis9, corr_mean, 5, 300, da3_a6_Apo_open)
corr_array(da3_a6_dis11, corr_mean, 75, 300, da3_a6_Apo_open)
corr_array(da3_a6_1sug, corr_mean, 5, 300, da3_a6_Apo_close)
corr_array(da3_a6_1sug2, corr_mean, 5, 300, da3_a6_Apo_close)
corr_array(da3_a6_1sug3, corr_mean, 5, 300, da3_a6_Apo_open)
corr_array(da3_a6_BBR_a7, corr_mean, 70, 300, da3_a6_BBR)
corr_array(da3_a6_AD_dis11, corr_mean, 5, 300, da3_a6_AD)
corr_array(da3_a6_AD_alt, corr_mean, 60, 300, da3_a6_AD)
corr_array(da3_a6_AD_alt2, corr_mean, 30, 300, da3_a6_AD)

#Input Data for a3 and a7 residue inter
for i in open("../../../mutate/WT/AD/analysis/a7_a3_inter.txt").readlines():
    da3_a7_WT_AD.append(float(i))
for i in open("../../../mutate/WT/BBR/analysis/a7_a3_inter.txt").readlines():
    da3_a7_WT_BBR.append(float(i))
for i in open("../../../rebuild_a7/analysis/a7_a3_inter.txt").readlines():
    da3_a7_a7.append(float(i))
for i in open("../../../rebuild_a7_high/config9/analysis/a7_a3_inter.txt").readlines():
    da3_a7_dis9.append(float(i))
for i in open("../../../rebuild_a7_high/config11/analysis/a7_a3_inter.txt").readlines():
    da3_a7_dis11.append(float(i))
for i in open("../../../Apo_1SUG/analysis/1sug/a7_a3_inter.txt").readlines():
    da3_a7_1sug.append(float(i))
for i in open("../../../Apo_1SUG/analysis/1sug2/a7_a3_inter.txt").readlines():
    da3_a7_1sug2.append(float(i))
for i in open("../../../Apo_1SUG/analysis/1sug3/a7_a3_inter.txt").readlines():
    da3_a7_1sug3.append(float(i))
for i in open("../../../BBR_a7/analysis/a7_a3_inter.txt").readlines():
    da3_a7_BBR_a7.append(float(i))
for i in open("../../../1sug_dis_AD/analysis/config11/a7_a3_inter.txt").readlines():
    da3_a7_AD_dis11.append(float(i))
for i in open("../../../1sug_dis_AD/analysis/config_alt/a7_a3_inter.txt").readlines():
    da3_a7_AD_alt.append(float(i))
for i in open("../../../1sug_dis_AD/analysis/config_alt2/a7_a3_inter.txt").readlines():
    da3_a7_AD_alt2.append(float(i))

#Determine approximate time scale for correlated motions
corr_time = []
corr_time.append(deter_corr_time(da3_a7_WT_AD, 5, 200))
corr_time.append(deter_corr_time(da3_a7_WT_BBR, 5, 200))
corr_time.append(deter_corr_time(da3_a7_a7, 50, 300))
corr_time.append(deter_corr_time(da3_a7_dis9, 5, 300))
corr_time.append(deter_corr_time(da3_a7_dis11, 75, 300))
corr_time.append(deter_corr_time(da3_a7_1sug, 5, 300))
corr_time.append(deter_corr_time(da3_a7_1sug2, 5, 300))
corr_time.append(deter_corr_time(da3_a7_1sug3, 5, 300))
corr_time.append(deter_corr_time(da3_a7_BBR_a7, 70, 300))
corr_time.append(deter_corr_time(da3_a7_AD_dis11, 5, 300))
corr_time.append(deter_corr_time(da3_a7_AD_alt, 60, 300))
corr_time.append(deter_corr_time(da3_a7_AD_alt2, 30, 300))

corr_mean = np.mean(corr_time)

#Limit arrays to only samples within the correlation time and convert to Angstrom
corr_array(da3_a7_WT_AD, corr_mean, 5, 200, da3_a7_AD)
corr_array(da3_a7_WT_BBR, corr_mean, 5, 200, da3_a7_BBR)
corr_array(da3_a7_a7, corr_mean, 50, 300, da3_a7_Apo_open)
corr_array(da3_a7_dis9, corr_mean, 5, 300, da3_a7_Apo_open)
corr_array(da3_a7_dis11, corr_mean, 75, 300, da3_a7_Apo_open)
corr_array(da3_a7_1sug, corr_mean, 5, 300, da3_a7_Apo_close)
corr_array(da3_a7_1sug2, corr_mean, 5, 300, da3_a7_Apo_close)
corr_array(da3_a7_1sug3, corr_mean, 5, 300, da3_a7_Apo_open)
corr_array(da3_a7_BBR_a7, corr_mean, 70, 300, da3_a7_BBR)
corr_array(da3_a7_AD_dis11, corr_mean, 5, 300, da3_a7_AD)
corr_array(da3_a7_AD_alt, corr_mean, 60, 300, da3_a7_AD)
corr_array(da3_a7_AD_alt2, corr_mean, 30, 300, da3_a7_AD)

#Input Data for a6 and a7 residue inter
for i in open("../../../mutate/WT/AD/analysis/a7_a6_inter.txt").readlines():
    da6_a7_WT_AD.append(float(i))
for i in open("../../../mutate/WT/BBR/analysis/a7_a6_inter.txt").readlines():
    da6_a7_WT_BBR.append(float(i))
for i in open("../../../rebuild_a7/analysis/a7_a6_inter.txt").readlines():
    da6_a7_a7.append(float(i))
for i in open("../../../rebuild_a7_high/config9/analysis/a7_a6_inter.txt").readlines():
    da6_a7_dis9.append(float(i))
for i in open("../../../rebuild_a7_high/config11/analysis/a7_a6_inter.txt").readlines():
    da6_a7_dis11.append(float(i))
for i in open("../../../Apo_1SUG/analysis/1sug/a7_a6_inter.txt").readlines():
    da6_a7_1sug.append(float(i))
for i in open("../../../Apo_1SUG/analysis/1sug2/a7_a6_inter.txt").readlines():
    da6_a7_1sug2.append(float(i))
for i in open("../../../Apo_1SUG/analysis/1sug3/a7_a6_inter.txt").readlines():
    da6_a7_1sug3.append(float(i))
for i in open("../../../BBR_a7/analysis/a7_a6_inter.txt").readlines():
    da6_a7_BBR_a7.append(float(i))
for i in open("../../../1sug_dis_AD/analysis/config11/a7_a6_inter.txt").readlines():
    da6_a7_AD_dis11.append(float(i))
for i in open("../../../1sug_dis_AD/analysis/config_alt/a7_a6_inter.txt").readlines():
    da6_a7_AD_alt.append(float(i))
for i in open("../../../1sug_dis_AD/analysis/config_alt2/a7_a6_inter.txt").readlines():
    da6_a7_AD_alt2.append(float(i))

#Determine approximate time scale for correlated motions
corr_time = []
corr_time.append(deter_corr_time(da6_a7_WT_AD, 5, 200))
corr_time.append(deter_corr_time(da6_a7_WT_BBR, 5, 200))
corr_time.append(deter_corr_time(da6_a7_a7, 50, 300))
corr_time.append(deter_corr_time(da6_a7_dis9, 5, 300))
corr_time.append(deter_corr_time(da6_a7_dis11, 75, 300))
corr_time.append(deter_corr_time(da6_a7_1sug, 5, 300))
corr_time.append(deter_corr_time(da6_a7_1sug2, 5, 300))
corr_time.append(deter_corr_time(da6_a7_1sug3, 5, 300))
corr_time.append(deter_corr_time(da6_a7_BBR_a7, 70, 300))
corr_time.append(deter_corr_time(da6_a7_AD_dis11, 5, 300))
corr_time.append(deter_corr_time(da6_a7_AD_alt, 60, 300))
corr_time.append(deter_corr_time(da6_a7_AD_alt2, 30, 300))

corr_mean = np.mean(corr_time)

#Limit arrays to only samples within the correlation time and convert to Angstrom
corr_array(da6_a7_WT_AD, corr_mean, 5, 200, da6_a7_AD)
corr_array(da6_a7_WT_BBR, corr_mean, 5, 200, da6_a7_BBR)
corr_array(da6_a7_a7, corr_mean, 50, 300, da6_a7_Apo_open)
corr_array(da6_a7_dis9, corr_mean, 5, 300, da6_a7_Apo_open)
corr_array(da6_a7_dis11, corr_mean, 75, 300, da6_a7_Apo_open)
corr_array(da6_a7_1sug, corr_mean, 5, 300, da6_a7_Apo_close)
corr_array(da6_a7_1sug2, corr_mean, 5, 300, da6_a7_Apo_close)
corr_array(da6_a7_1sug3, corr_mean, 5, 300, da6_a7_Apo_open)
corr_array(da6_a7_BBR_a7, corr_mean, 70, 300, da6_a7_BBR)
corr_array(da6_a7_AD_dis11, corr_mean, 5, 300, da6_a7_AD)
corr_array(da6_a7_AD_alt, corr_mean, 60, 300, da6_a7_AD)
corr_array(da6_a7_AD_alt2, corr_mean, 30, 300, da6_a7_AD)

#Input Data for a3 and WPD residue inter
for i in open("../../../mutate/WT/AD/analysis/a7_L11_inter.txt").readlines():
    dL11_a7_WT_AD.append(float(i))
for i in open("../../../mutate/WT/BBR/analysis/a7_L11_inter.txt").readlines():
    dL11_a7_WT_BBR.append(float(i))
for i in open("../../../rebuild_a7/analysis/a7_L11_inter.txt").readlines():
    dL11_a7_a7.append(float(i))
for i in open("../../../rebuild_a7_high/config9/analysis/a7_L11_inter.txt").readlines():
    dL11_a7_dis9.append(float(i))
for i in open("../../../rebuild_a7_high/config11/analysis/a7_L11_inter.txt").readlines():
    dL11_a7_dis11.append(float(i))
for i in open("../../../Apo_1SUG/analysis/1sug/a7_L11_inter.txt").readlines():
    dL11_a7_1sug.append(float(i))
for i in open("../../../Apo_1SUG/analysis/1sug2/a7_L11_inter.txt").readlines():
    dL11_a7_1sug2.append(float(i))
for i in open("../../../Apo_1SUG/analysis/1sug3/a7_L11_inter.txt").readlines():
    dL11_a7_1sug3.append(float(i))
for i in open("../../../BBR_a7/analysis/a7_L11_inter.txt").readlines():
    dL11_a7_BBR_a7.append(float(i))
for i in open("../../../1sug_dis_AD/analysis/config11/a7_L11_inter.txt").readlines():
    dL11_a7_AD_dis11.append(float(i))
for i in open("../../../1sug_dis_AD/analysis/config_alt/a7_L11_inter.txt").readlines():
    dL11_a7_AD_alt.append(float(i))
for i in open("../../../1sug_dis_AD/analysis/config_alt2/a7_L11_inter.txt").readlines():
    dL11_a7_AD_alt2.append(float(i))

#Determine approximate time scale for correlated motions
corr_time = []
corr_time.append(deter_corr_time(dL11_a7_WT_AD, 5, 200))
corr_time.append(deter_corr_time(dL11_a7_WT_BBR, 5, 200))
corr_time.append(deter_corr_time(dL11_a7_a7, 50, 300))
corr_time.append(deter_corr_time(dL11_a7_dis9, 5, 300))
corr_time.append(deter_corr_time(dL11_a7_dis11, 75, 300))
corr_time.append(deter_corr_time(dL11_a7_1sug, 5, 300))
corr_time.append(deter_corr_time(dL11_a7_1sug2, 5, 300))
corr_time.append(deter_corr_time(dL11_a7_1sug3, 5, 300))
corr_time.append(deter_corr_time(dL11_a7_BBR_a7, 70, 300))
corr_time.append(deter_corr_time(dL11_a7_AD_dis11, 5, 300))
corr_time.append(deter_corr_time(dL11_a7_AD_alt, 60, 300))
corr_time.append(deter_corr_time(dL11_a7_AD_alt2, 30, 300))

corr_mean = np.mean(corr_time)

#Limit arrays to only samples within the correlation time and convert to Angstrom
corr_array(dL11_a7_WT_AD, corr_mean, 5, 200, dL11_a7_AD)
corr_array(dL11_a7_WT_BBR, corr_mean, 5, 200, dL11_a7_BBR)
corr_array(dL11_a7_a7, corr_mean, 50, 300, dL11_a7_Apo_open)
corr_array(dL11_a7_dis9, corr_mean, 5, 300, dL11_a7_Apo_open)
corr_array(dL11_a7_dis11, corr_mean, 75, 300, dL11_a7_Apo_open)
corr_array(dL11_a7_1sug, corr_mean, 5, 300, dL11_a7_Apo_close)
corr_array(dL11_a7_1sug2, corr_mean, 5, 300, dL11_a7_Apo_close)
corr_array(dL11_a7_1sug3, corr_mean, 5, 300, dL11_a7_Apo_open)
corr_array(dL11_a7_BBR_a7, corr_mean, 70, 300, dL11_a7_BBR)
corr_array(dL11_a7_AD_dis11, corr_mean, 5, 300, dL11_a7_AD)
corr_array(dL11_a7_AD_alt, corr_mean, 60, 300, dL11_a7_AD)
corr_array(dL11_a7_AD_alt2, corr_mean, 30, 300, dL11_a7_AD)

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

plot_func(da3_a6, da3_a6_err, 'a3', 'a6', p1, p2, p3)
plot_func(da3_a7, da3_a7_err, 'a3', 'a7', p4, p5, p6)
plot_func(da6_a7, da6_a7_err, 'a6', 'a7', p7, p8, p9)
plot_func(dL11_a7, dL11_a7_err, 'L11', 'a7', p10, p11, p12)

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
label = ['Apo Open', 'AD', 'BBR']
helices = ['a3 and a6', 'a3 and a7', 'a6 and a7', 'L-11 and a7']
per_diff_all = np.zeros((len(helices), len(label)))
index = [0, 2, 3]
for i in range(len(index)):
    n = index[i]
    per_diff_all[0][i] = ((da3_a6[n] - da3_a6[1]) / ((da3_a6[n] + da3_a6[1])/2)) * 100
    per_diff_all[1][i] = ((da3_a7[n] - da3_a7[1]) / ((da3_a7[n] + da3_a7[1])/2)) * 100
    per_diff_all[2][i] = ((da6_a7[n] - da6_a7[1]) / ((da6_a7[n] + da6_a7[1])/2)) * 100
    per_diff_all[3][i] = ((dL11_a7[n] - dL11_a7[1]) / ((dL11_a7[n] + dL11_a7[1])/2)) * 100

#Make heatmap to show the changes in helical interactions relative to the Apo Closed State
ax = plt.figure(figsize=(8, 10), frameon=True) # no visible frame
ax = sns.heatmap(per_diff_all, annot=False, cmap = 'bwr_r', cbar = True, vmin = -100, vmax = 100, cbar_kws={'label': 'Percentage Difference from Apo Closed'}, xticklabels = label, yticklabels = helices)
#ax.add_artist(lines.Line2D([0, 20], [7, 7], color = 'black', linestyle= '--', linewidth = 4))
plt.title('Helical Interactions Compared to Apo Closed WPD Loop')
plt.savefig('Hel_inter_cmpr.png')
plt.close()

