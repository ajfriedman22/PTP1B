#Import Necessary Packages
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from itertools import product
import seaborn as sns
#Load data, determine correlated samples and caculate mean and error
def load_data(folder):
    a3_a7_pt1, a3_a7_pt2, a6_a7_pt1, a6_a7_pt2, a6_a7_pt3 = [],[],[],[],[]
    for i in open('../../../' + folder + '/analysis/a3_a7_pt1_tot_inter.txt'):
        a3_a7_pt1.append(float(i))
    for i in open('../../../' + folder + '/analysis/a3_a7_pt2_tot_inter.txt'):
        a3_a7_pt2.append(float(i))
    for i in open('../../../' + folder + '/analysis/a6_a7_pt1_tot_inter.txt'):
        a6_a7_pt1.append(float(i))
    for i in open('../../../' + folder + '/analysis/a6_a7_pt2_tot_inter.txt'):
        a6_a7_pt2.append(float(i))
    for i in open('../../../' + folder + '/analysis/a6_a7_pt3_tot_inter.txt'):
        a6_a7_pt3.append(float(i))


    tot_data = [a3_a7_pt1, a3_a7_pt2, a6_a7_pt1, a6_a7_pt2, a6_a7_pt3]
    
    #Determine correlation time
    corr_time = []
    for i in range(5):
        ref = tot_data[i][0]
        t_ref = 0
        threshold = 2
        for j in range(len(tot_data[i][:])):
            if abs(tot_data[i][j] - ref) > threshold:
                corr_time.append(j - t_ref)
                ref = tot_data[i][j]
                t_ref = j
    mean_corr_time = int(np.mean(corr_time))
    num_ucorr = int(np.round(len(a3_a7_pt1)/mean_corr_time) + 1)
    
    #Seperate Uncorrelated data
    uncorr_data = np.zeros((5, num_ucorr))
    n = 0
    for i in range(len(tot_data[0][:])):
        if i % mean_corr_time == 0:
            uncorr_data[0,n] = tot_data[0][i]
            uncorr_data[1,n] = tot_data[1][i]
            uncorr_data[2,n] = tot_data[2][i]
            uncorr_data[3,n] = tot_data[3][i]
            uncorr_data[4,n] = tot_data[4][i]
            n += 1
    
    #Seperate arrays for different measurements
    a3_a7_pt1 = uncorr_data[0,:]
    a3_a7_pt2 = uncorr_data[1,:]
    a6_a7_pt1 = uncorr_data[2,:]
    a6_a7_pt2 = uncorr_data[3,:]
    a6_a7_pt3 = uncorr_data[4,:]

    return a3_a7_pt1, a3_a7_pt2, a6_a7_pt1, a6_a7_pt2, a6_a7_pt3

def plot_all(inter, err, label_inter, P, label_x, j):
    num = np.linspace(0, len(label_x), num = len(label_x))
    inter_j = inter[:][j]
    err_j = err[:][j]
    p1 = P[0][j]
    p2 = P[1][j]
    Color = ['gray', 'black', 'blue', 'red']
    #Plot Function
    fig = plt.figure(figsize=(12, 8))
    ax1 = fig.add_subplot(111)
    ax1.set_title('Comparison of ' + label_inter[j] + ' Helix Interactions')
    ax1.set_ylabel('Mean Number of Interactions')
    ax1.bar(num, inter_j, color = Color, width=0.9)
    
    plt.errorbar(num, inter_j, yerr= err_j, fmt='o', color='black')
    plt.xticks(num, label_x, fontsize=8)
    fig.savefig(label_inter[j] + '_cmpr_inter.png')
    plt.close(fig)
    if p < 0.01 and p > 0.001:
        x1, x2 = 5, 15 #Columns for Apo and AD
        y, h, col = (1.1*inter_j[[0, 2]].max()), 1, 'b'
        plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
        plt.text((x1+x2)*0.5, y+h, "**" , ha='center', va='bottom', color=col)
    if p < 0.001:
        x1, x2 = 5, 15 #Columns for Apo and AD
        y, h, col = (1.1*inter_j[[0, 2]].max()), 1, 'b'
        plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
        plt.text((x1+x2)*0.5, y+h, "***" , ha='center', va='bottom', color=col)
    if p1 < 0.05 and p1 > 0.01:
        x1, x2 = 5, 20 #Columns for Apo and BBR
        y, h, col = (1.1*inter_j[[0, 3]].max()), 1, 'r'
        plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
        plt.text((x1+x2)*0.5, y+h, "*" , ha='center', va='bottom', color=col)
    if p1 < 0.01 and p1 > 0.001:
        x1, x2 = 5, 20 #Columns for Apo and BBR
        y, h, col = (1.1*inter_j[[0, 3]].max()), 1, 'r'
        plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
        plt.text((x1+x2)*0.5, y+h, "**" , ha='center', va='bottom', color=col)
    if p1 < 0.001:
        x1, x2 = 5, 20 #Columns for Apo and BBR
        y, h, col = (1.1*inter_j[[0, 3]].max()), 1, 'r'
        plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
        plt.text((x1+x2)*0.5, y+h, "***" , ha='center', va='bottom', color=col)
    if p2 < 0.05 and p2 > 0.01:
        x1, x2 = 5, 10 #Columns for Apo and AD
        y, h, col = (1.1*inter_j[[0, 1]].max()), 1, 'k'
        plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
        plt.text((x1+x2)*0.5, y+h, "*" , ha='center', va='bottom', color=col)
    if p2 < 0.01 and p2 > 0.001:
        x1, x2 = 5, 10 #Columns for Apo and AD
        y, h, col = (1.1*inter_j[[0, 1]].max()), 1, 'k'
        plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
        plt.text((x1+x2)*0.5, y+h, "**" , ha='center', va='bottom', color=col)
    if p2 < 0.001:
        x1, x2 = 5, 10 #Columns for Apo and AD
        y, h, col = (1.1*inter_j[[0, 1]].max()), 1, 'k'
        plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
        plt.text((x1+x2)*0.5, y+h, "***" , ha='center', va='bottom', color=col)
    #Write to file
    output = open(label[j] + '_cmpr_inter.txt', 'w')
    for k in range(len(label_x)):
        output.write(str(label_x[k]) + ': ' + str(inter_j[k]) + ' +/- ' + str(err_j[k]) + '\n')
        output.write('P-value:' + str(P[k,j]) + '\n')
    output.close()

#Set Array of file folders
Folders_Apo_open = ['rebuild_a7/analysis', 'rebuild_a7_high/config9/analysis', 'rebuild_a7_high/config11/analysis']
Folders_Apo_closed = ['Apo_1SUG/analysis/1sug', 'Apo_1SUG/analysis/1sug2', 'Apo_1SUG/analysis/1sug3']
Folders_AD = ['mutate/WT/AD/analysis', '1sug_dis_AD/analysis/config11', '1sug_dis_AD/analysis/config_alt', '1sug_dis_AD/analysis/config_alt2']
Folders_BBR = ['mutate/WT/BBR/analysis', 'BBR_a7/analysis']

#Load Data for all Folders
ApoO_a3_a7_pt1, ApoO_a3_a7_pt2, ApoO_a6_a7_pt1, ApoO_a6_a7_pt2, ApoO_a6_a7_pt3 = [],[],[],[],[]
for i in Folders_Apo_open:
    a3_a7_pt1, a3_a7_pt2, a6_a7_pt1, a6_a7_pt2, a6_a7_pt3 = load_data(i)
    ApoO_a3_a7_pt1.extend(a3_a7_pt1)
    ApoO_a3_a7_pt2.extend(a3_a7_pt2)
    ApoO_a6_a7_pt1.extend(a6_a7_pt1)
    ApoO_a6_a7_pt2.extend(a6_a7_pt2)
    ApoO_a6_a7_pt3.extend(a6_a7_pt3)

ApoC_a3_a7_pt1, ApoC_a3_a7_pt2, ApoC_a6_a7_pt1, ApoC_a6_a7_pt2, ApoC_a6_a7_pt3 = [],[],[],[],[]
for i in Folders_Apo_open:
    a3_a7_pt1, a3_a7_pt2, a6_a7_pt1, a6_a7_pt2, a6_a7_pt3 = load_data(i)
    ApoC_a3_a7_pt1.extend(a3_a7_pt1)
    ApoC_a3_a7_pt2.extend(a3_a7_pt2)
    ApoC_a6_a7_pt1.extend(a6_a7_pt1)
    ApoC_a6_a7_pt2.extend(a6_a7_pt2)
    ApoC_a6_a7_pt3.extend(a6_a7_pt3)

AD_a3_a7_pt1, AD_a3_a7_pt2, AD_a6_a7_pt1, AD_a6_a7_pt2, AD_a6_a7_pt3 = [],[],[],[],[]
for i in Folders_Apo_open:
    a3_a7_pt1, a3_a7_pt2, a6_a7_pt1, a6_a7_pt2, a6_a7_pt3 = load_data(i)
    AD_a3_a7_pt1.extend(a3_a7_pt1)
    AD_a3_a7_pt2.extend(a3_a7_pt2)
    AD_a6_a7_pt1.extend(a6_a7_pt1)
    AD_a6_a7_pt2.extend(a6_a7_pt2)
    AD_a6_a7_pt3.extend(a6_a7_pt3)

BBR_a3_a7_pt1, BBR_a3_a7_pt2, BBR_a6_a7_pt1, BBR_a6_a7_pt2, BBR_a6_a7_pt3 = [],[],[],[],[]
for i in Folders_Apo_open:
    a3_a7_pt1, a3_a7_pt2, a6_a7_pt1, a6_a7_pt2, a6_a7_pt3 = load_data(i)
    BBR_a3_a7_pt1.extend(a3_a7_pt1)
    BBR_a3_a7_pt2.extend(a3_a7_pt2)
    BBR_a6_a7_pt1.extend(a6_a7_pt1)
    BBR_a6_a7_pt2.extend(a6_a7_pt2)
    BBR_a6_a7_pt3.extend(a6_a7_pt3)

#Caculate P-value relative to Apo Open
P_rel_ApoO = np.zeros((2, 5))
P_rel_ApoO[0][0] = stats.ttest_ind(ADa3_a7_pt1, ApoOa3_a7_pt1, equal_var = False) #Welch's t-test
P_rel_ApoO[0][1] = stats.ttest_ind(ADa3_a7_pt2, ApoOa3_a7_pt2, equal_var = False) #Welch's t-test
P_rel_ApoO[0][2] = stats.ttest_ind(ADa6_a7_pt1, ApoOa6_a7_pt1, equal_var = False) #Welch's t-test
P_rel_ApoO[0][3] = stats.ttest_ind(ADa6_a7_pt2, ApoOa6_a7_pt2, equal_var = False) #Welch's t-test
P_rel_ApoO[0][4] = stats.ttest_ind(ADa6_a7_pt3, ApoOa6_a7_pt3, equal_var = False) #Welch's t-test

P_rel_ApoO[1][0] = stats.ttest_ind(BBRa3_a7_pt1, ApoOa3_a7_pt1, equal_var = False) #Welch's t-test
P_rel_ApoO[1][1] = stats.ttest_ind(BBRa3_a7_pt2, ApoOa3_a7_pt2, equal_var = False) #Welch's t-test
P_rel_ApoO[1][2] = stats.ttest_ind(BBRa6_a7_pt1, ApoOa6_a7_pt1, equal_var = False) #Welch's t-test
P_rel_ApoO[1][3] = stats.ttest_ind(BBRa6_a7_pt2, ApoOa6_a7_pt2, equal_var = False) #Welch's t-test
P_rel_ApoO[1][4] = stats.ttest_ind(BBRa6_a7_pt3, ApoOa6_a7_pt3, equal_var = False) #Welch's t-test

#Mean and Standard Error into one vector
label_all = ['Apo Open', 'Apo Closed',  'AD', 'BBR']
inter = ['a3-a7-pt1', 'a3-a7-pt2', 'a6-a7-pt1', 'a6-a7-pt2', 'a6-a7-pt3']
all_mean = np.zeros((len(label_all), len(inter)))
all_err = np.zeros((len(label_all), len(inter)))

all_mean[0][0] = np.mean(ApoO_a3_a7_pt1)
all_mean[0][1] = np.mean(ApoO_a3_a7_pt2)
all_mean[0][2] = np.mean(ApoO_a6_a7_pt1)
all_mean[0][3] = np.mean(ApoO_a6_a7_pt2)
all_mean[0][4] = np.mean(ApoO_a6_a7_pt3)

all_mean[1][0] = np.mean(ApoC_a3_a7_pt1)
all_mean[1][1] = np.mean(ApoC_a3_a7_pt2)
all_mean[1][2] = np.mean(ApoC_a6_a7_pt1)
all_mean[1][3] = np.mean(ApoC_a6_a7_pt2)
all_mean[1][4] = np.mean(ApoC_a6_a7_pt3)

all_mean[2][0] = np.mean(AD_a3_a7_pt1)
all_mean[2][1] = np.mean(AD_a3_a7_pt2)
all_mean[2][2] = np.mean(AD_a6_a7_pt1)
all_mean[2][3] = np.mean(AD_a6_a7_pt2)
all_mean[2][4] = np.mean(AD_a6_a7_pt3)

all_mean[3][0] = np.mean(BBR_a3_a7_pt1)
all_mean[3][1] = np.mean(BBR_a3_a7_pt2)
all_mean[3][2] = np.mean(BBR_a6_a7_pt1)
all_mean[3][3] = np.mean(BBR_a6_a7_pt2)
all_mean[3][4] = np.mean(BBR_a6_a7_pt3)

all_err[0][0] = stats.sem(ApoO_a3_a7_pt1)
all_err[0][1] = stats.sem(ApoO_a3_a7_pt2)
all_err[0][2] = stats.sem(ApoO_a6_a7_pt1)
all_err[0][3] = stats.sem(ApoO_a6_a7_pt2)
all_err[0][4] = stats.sem(ApoO_a6_a7_pt3)

all_err[1][0] = stats.sem(ApoC_a3_a7_pt1)
all_err[1][1] = stats.sem(ApoC_a3_a7_pt2)
all_err[1][2] = stats.sem(ApoC_a6_a7_pt1)
all_err[1][3] = stats.sem(ApoC_a6_a7_pt2)
all_err[1][4] = stats.sem(ApoC_a6_a7_pt3)

all_err[2][0] = stats.sem(AD_a3_a7_pt1)
all_err[2][1] = stats.sem(AD_a3_a7_pt2)
all_err[2][2] = stats.sem(AD_a6_a7_pt1)
all_err[2][3] = stats.sem(AD_a6_a7_pt2)
all_err[2][4] = stats.sem(AD_a6_a7_pt3)

all_err[3][0] = stats.sem(BBR_a3_a7_pt1)
all_err[3][1] = stats.sem(BBR_a3_a7_pt2)
all_err[3][2] = stats.sem(BBR_a6_a7_pt1)
all_err[3][3] = stats.sem(BBR_a6_a7_pt2)
all_err[3][4] = stats.sem(BBR_a6_a7_pt3)

#Plot individual graphs
for i in range(len(inter)):
    plot_all(all_mean, all_err, inter, P_rel_ApoO, label_all, i):

#Determine percent difference from Apo Closed
label = ['Apo Open', 'AD', 'BBR']
per_diff = np.zeros((len(inter), len(label)))
index = [0, 2, 3]
for j in range(len(inter)):
    for i in range(len(label)):
        n = index[i]
        per_diff[j][i] = (all_mean[n][j] - all_mean[1][j])/((all_mean[n][j] + all_mean[1][j])/2) * 100
   
#Plot table comparing residue interactions to WT
ax = plt.figure(figsize=(12, 6), frameon=False) # no visible frame
ax = sns.heatmap(per_diff, annot=False, cmap = 'bwr', cbar = True, cbar_kws={'label': 'Percentage Difference from WT'}, vmin = -200, vmax = 200, xticklabels = mut_only, yticklabels = inter)
ax.add_artist(lines.Line2D([0, 20], [7, 7], color = 'black', linestyle= '--', linewidth = 4))
plt.title('Helical Distance Compared to WT for AD')
plt.savefig('mutate_AD_helix_dist.png')
plt.close()

