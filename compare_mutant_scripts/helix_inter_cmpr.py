#Import Necessary Packages
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as lines
from scipy import stats
from itertools import product
import seaborn as sns

def plot_func(inter, err, hel1, hel2, mut, p, p1, j, k):    
    num = [5, 10, 15, 20]
    Method = ['WT AD', 'WT BBR', mut + ' AD', mut + ' BBR']
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.set_title('Comparison of' + hel1 + ' ' + hel2 + 'Helix Interactions ' + mut)    
    ax1.set_ylabel('Mean Number of Interactions')
    ax1.bar(num, inter[[0, 1, j, k]], color = ['darkblue', 'darkred', 'blue', 'red'], width=4.5)
    plt.errorbar(num, inter[[0, 1, j, k]], yerr= err[[0, 1, j, k]], fmt='o', color='black')
    plt.xticks(num, Method, fontsize=8)
    if p < 0.05 and p > 0.01:
        x1, x2 = 5, 15 #Columns for WT and mutant AD
        y, h, col = (1.1*inter[[0,2]].max()) + 2, 2, 'b'
        plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
        plt.text((x1+x2)*0.5, y+h, "*" , ha='center', va='bottom', color=col)
    if p < 0.01 and p > 0.001:
        x1, x2 = 5, 15 #Columns for WT and mutant AD
        y, h, col = (1.1*inter[[0,2]].max()) + 2, 2, 'b'
        plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
        plt.text((x1+x2)*0.5, y+h, "**" , ha='center', va='bottom', color=col)
    if p < 0.001:
        x1, x2 = 5, 15 #Columns for WT and mutant AD
        y, h, col = (1.1*inter[[0,2]].max()) + 2, 2, 'b'
        plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
        plt.text((x1+x2)*0.5, y+h, "***" , ha='center', va='bottom', color=col)
    if p1 < 0.05 and p1 > 0.01:
        x1, x2 = 10, 20 #Columns for WT and mutant BBR
        y, h, col = (1.1*inter[[1,3]].max()) + 2, 2, 'r'
        plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
        plt.text((x1+x2)*0.5, y+h, "*" , ha='center', va='bottom', color=col)
    if p1 < 0.01 and p1 > 0.001:
        x1, x2 = 10, 20 #Columns for WT and mutant BBR
        y, h, col = (1.1*inter[[1,3]].max()) + 2, 2, 'r'
        plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
        plt.text((x1+x2)*0.5, y+h, "**" , ha='center', va='bottom', color=col)
    if p1 < 0.001:
        x1, x2 = 10, 20 #Columns for WT and mutant BBR
        y, h, col = (1.1*inter[[1,3]].max()) + 2, 2, 'r'
        plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
        plt.text((x1+x2)*0.5, y+h, "***" , ha='center', va='bottom', color=col)
    fig.savefig(mut + '/' + hel1 + '_' + hel2 + '_' + mut + '_inter.png')
    plt.close(fig)

def mut_inter(hel1, hel2, mut, res):
    mut_array = np.zeros([len(res), 2])
    num = 0
    for i in open('../../../' + mut + '/AD/analysis/' + hel1 + '_' + hel2 + '_inter_all.txt').readlines():
        mut_array[num][0] = 100 - float(i)
        num += 1
    num = 0
    for i in open('../../../' + mut + '/BBR/analysis/' + hel1 + '_' + hel2 + '_inter_all.txt').readlines():
        mut_array[num][1] = 100 - float(i)
        num += 1
    return mut_array

def corr_inter(WT_array, mut_array, mut, i, res):
    #Print % time interaction formed for those correlated to WT ligand binding
    file_mut = open('./' + mut + '/inter_cmpr.txt', 'a')
    file_mut.write(str(res[i]) + '\n')
    file_mut.write('WT AD: ' + str(WT_array[i][0]) + '\n')
    file_mut.write('WT BBR: ' + str(WT_array[i][1]) + '\n')
    file_mut.write( mut + ' AD: ' + str(mut_array[i][0]) + '\n')
    file_mut.write( mut + ' BBR: ' + str(mut_array[i][1]) + '\n')
    #Save percentage of WT value
    corr_per_AD = np.round(((mut_array[i][0]-WT_array[i][0])/((WT_array[i][0] + mut_array[i][0])/2)) * 100, decimals=0)
    corr_per_BBR = np.round(((mut_array[i][1]-WT_array[i][1])/((WT_array[i][1] + mut_array[i][1])/2)) * 100, decimals=0)
    return corr_per_AD, corr_per_BBR
   
def plot_inter(WT_array, mut_array, mut, hel1, hel2, res): 
    if abs(WT_array[i][0]-mut_array[i][0]) > 15 or abs(WT_array[i][1] - mut_array[i][1]) > 15:
        num = [1, 2, 3, 4]
        Method = ['WT AD', 'WT BBR', mut + ' AD', mut + ' BBR']
        inter_mean = np.array([WT_array[i][0], WT_array[i][1], mut_array[i][0], mut_array[i][1]])
        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        ax1.set_title('Comparison of Interactions between Residues ' + str(res[i]))    
        ax1.set_ylabel('Percent Time Interaction was Formed')
        ax1.bar(num, inter_mean, color = ['darkblue', 'darkred', 'blue', 'red'], width=4.5)
        plt.xticks(num, Method, fontsize=8)
        fig.savefig(mut + '/' + hel1 + '_' + hel2 + '_inter_' + str(res[i]) + '.png')
        plt.close(fig)

def plot_mult_cmpr(Label, index, lig, file_label, hel, hel_err, p):
    #Set num and inter
    num, inter_mean,inter_mean_err = [],[],[]
    for i in range(len(index)):
        num.append(i)
        n = index[i]
        inter_mean.append(hel[n])
        inter_mean_err.append(hel_err[n])

    #Set Colors
    color_bar = ['black', 'gray', 'blue']
    for i in range(3, len(index)):
        diff = inter_mean[i] - inter_mean[2] 
        if p[i-3] < 0.05 and diff < 0:
            color_bar.append('green')
        elif p[i-3] < 0.05 and diff > 0:
            color_bar.append('red')
        else:
            color_bar.append('lightblue')

    fig = plt.figure(figsize=(14,8))
    ax1 = fig.add_subplot(111)
    ax1.set_title('Comparison of Interactions between ' + str(file_label) + ' for ' + str(lig)) 
    ax1.set_ylabel('Mean Number of Interactions')
    ax1.bar(num, inter_mean, color = color_bar, width=0.9)
    plt.errorbar(num, inter_mean, yerr= inter_mean_err, fmt='o', color='black')
    plt.xticks(num, Label, fontsize=8)
    fig.savefig(file_label + '_' + lig + '.png')
    plt.close(fig)         


#Make open arrays for time and atomic distances
a3_a6_Apo_open, a3_a6_Apo_closed = [],[]
a3_a6_WT_AD, a3_a6_WT_BBR, a3_a6_F196A_AD, a3_a6_F196A_BBR = [],[],[],[]
a3_a6_L192A_AD, a3_a6_L192A_BBR, a3_a6_L192F_AD, a3_a6_L192F_BBR, a3_a6_L192N_AD, a3_a6_L192N_BBR = [],[],[],[],[],[]
a3_a6_L195A_AD, a3_a6_L195A_BBR, a3_a6_L195F_AD, a3_a6_L195F_BBR, a3_a6_L195N_AD, a3_a6_L195N_BBR = [],[],[],[],[],[]
a3_a6_S286A_AD, a3_a6_F280Y_AD, a3_a6_E276L_AD, a3_a6_E276F_AD, a3_a6_K279M_AD, a3_a6_K279W_AD, a3_a6_V287T_AD = [],[],[],[],[],[],[]
a3_a6_S286A_BBR, a3_a6_F280Y_BBR, a3_a6_E276L_BBR, a3_a6_E276F_BBR, a3_a6_K279M_BBR, a3_a6_K279W_BBR, a3_a6_V287T_BBR = [],[],[],[],[],[],[]

a3_a7_Apo_open, a3_a7_Apo_closed = [],[]
a3_a7_WT_AD, a3_a7_WT_BBR, a3_a7_F196A_AD, a3_a7_F196A_BBR = [],[],[],[]
a3_a7_L192A_AD, a3_a7_L192A_BBR, a3_a7_L192F_AD, a3_a7_L192F_BBR, a3_a7_L192N_AD, a3_a7_L192N_BBR = [],[],[],[],[],[]
a3_a7_L195A_AD, a3_a7_L195A_BBR, a3_a7_L195F_AD, a3_a7_L195F_BBR, a3_a7_L195N_AD, a3_a7_L195N_BBR = [],[],[],[],[],[]
a3_a7_S286A_AD, a3_a7_F280Y_AD, a3_a7_E276L_AD, a3_a7_E276F_AD, a3_a7_K279M_AD, a3_a7_K279W_AD, a3_a7_V287T_AD = [],[],[],[],[],[],[]
a3_a7_S286A_BBR, a3_a7_F280Y_BBR, a3_a7_E276L_BBR, a3_a7_E276F_BBR, a3_a7_K279M_BBR, a3_a7_K279W_BBR, a3_a7_V287T_BBR = [],[],[],[],[],[],[]

a6_a7_Apo_open, a6_a7_Apo_closed = [],[]
a6_a7_WT_AD, a6_a7_WT_BBR, a6_a7_F196A_AD, a6_a7_F196A_BBR = [],[],[],[]
a6_a7_L192A_AD, a6_a7_L192A_BBR, a6_a7_L192F_AD, a6_a7_L192F_BBR, a6_a7_L192N_AD, a6_a7_L192N_BBR = [],[],[],[],[],[]
a6_a7_L195A_AD, a6_a7_L195A_BBR, a6_a7_L195F_AD, a6_a7_L195F_BBR, a6_a7_L195N_AD, a6_a7_L195N_BBR = [],[],[],[],[],[]
a6_a7_S286A_AD, a6_a7_F280Y_AD, a6_a7_E276L_AD, a6_a7_E276F_AD, a6_a7_K279M_AD, a6_a7_K279W_AD, a6_a7_V287T_AD = [],[],[],[],[],[],[]
a6_a7_S286A_BBR, a6_a7_F280Y_BBR, a6_a7_E276L_BBR, a6_a7_E276F_BBR, a6_a7_K279M_BBR, a6_a7_K279W_BBR, a6_a7_V287T_BBR = [],[],[],[],[],[],[]

#Input Data for a3 and a6 interactions
for i in open("../../../WT/AD/analysis/a3_a6_inter.txt").readlines():
    a3_a6_WT_AD.append(float(i))
for i in open("../../../WT/BBR/analysis/a3_a6_inter.txt").readlines():
    a3_a6_WT_BBR.append(float(i))
for i in open("../../../F196A/AD/analysis/a3_a6_inter.txt").readlines():
    a3_a6_F196A_AD.append(float(i))
for i in open("../../../F196A/BBR/analysis/a3_a6_inter.txt").readlines():
    a3_a6_F196A_BBR.append(float(i))
for i in open("../../../L192A/AD/analysis/a3_a6_inter.txt").readlines():
    a3_a6_L192A_AD.append(float(i))
for i in open("../../../L192A/BBR/analysis/a3_a6_inter.txt").readlines():
    a3_a6_L192A_BBR.append(float(i))
for i in open("../../../L192F/AD/analysis/a3_a6_inter.txt").readlines():
    a3_a6_L192F_AD.append(float(i))
for i in open("../../../L192F/BBR/analysis/a3_a6_inter.txt").readlines():
    a3_a6_L192F_BBR.append(float(i))
for i in open("../../../L195A/AD/analysis/a3_a6_inter.txt").readlines():
    a3_a6_L195A_AD.append(float(i))
for i in open("../../../L195A/BBR/analysis/a3_a6_inter.txt").readlines():
    a3_a6_L195A_BBR.append(float(i))
for i in open("../../../L195F/AD/analysis/a3_a6_inter.txt").readlines():
    a3_a6_L195F_AD.append(float(i))
for i in open("../../../L195F/BBR/analysis/a3_a6_inter.txt").readlines():
    a3_a6_L195F_BBR.append(float(i))
for i in open("../../../L195N/AD/analysis/a3_a6_inter.txt").readlines():
    a3_a6_L195N_AD.append(float(i))
for i in open("../../../L195N/BBR/analysis/a3_a6_inter.txt").readlines():
    a3_a6_L195N_BBR.append(float(i))
for i in open("../../../S286A/AD/analysis/a3_a6_inter.txt").readlines():
    a3_a6_S286A_AD.append(float(i))
for i in open("../../../S286A/BBR/analysis/a3_a6_inter.txt").readlines():
    a3_a6_S286A_BBR.append(float(i))
for i in open("../../../F280Y/AD/analysis/a3_a6_inter.txt").readlines():
    a3_a6_F280Y_AD.append(float(i))
for i in open("../../../F280Y/BBR/analysis/a3_a6_inter.txt").readlines():
    a3_a6_F280Y_BBR.append(float(i))
for i in open("../../../E276L/AD/analysis/a3_a6_inter.txt").readlines():
    a3_a6_E276L_AD.append(float(i))
for i in open("../../../E276L/BBR/analysis/a3_a6_inter.txt").readlines():
    a3_a6_E276L_BBR.append(float(i))
for i in open("../../../E276F/AD/analysis/a3_a6_inter.txt").readlines():
    a3_a6_E276F_AD.append(float(i))
for i in open("../../../E276F/BBR/analysis/a3_a6_inter.txt").readlines():
    a3_a6_E276F_BBR.append(float(i))
for i in open("../../../K279M/AD/analysis/a3_a6_inter.txt").readlines():
    a3_a6_K279M_AD.append(float(i))
for i in open("../../../K279M/BBR/analysis/a3_a6_inter.txt").readlines():
    a3_a6_K279M_BBR.append(float(i))
for i in open("../../../K279W/AD/analysis/a3_a6_inter.txt").readlines():
    a3_a6_K279W_AD.append(float(i))
for i in open("../../../K279W/BBR/analysis/a3_a6_inter.txt").readlines():
    a3_a6_K279W_BBR.append(float(i))
for i in open("../../../V287T/AD/analysis/a3_a6_inter.txt").readlines():
    a3_a6_V287T_AD.append(float(i))
for i in open("../../../V287T/BBR/analysis/a3_a6_inter.txt").readlines():
    a3_a6_V287T_BBR.append(float(i))
for i in open("../../../../rebuild_a7/analysis/a3_a6_inter.txt").readlines():
    a3_a6_Apo_open.append(float(i))
for i in open("../../../../rebuild_a7_high/config9/analysis/a3_a6_inter.txt").readlines():
    a3_a6_Apo_open.append(float(i))
for i in open("../../../../rebuild_a7_high/config11/analysis/a3_a6_inter.txt").readlines():
    a3_a6_Apo_open.append(float(i))
for i in open("../../../../Apo_1SUG/analysis/1sug/a3_a6_inter.txt").readlines():
    a3_a6_Apo_closed.append(float(i))
for i in open("../../../../Apo_1SUG/analysis/1sug2/a3_a6_inter.txt").readlines():
    a3_a6_Apo_closed.append(float(i))
for i in open("../../../../Apo_1SUG/analysis/1sug3/a3_a6_inter.txt").readlines():
    a3_a6_Apo_open.append(float(i))


#Input Data for a3 and a7 interactions
for i in open("../../../WT/AD/analysis/a7_a3_inter.txt").readlines():
    a3_a7_WT_AD.append(float(i))
for i in open("../../../WT/BBR/analysis/a7_a3_inter.txt").readlines():
    a3_a7_WT_BBR.append(float(i))
for i in open("../../../F196A/AD/analysis/a7_a3_inter.txt").readlines():
    a3_a7_F196A_AD.append(float(i))
for i in open("../../../F196A/BBR/analysis/a7_a3_inter.txt").readlines():
    a3_a7_F196A_BBR.append(float(i))
for i in open("../../../L192A/AD/analysis/a7_a3_inter.txt").readlines():
    a3_a7_L192A_AD.append(float(i))
for i in open("../../../L192A/BBR/analysis/a7_a3_inter.txt").readlines():
    a3_a7_L192A_BBR.append(float(i))
for i in open("../../../L192F/AD/analysis/a7_a3_inter.txt").readlines():
    a3_a7_L192F_AD.append(float(i))
for i in open("../../../L192F/BBR/analysis/a7_a3_inter.txt").readlines():
    a3_a7_L192F_BBR.append(float(i))
for i in open("../../../L195A/AD/analysis/a7_a3_inter.txt").readlines():
    a3_a7_L195A_AD.append(float(i))
for i in open("../../../L195A/BBR/analysis/a7_a3_inter.txt").readlines():
    a3_a7_L195A_BBR.append(float(i))
for i in open("../../../L195F/AD/analysis/a7_a3_inter.txt").readlines():
    a3_a7_L195F_AD.append(float(i))
for i in open("../../../L195F/BBR/analysis/a7_a3_inter.txt").readlines():
    a3_a7_L195F_BBR.append(float(i))
for i in open("../../../L195N/AD/analysis/a7_a3_inter.txt").readlines():
    a3_a7_L195N_AD.append(float(i))
for i in open("../../../L195N/BBR/analysis/a7_a3_inter.txt").readlines():
    a3_a7_L195N_BBR.append(float(i))
for i in open("../../../S286A/AD/analysis/a7_a3_inter.txt").readlines():
    a3_a7_S286A_AD.append(float(i))
for i in open("../../../S286A/BBR/analysis/a7_a3_inter.txt").readlines():
    a3_a7_S286A_BBR.append(float(i))
for i in open("../../../F280Y/AD/analysis/a7_a3_inter.txt").readlines():
    a3_a7_F280Y_AD.append(float(i))
for i in open("../../../F280Y/BBR/analysis/a7_a3_inter.txt").readlines():
    a3_a7_F280Y_BBR.append(float(i))
for i in open("../../../E276L/AD/analysis/a7_a3_inter.txt").readlines():
    a3_a7_E276L_AD.append(float(i))
for i in open("../../../E276L/BBR/analysis/a7_a3_inter.txt").readlines():
    a3_a7_E276L_BBR.append(float(i))
for i in open("../../../E276F/AD/analysis/a7_a3_inter.txt").readlines():
    a3_a7_E276F_AD.append(float(i))
for i in open("../../../E276F/BBR/analysis/a7_a3_inter.txt").readlines():
    a3_a7_E276F_BBR.append(float(i))
for i in open("../../../K279M/AD/analysis/a7_a3_inter.txt").readlines():
    a3_a7_K279M_AD.append(float(i))
for i in open("../../../K279M/BBR/analysis/a7_a3_inter.txt").readlines():
    a3_a7_K279M_BBR.append(float(i))
for i in open("../../../K279W/AD/analysis/a7_a3_inter.txt").readlines():
    a3_a7_K279W_AD.append(float(i))
for i in open("../../../K279W/BBR/analysis/a7_a3_inter.txt").readlines():
    a3_a7_K279W_BBR.append(float(i))
for i in open("../../../V287T/AD/analysis/a7_a3_inter.txt").readlines():
    a3_a7_V287T_AD.append(float(i))
for i in open("../../../V287T/BBR/analysis/a7_a3_inter.txt").readlines():
    a3_a7_V287T_BBR.append(float(i))
for i in open("../../../../rebuild_a7/analysis/a7_a3_inter.txt").readlines():
    a3_a7_Apo_open.append(float(i))
for i in open("../../../../rebuild_a7_high/config9/analysis/a7_a3_inter.txt").readlines():
    a3_a7_Apo_open.append(float(i))
for i in open("../../../../rebuild_a7_high/config11/analysis/a7_a3_inter.txt").readlines():
    a3_a7_Apo_open.append(float(i))
for i in open("../../../../Apo_1SUG/analysis/1sug/a7_a3_inter.txt").readlines():
    a3_a7_Apo_closed.append(float(i))
for i in open("../../../../Apo_1SUG/analysis/1sug2/a7_a3_inter.txt").readlines():
    a3_a7_Apo_closed.append(float(i))
for i in open("../../../../Apo_1SUG/analysis/1sug3/a7_a3_inter.txt").readlines():
    a3_a7_Apo_open.append(float(i))

#Input Data for a6 and a7 interactions
for i in open("../../../WT/AD/analysis/a7_a6_inter.txt").readlines():
    a6_a7_WT_AD.append(float(i))
for i in open("../../../WT/BBR/analysis/a7_a6_inter.txt").readlines():
    a6_a7_WT_BBR.append(float(i))
for i in open("../../../F196A/AD/analysis/a7_a6_inter.txt").readlines():
    a6_a7_F196A_AD.append(float(i))
for i in open("../../../F196A/BBR/analysis/a7_a6_inter.txt").readlines():
    a6_a7_F196A_BBR.append(float(i))
for i in open("../../../L192A/AD/analysis/a7_a6_inter.txt").readlines():
    a6_a7_L192A_AD.append(float(i))
for i in open("../../../L192A/BBR/analysis/a7_a6_inter.txt").readlines():
    a6_a7_L192A_BBR.append(float(i))
for i in open("../../../L192F/AD/analysis/a7_a6_inter.txt").readlines():
    a6_a7_L192F_AD.append(float(i))
for i in open("../../../L192F/BBR/analysis/a7_a6_inter.txt").readlines():
    a6_a7_L192F_BBR.append(float(i))
for i in open("../../../L195A/AD/analysis/a7_a6_inter.txt").readlines():
    a6_a7_L195A_AD.append(float(i))
for i in open("../../../L195A/BBR/analysis/a7_a6_inter.txt").readlines():
    a6_a7_L195A_BBR.append(float(i))
for i in open("../../../L195F/AD/analysis/a7_a6_inter.txt").readlines():
    a6_a7_L195F_AD.append(float(i))
for i in open("../../../L195F/BBR/analysis/a7_a6_inter.txt").readlines():
    a6_a7_L195F_BBR.append(float(i))
for i in open("../../../L195N/AD/analysis/a7_a6_inter.txt").readlines():
    a6_a7_L195N_AD.append(float(i))
for i in open("../../../L195N/BBR/analysis/a7_a6_inter.txt").readlines():
    a6_a7_L195N_BBR.append(float(i))
for i in open("../../../S286A/AD/analysis/a7_a6_inter.txt").readlines():
    a6_a7_S286A_AD.append(float(i))
for i in open("../../../S286A/BBR/analysis/a7_a6_inter.txt").readlines():
    a6_a7_S286A_BBR.append(float(i))
for i in open("../../../F280Y/AD/analysis/a7_a6_inter.txt").readlines():
    a6_a7_F280Y_AD.append(float(i))
for i in open("../../../F280Y/BBR/analysis/a7_a6_inter.txt").readlines():
    a6_a7_F280Y_BBR.append(float(i))
for i in open("../../../E276L/AD/analysis/a7_a6_inter.txt").readlines():
    a6_a7_E276L_AD.append(float(i))
for i in open("../../../E276L/BBR/analysis/a7_a6_inter.txt").readlines():
    a6_a7_E276L_BBR.append(float(i))
for i in open("../../../E276F/AD/analysis/a7_a6_inter.txt").readlines():
    a6_a7_E276F_AD.append(float(i))
for i in open("../../../E276F/BBR/analysis/a7_a6_inter.txt").readlines():
    a6_a7_E276F_BBR.append(float(i))
for i in open("../../../K279M/AD/analysis/a7_a6_inter.txt").readlines():
    a6_a7_K279M_AD.append(float(i))
for i in open("../../../K279M/BBR/analysis/a7_a6_inter.txt").readlines():
    a6_a7_K279M_BBR.append(float(i))
for i in open("../../../K279W/AD/analysis/a7_a6_inter.txt").readlines():
    a6_a7_K279W_AD.append(float(i))
for i in open("../../../K279W/BBR/analysis/a7_a6_inter.txt").readlines():
    a6_a7_K279W_BBR.append(float(i))
for i in open("../../../V287T/AD/analysis/a7_a6_inter.txt").readlines():
    a6_a7_V287T_AD.append(float(i))
for i in open("../../../V287T/BBR/analysis/a7_a6_inter.txt").readlines():
    a6_a7_V287T_BBR.append(float(i))
for i in open("../../../../rebuild_a7/analysis/a7_a6_inter.txt").readlines():
    a6_a7_Apo_open.append(float(i))
for i in open("../../../../rebuild_a7_high/config9/analysis/a7_a6_inter.txt").readlines():
    a6_a7_Apo_open.append(float(i))
for i in open("../../../../rebuild_a7_high/config11/analysis/a7_a6_inter.txt").readlines():
    a6_a7_Apo_open.append(float(i))
for i in open("../../../../Apo_1SUG/analysis/1sug/a7_a6_inter.txt").readlines():
    a6_a7_Apo_closed.append(float(i))
for i in open("../../../../Apo_1SUG/analysis/1sug2/a7_a6_inter.txt").readlines():
    a6_a7_Apo_closed.append(float(i))
for i in open("../../../../Apo_1SUG/analysis/1sug3/a7_a6_inter.txt").readlines():
    a6_a7_Apo_open.append(float(i))

#Calculate mean and sem for interactions
a3_a6 = np.array([np.mean(a3_a6_WT_AD), np.mean(a3_a6_WT_BBR), np.mean(a3_a6_F196A_AD), np.mean(a3_a6_F196A_BBR), np.mean(a3_a6_L192A_AD), 
    np.mean(a3_a6_L192A_BBR), np.mean(a3_a6_L192F_AD), np.mean(a3_a6_L192F_BBR), np.mean(a3_a6_L192N_AD), np.mean(a3_a6_L192N_BBR), 
    np.mean(a3_a6_L195A_AD), np.mean(a3_a6_L195A_BBR), np.mean(a3_a6_L195F_AD), np.mean(a3_a6_L195F_BBR), np.mean(a3_a6_L195N_AD), 
    np.mean(a3_a6_L195N_BBR),
    np.mean(a3_a6_S286A_AD),  np.mean(a3_a6_S286A_BBR), np.mean(a3_a6_F280Y_AD), 
    np.mean(a3_a6_F280Y_BBR), np.mean(a3_a6_E276L_AD), np.mean(a3_a6_E276L_BBR),
    np.mean(a3_a6_E276F_AD), np.mean(a3_a6_E276F_BBR), np.mean(a3_a6_K279M_AD),
    np.mean(a3_a6_K279M_BBR), np.mean(a3_a6_K279W_AD), np.mean(a3_a6_K279W_BBR),
    np.mean(a3_a6_V287T_AD), np.mean(a3_a6_V287T_BBR), np.mean(a3_a6_Apo_open), np.mean(a3_a6_Apo_closed)])

a3_a7 = np.array([np.mean(a3_a7_WT_AD), np.mean(a3_a7_WT_BBR), np.mean(a3_a7_F196A_AD), np.mean(a3_a7_F196A_BBR), np.mean(a3_a7_L192A_AD), 
    np.mean(a3_a7_L192A_BBR), np.mean(a3_a7_L192F_AD), np.mean(a3_a7_L192F_BBR), np.mean(a3_a7_L192N_AD), np.mean(a3_a7_L192N_BBR), 
    np.mean(a3_a7_L195A_AD), np.mean(a3_a7_L195A_BBR), np.mean(a3_a7_L195F_AD), np.mean(a3_a7_L195F_BBR), np.mean(a3_a7_L195N_AD), 
    np.mean(a3_a7_L195N_BBR),
    np.mean(a3_a7_S286A_AD), np.mean(a3_a7_S286A_BBR), np.mean(a3_a7_F280Y_AD), 
    np.mean(a3_a7_F280Y_BBR), np.mean(a3_a7_E276L_AD), np.mean(a3_a7_E276L_BBR),
    np.mean(a3_a7_E276F_AD), np.mean(a3_a7_E276F_BBR), np.mean(a3_a7_K279M_AD),
    np.mean(a3_a7_K279M_BBR), np.mean(a3_a7_K279W_AD), np.mean(a3_a7_K279W_BBR),
    np.mean(a3_a7_V287T_AD), np.mean(a3_a7_V287T_BBR), np.mean(a3_a7_Apo_open), np.mean(a3_a7_Apo_closed)])

a6_a7 = np.array([np.mean(a6_a7_WT_AD), np.mean(a6_a7_WT_BBR), np.mean(a6_a7_F196A_AD), np.mean(a6_a7_F196A_BBR), np.mean(a6_a7_L192A_AD), 
    np.mean(a6_a7_L192A_BBR), np.mean(a6_a7_L192F_AD), np.mean(a6_a7_L192F_BBR), np.mean(a6_a7_L192N_AD), np.mean(a6_a7_L192N_BBR), 
    np.mean(a6_a7_L195A_AD), np.mean(a6_a7_L195A_BBR), np.mean(a6_a7_L195F_AD), np.mean(a6_a7_L195F_BBR), np.mean(a6_a7_L195N_AD), 
    np.mean(a6_a7_L195N_BBR),
    np.mean(a6_a7_S286A_AD),  np.mean(a6_a7_S286A_BBR), np.mean(a6_a7_F280Y_AD), 
    np.mean(a6_a7_F280Y_BBR), np.mean(a6_a7_E276L_AD), np.mean(a6_a7_E276L_BBR),
    np.mean(a6_a7_E276F_AD), np.mean(a6_a7_E276F_BBR), np.mean(a6_a7_K279M_AD),
    np.mean(a6_a7_K279M_BBR), np.mean(a6_a7_K279W_AD), np.mean(a6_a7_K279W_BBR),
    np.mean(a6_a7_V287T_AD), np.mean(a6_a7_V287T_BBR), np.mean(a6_a7_Apo_open), np.mean(a6_a7_Apo_closed)])

a3_a6_err = np.array([stats.sem(a3_a6_WT_AD), stats.sem(a3_a6_WT_BBR), stats.sem(a3_a6_F196A_AD), stats.sem(a3_a6_F196A_BBR), stats.sem(a3_a6_L192A_AD), 
    stats.sem(a3_a6_L192A_BBR), stats.sem(a3_a6_L192F_AD), stats.sem(a3_a6_L192F_BBR), stats.sem(a3_a6_L192N_AD), stats.sem(a3_a6_L192N_BBR), 
    stats.sem(a3_a6_L195A_AD), stats.sem(a3_a6_L195A_BBR), stats.sem(a3_a6_L195F_AD), stats.sem(a3_a6_L195F_BBR), stats.sem(a3_a6_L195N_AD), 
    stats.sem(a3_a6_L195N_BBR),
    stats.sem(a3_a6_S286A_AD),  stats.sem(a3_a6_S286A_BBR), stats.sem(a3_a6_F280Y_AD), 
    stats.sem(a3_a6_F280Y_BBR), stats.sem(a3_a6_E276L_AD), stats.sem(a3_a6_E276L_BBR),
    stats.sem(a3_a6_E276F_AD), stats.sem(a3_a6_E276F_BBR), stats.sem(a3_a6_K279M_AD),
    stats.sem(a3_a6_K279M_BBR), stats.sem(a3_a6_K279W_AD), stats.sem(a3_a6_K279W_BBR),
    stats.sem(a3_a6_V287T_AD), stats.sem(a3_a6_V287T_BBR), stats.sem(a3_a6_Apo_open), stats.sem(a3_a6_Apo_closed)])

a3_a7_err = np.array([stats.sem(a3_a7_WT_AD), stats.sem(a3_a7_WT_BBR), stats.sem(a3_a7_F196A_AD), stats.sem(a3_a7_F196A_BBR), stats.sem(a3_a7_L192A_AD), 
    stats.sem(a3_a7_L192A_BBR), stats.sem(a3_a7_L192F_AD), stats.sem(a3_a7_L192F_BBR), stats.sem(a3_a7_L192N_AD), stats.sem(a3_a7_L192N_BBR), 
    stats.sem(a3_a7_L195A_AD), stats.sem(a3_a7_L195A_BBR), stats.sem(a3_a7_L195F_AD), stats.sem(a3_a7_L195F_BBR), stats.sem(a3_a7_L195N_AD), 
    stats.sem(a3_a7_L195N_BBR),
    stats.sem(a3_a7_S286A_AD),  stats.sem(a3_a7_S286A_BBR), stats.sem(a3_a7_F280Y_AD), 
    stats.sem(a3_a7_F280Y_BBR), stats.sem(a3_a7_E276L_AD), stats.sem(a3_a7_E276L_BBR),
    stats.sem(a3_a7_E276F_AD), stats.sem(a3_a7_E276F_BBR), stats.sem(a3_a7_K279M_AD),
    stats.sem(a3_a7_K279M_BBR), stats.sem(a3_a7_K279W_AD), stats.sem(a3_a7_K279W_BBR),
    stats.sem(a3_a7_V287T_AD), stats.sem(a3_a7_V287T_BBR), stats.sem(a3_a7_Apo_open), stats.sem(a3_a7_Apo_closed)])

a6_a7_err = np.array([stats.sem(a6_a7_WT_AD), stats.sem(a6_a7_WT_BBR), stats.sem(a6_a7_F196A_AD), stats.sem(a6_a7_F196A_BBR), stats.sem(a6_a7_L192A_AD), 
    stats.sem(a6_a7_L192A_BBR), stats.sem(a6_a7_L192F_AD), stats.sem(a6_a7_L192F_BBR), stats.sem(a6_a7_L192N_AD), stats.sem(a6_a7_L192N_BBR), 
    stats.sem(a6_a7_L195A_AD), stats.sem(a6_a7_L195A_BBR), stats.sem(a6_a7_L195F_AD), stats.sem(a6_a7_L195F_BBR), stats.sem(a6_a7_L195N_AD), 
    stats.sem(a6_a7_L195N_BBR),
    stats.sem(a6_a7_S286A_AD),  stats.sem(a6_a7_S286A_BBR), stats.sem(a6_a7_F280Y_AD), 
    stats.sem(a6_a7_F280Y_BBR), stats.sem(a6_a7_E276L_AD), stats.sem(a6_a7_E276L_BBR),
    stats.sem(a6_a7_E276F_AD), stats.sem(a6_a7_E276F_BBR), stats.sem(a6_a7_K279M_AD),
    stats.sem(a6_a7_K279M_BBR), stats.sem(a6_a7_K279W_AD), stats.sem(a6_a7_K279W_BBR),
    stats.sem(a6_a7_V287T_AD), stats.sem(a6_a7_V287T_BBR), stats.sem(a6_a7_Apo_open), stats.sem(a6_a7_Apo_closed)])

#Run t-test between groups
st, p1 = stats.ttest_ind(a3_a6_WT_AD, a3_a6_F196A_AD, equal_var = False) #Welch's t-test b/w WT and F196A AD
st, p2 = stats.ttest_ind(a3_a6_WT_AD, a3_a6_L192A_AD, equal_var = False) #Welch's t-test b/w WT and L192A AD
st, p3 = stats.ttest_ind(a3_a6_WT_AD, a3_a6_L192F_AD, equal_var = False) #Welch's t-test b/w WT and L192F AD
st, p4 = stats.ttest_ind(a3_a6_WT_AD, a3_a6_L192N_AD, equal_var = False) #Welch's t-test b/w WT and L192N AD
st, p5 = stats.ttest_ind(a3_a6_WT_AD, a3_a6_L195A_AD, equal_var = False) #Welch's t-test b/w WT and L195A AD
st, p6 = stats.ttest_ind(a3_a6_WT_AD, a3_a6_L195F_AD, equal_var = False) #Welch's t-test b/w WT and L195F AD
st, p7 = stats.ttest_ind(a3_a6_WT_AD, a3_a6_L195N_AD, equal_var = False) #Welch's t-test b/w WT and L195N AD

st, p1B = stats.ttest_ind(a3_a6_WT_BBR, a3_a6_F196A_BBR, equal_var = False) #Welch's t-test b/w WT and F196A AD
st, p2B = stats.ttest_ind(a3_a6_WT_BBR, a3_a6_L192A_BBR, equal_var = False) #Welch's t-test b/w WT and L192A AD
st, p3B = stats.ttest_ind(a3_a6_WT_BBR, a3_a6_L192F_BBR, equal_var = False) #Welch's t-test b/w WT and L192F AD
st, p4B = stats.ttest_ind(a3_a6_WT_BBR, a3_a6_L192N_BBR, equal_var = False) #Welch's t-test b/w WT and L192N AD
st, p5B = stats.ttest_ind(a3_a6_WT_BBR, a3_a6_L195A_BBR, equal_var = False) #Welch's t-test b/w WT and L195A AD
st, p6B = stats.ttest_ind(a3_a6_WT_BBR, a3_a6_L195F_BBR, equal_var = False) #Welch's t-test b/w WT and L195F AD
st, p7B = stats.ttest_ind(a3_a6_WT_BBR, a3_a6_L195N_BBR, equal_var = False) #Welch's t-test b/w WT and L195N AD

st, p8 = stats.ttest_ind(a3_a7_WT_AD, a3_a7_F196A_AD, equal_var = False) #Welch's t-test b/w WT and F196A AD
st, p9 = stats.ttest_ind(a3_a7_WT_AD, a3_a7_L192A_AD, equal_var = False) #Welch's t-test b/w WT and L192A AD
st, p10 = stats.ttest_ind(a3_a7_WT_AD, a3_a7_L192F_AD, equal_var = False) #Welch's t-test b/w WT and L192F AD
st, p11 = stats.ttest_ind(a3_a7_WT_AD, a3_a7_L192N_AD, equal_var = False) #Welch's t-test b/w WT and L192N AD
st, p12 = stats.ttest_ind(a3_a7_WT_AD, a3_a7_L195A_AD, equal_var = False) #Welch's t-test b/w WT and L195A AD
st, p13 = stats.ttest_ind(a3_a7_WT_AD, a3_a7_L195F_AD, equal_var = False) #Welch's t-test b/w WT and L195F AD
st, p14 = stats.ttest_ind(a3_a7_WT_AD, a3_a7_L195N_AD, equal_var = False) #Welch's t-test b/w WT and L195N AD

st, p8B = stats.ttest_ind(a3_a7_WT_BBR, a3_a7_F196A_BBR, equal_var = False) #Welch's t-test b/w WT and F196A AD
st, p9B = stats.ttest_ind(a3_a7_WT_BBR, a3_a7_L192A_BBR, equal_var = False) #Welch's t-test b/w WT and L192A AD
st, p10B = stats.ttest_ind(a3_a7_WT_BBR, a3_a7_L192F_BBR, equal_var = False) #Welch's t-test b/w WT and L192F AD
st, p11B = stats.ttest_ind(a3_a7_WT_BBR, a3_a7_L192N_BBR, equal_var = False) #Welch's t-test b/w WT and L192N AD
st, p12B = stats.ttest_ind(a3_a7_WT_BBR, a3_a7_L195A_BBR, equal_var = False) #Welch's t-test b/w WT and L195A AD
st, p13B = stats.ttest_ind(a3_a7_WT_BBR, a3_a7_L195F_BBR, equal_var = False) #Welch's t-test b/w WT and L195F AD
st, p14B = stats.ttest_ind(a3_a7_WT_BBR, a3_a7_L195N_BBR, equal_var = False) #Welch's t-test b/w WT and L195N AD

st, p15 = stats.ttest_ind(a6_a7_WT_AD, a6_a7_F196A_AD, equal_var = False) #Welch's t-test b/w WT and F196A AD
st, p16 = stats.ttest_ind(a6_a7_WT_AD, a6_a7_L192A_AD, equal_var = False) #Welch's t-test b/w WT and L192A AD
st, p17 = stats.ttest_ind(a6_a7_WT_AD, a6_a7_L192F_AD, equal_var = False) #Welch's t-test b/w WT and L192F AD
st, p18 = stats.ttest_ind(a6_a7_WT_AD, a6_a7_L192N_AD, equal_var = False) #Welch's t-test b/w WT and L192N AD
st, p19 = stats.ttest_ind(a6_a7_WT_AD, a6_a7_L195A_AD, equal_var = False) #Welch's t-test b/w WT and L195A AD
st, p20 = stats.ttest_ind(a6_a7_WT_AD, a6_a7_L195F_AD, equal_var = False) #Welch's t-test b/w WT and L195F AD
st, p21 = stats.ttest_ind(a6_a7_WT_AD, a6_a7_L195N_AD, equal_var = False) #Welch's t-test b/w WT and L195N AD

st, p15B = stats.ttest_ind(a6_a7_WT_BBR, a6_a7_F196A_BBR, equal_var = False) #Welch's t-test b/w WT and F196A AD
st, p16B = stats.ttest_ind(a6_a7_WT_BBR, a6_a7_L192A_BBR, equal_var = False) #Welch's t-test b/w WT and L192A AD
st, p17B = stats.ttest_ind(a6_a7_WT_BBR, a6_a7_L192F_BBR, equal_var = False) #Welch's t-test b/w WT and L192F AD
st, p18B = stats.ttest_ind(a6_a7_WT_BBR, a6_a7_L192N_BBR, equal_var = False) #Welch's t-test b/w WT and L192N AD
st, p19B = stats.ttest_ind(a6_a7_WT_BBR, a6_a7_L195A_BBR, equal_var = False) #Welch's t-test b/w WT and L195A AD
st, p20B = stats.ttest_ind(a6_a7_WT_BBR, a6_a7_L195F_BBR, equal_var = False) #Welch's t-test b/w WT and L195F AD
st, p21B = stats.ttest_ind(a6_a7_WT_BBR, a6_a7_L195N_BBR, equal_var = False) #Welch's t-test b/w WT and L195N AD

st, p19 = stats.ttest_ind(a3_a6_WT_AD, a3_a6_S286A_AD, equal_var = False) #Welch's t-test b/w WT and F196A AD
st, p20 = stats.ttest_ind(a3_a6_WT_AD, a3_a6_F280Y_AD, equal_var = False) #Welch's t-test b/w WT and L192A AD
st, p21 = stats.ttest_ind(a3_a6_WT_AD, a3_a6_E276L_AD, equal_var = False) #Welch's t-test b/w WT and L192F AD
st, p22 = stats.ttest_ind(a3_a6_WT_AD, a3_a6_E276F_AD, equal_var = False) #Welch's t-test b/w WT and L195A AD
st, p23 = stats.ttest_ind(a3_a6_WT_AD, a3_a6_K279M_AD, equal_var = False) #Welch's t-test b/w WT and L195F AD
st, p24 = stats.ttest_ind(a3_a6_WT_AD, a3_a6_K279W_AD, equal_var = False) #Welch's t-test b/w WT and L195N AD
st, p25 = stats.ttest_ind(a3_a6_WT_AD, a3_a6_V287T_AD, equal_var = False) #Welch's t-test b/w WT and L195N AD

st, p26 = stats.ttest_ind(a3_a7_WT_AD, a3_a7_S286A_AD, equal_var = False) #Welch's t-test b/w WT and F196A AD
st, p27 = stats.ttest_ind(a3_a7_WT_AD, a3_a7_F280Y_AD, equal_var = False) #Welch's t-test b/w WT and L192A AD
st, p28 = stats.ttest_ind(a3_a7_WT_AD, a3_a7_E276L_AD, equal_var = False) #Welch's t-test b/w WT and L192F AD
st, p29 = stats.ttest_ind(a3_a7_WT_AD, a3_a7_E276F_AD, equal_var = False) #Welch's t-test b/w WT and L195A AD
st, p30 = stats.ttest_ind(a3_a7_WT_AD, a3_a7_K279M_AD, equal_var = False) #Welch's t-test b/w WT and L195F AD
st, p31 = stats.ttest_ind(a3_a7_WT_AD, a3_a7_K279W_AD, equal_var = False) #Welch's t-test b/w WT and L195N AD
st, p32 = stats.ttest_ind(a3_a7_WT_AD, a3_a7_V287T_AD, equal_var = False) #Welch's t-test b/w WT and L195N AD

st, p33 = stats.ttest_ind(a6_a7_WT_AD, a6_a7_S286A_AD, equal_var = False) #Welch's t-test b/w WT and F196A AD
st, p34 = stats.ttest_ind(a6_a7_WT_AD, a6_a7_F280Y_AD, equal_var = False) #Welch's t-test b/w WT and L192A AD
st, p35 = stats.ttest_ind(a6_a7_WT_AD, a6_a7_E276L_AD, equal_var = False) #Welch's t-test b/w WT and L192F AD
st, p36 = stats.ttest_ind(a6_a7_WT_AD, a6_a7_E276F_AD, equal_var = False) #Welch's t-test b/w WT and L195A AD
st, p37 = stats.ttest_ind(a6_a7_WT_AD, a6_a7_K279M_AD, equal_var = False) #Welch's t-test b/w WT and L195F AD
st, p38 = stats.ttest_ind(a6_a7_WT_AD, a6_a7_K279W_AD, equal_var = False) #Welch's t-test b/w WT and L195N AD
st, p39 = stats.ttest_ind(a6_a7_WT_AD, a6_a7_V287T_AD, equal_var = False) #Welch's t-test b/w WT and L195N AD

st, p19B = stats.ttest_ind(a3_a6_WT_BBR, a3_a6_S286A_BBR, equal_var = False) #Welch's t-test b/w WT and F196A AD
st, p20B = stats.ttest_ind(a3_a6_WT_BBR, a3_a6_F280Y_BBR, equal_var = False) #Welch's t-test b/w WT and L192A AD
st, p21B = stats.ttest_ind(a3_a6_WT_BBR, a3_a6_E276L_BBR, equal_var = False) #Welch's t-test b/w WT and L192F AD
st, p22B = stats.ttest_ind(a3_a6_WT_BBR, a3_a6_E276F_BBR, equal_var = False) #Welch's t-test b/w WT and L195A AD
st, p23B = stats.ttest_ind(a3_a6_WT_BBR, a3_a6_K279M_BBR, equal_var = False) #Welch's t-test b/w WT and L195F AD
st, p24B = stats.ttest_ind(a3_a6_WT_BBR, a3_a6_K279W_BBR, equal_var = False) #Welch's t-test b/w WT and L195N AD
st, p25B = stats.ttest_ind(a3_a6_WT_BBR, a3_a6_V287T_BBR, equal_var = False) #Welch's t-test b/w WT and L195N AD

st, p26B = stats.ttest_ind(a3_a7_WT_BBR, a3_a7_S286A_BBR, equal_var = False) #Welch's t-test b/w WT and F196A AD
st, p27B = stats.ttest_ind(a3_a7_WT_BBR, a3_a7_F280Y_BBR, equal_var = False) #Welch's t-test b/w WT and L192A AD
st, p28B = stats.ttest_ind(a3_a7_WT_BBR, a3_a7_E276L_BBR, equal_var = False) #Welch's t-test b/w WT and L192F AD
st, p29B = stats.ttest_ind(a3_a7_WT_BBR, a3_a7_E276F_BBR, equal_var = False) #Welch's t-test b/w WT and L195A AD
st, p30B = stats.ttest_ind(a3_a7_WT_BBR, a3_a7_K279M_BBR, equal_var = False) #Welch's t-test b/w WT and L195F AD
st, p31B = stats.ttest_ind(a3_a7_WT_BBR, a3_a7_K279W_BBR, equal_var = False) #Welch's t-test b/w WT and L195N AD
st, p32B = stats.ttest_ind(a3_a7_WT_BBR, a3_a7_V287T_BBR, equal_var = False) #Welch's t-test b/w WT and L195N AD

st, p33B = stats.ttest_ind(a6_a7_WT_BBR, a6_a7_S286A_BBR, equal_var = False) #Welch's t-test b/w WT and F196A AD
st, p34B = stats.ttest_ind(a6_a7_WT_BBR, a6_a7_F280Y_BBR, equal_var = False) #Welch's t-test b/w WT and L192A AD
st, p35B = stats.ttest_ind(a6_a7_WT_BBR, a6_a7_E276L_BBR, equal_var = False) #Welch's t-test b/w WT and L192F AD
st, p36B = stats.ttest_ind(a6_a7_WT_BBR, a6_a7_E276F_BBR, equal_var = False) #Welch's t-test b/w WT and L195A AD
st, p37B = stats.ttest_ind(a6_a7_WT_BBR, a6_a7_K279M_BBR, equal_var = False) #Welch's t-test b/w WT and L195F AD
st, p38B = stats.ttest_ind(a6_a7_WT_BBR, a6_a7_K279W_BBR, equal_var = False) #Welch's t-test b/w WT and L195N AD
st, p39B = stats.ttest_ind(a6_a7_WT_BBR, a6_a7_V287T_BBR, equal_var = False) #Welch's t-test b/w WT and L195N AD

plot_func(a3_a6, a3_a6_err, 'a3', 'a6', 'F196A', p1, p1B, 2, 3)
plot_func(a3_a6, a3_a6_err, 'a3', 'a6', 'L192A', p2, p2B, 4, 5)
plot_func(a3_a6, a3_a6_err, 'a3', 'a6', 'L192F', p3, p3B, 6, 7)
plot_func(a3_a6, a3_a6_err, 'a3', 'a6', 'L192N', p4, p4B, 8, 9)
plot_func(a3_a6, a3_a6_err, 'a3', 'a6', 'L195A', p5, p5B, 10, 11)
plot_func(a3_a6, a3_a6_err, 'a3', 'a6', 'L195F', p6, p6B, 12, 13)
plot_func(a3_a6, a3_a6_err, 'a3', 'a6', 'L195N', p7, p7B, 14, 15)

plot_func(a3_a7, a3_a7_err, 'a3', 'a7', 'F196A', p8, p8B, 2, 3)
plot_func(a3_a7, a3_a7_err, 'a3', 'a7', 'L192A', p9, p9B, 4, 5)
plot_func(a3_a7, a3_a7_err, 'a3', 'a7', 'L192F', p10, p10B, 6, 7)
plot_func(a3_a7, a3_a7_err, 'a3', 'a7', 'L192N', p11, p11B, 8, 9)
plot_func(a3_a7, a3_a7_err, 'a3', 'a7', 'L195A', p12, p12B, 10, 11)
plot_func(a3_a7, a3_a7_err, 'a3', 'a7', 'L195F', p13, p13B, 12, 13)
plot_func(a3_a7, a3_a7_err, 'a3', 'a7', 'L195N', p14, p14B, 14, 15)

plot_func(a6_a7, a6_a7_err, 'a6', 'a7', 'F196A', p15, p15B, 2, 3)
plot_func(a6_a7, a6_a7_err, 'a6', 'a7', 'L192A', p16, p16B, 4, 5)
plot_func(a6_a7, a6_a7_err, 'a6', 'a7', 'L192F', p17, p17B, 6, 7)
plot_func(a6_a7, a6_a7_err, 'a6', 'a7', 'L192N', p18, p18B, 8, 9)
plot_func(a6_a7, a6_a7_err, 'a6', 'a7', 'L195A', p19, p19B, 10, 11)
plot_func(a6_a7, a6_a7_err, 'a6', 'a7', 'L195F', p20, p20B, 12, 13)
plot_func(a6_a7, a6_a7_err, 'a6', 'a7', 'L195N', p21, p21B, 14, 15)

plot_func(a3_a6, a3_a6_err, 'a3', 'a6', 'S286A', p19, p19B, 16, 17)
plot_func(a3_a6, a3_a6_err, 'a3', 'a6', 'F280Y', p20, p20B, 18, 19)
plot_func(a3_a6, a3_a6_err, 'a3', 'a6', 'E276L', p21, p21B, 20, 21)
plot_func(a3_a6, a3_a6_err, 'a3', 'a6', 'E276F', p22, p22B, 22, 23)
plot_func(a3_a6, a3_a6_err, 'a3', 'a6', 'K279M', p23, p23B, 24, 25)
plot_func(a3_a6, a3_a6_err, 'a3', 'a6', 'K279W', p24, p24B, 26, 27)
plot_func(a3_a6, a3_a6_err, 'a3', 'a6', 'V287T', p25, p25B, 28, 29)

plot_func(a3_a7, a3_a7_err, 'a3', 'a7', 'S286A', p26, p26B, 16, 17)
plot_func(a3_a7, a3_a7_err, 'a3', 'a7', 'F280Y', p27, p27B, 18, 19)
plot_func(a3_a7, a3_a7_err, 'a3', 'a7', 'E276L', p28, p28B, 20, 21)
plot_func(a3_a7, a3_a7_err, 'a3', 'a7', 'E276F', p29, p29B, 22, 23)
plot_func(a3_a7, a3_a7_err, 'a3', 'a7', 'K279M', p30, p30B, 24, 25)
plot_func(a3_a7, a3_a7_err, 'a3', 'a7', 'K279W', p31, p31B, 26, 27)
plot_func(a3_a7, a3_a7_err, 'a3', 'a7', 'V287T', p32, p32B, 28, 29)

plot_func(a6_a7, a6_a7_err, 'a6', 'a7', 'S286A', p33, p33B, 16, 17)
plot_func(a6_a7, a6_a7_err, 'a6', 'a7', 'F280Y', p34, p34B, 18, 19)
plot_func(a6_a7, a6_a7_err, 'a6', 'a7', 'E276L', p35, p35B, 20, 21)
plot_func(a6_a7, a6_a7_err, 'a6', 'a7', 'E276F', p36, p36B, 22, 23)
plot_func(a6_a7, a6_a7_err, 'a6', 'a7', 'K279M', p37, p37B, 24, 25)
plot_func(a6_a7, a6_a7_err, 'a6', 'a7', 'K279W', p38, p38B, 26, 27)
plot_func(a6_a7, a6_a7_err, 'a6', 'a7', 'V287T', p39, p39B, 28, 29)

#Make array of p-values for all mutants
p_a3_a7_AD = [p8, p9, p10, p12, p13, p14, p26, p27, p28, p29, p30, p31, p32]
p_a6_a7_AD = [p15, p16, p17, p19, p20, p21, p33, p34, p35, p36, p37, p38, p39]
p_a3_a6_AD = [p1, p2, p3, p5, p6, p7, p19, p20, p21, p22, p23, p24, p25]

p_a3_a7_BBR = [p8B, p9B, p10B, p12B, p13B, p14B, p26B, p27B, p28B, p29B, p30B, p31B, p32B]
p_a6_a7_BBR = [p15B, p16B, p17B, p19B, p20B, p21B, p33B, p34B, p35B, p36B, p37B, p38B, p39B]
p_a3_a6_BBR = [p1B, p2B, p3B, p5B, p6B, p7B, p19B, p20B, p21B, p22B, p23B, p24B, p25B]

#Save All Values to Files
Label = ['Apo Open', 'Apo Closed', 'WT AD', 'WT BBR', 'F196A AD', 'F196A BBR', 'L192A AD', 'L192A BBR', 'L192F AD', 'L192F BBR', 'L195A AD', 'L195A BBR', 'L195F AD', 'L195F BBR', 'L195N AD', 'L195N BBR', 'S286A AD', 'S286A BBR', 
        'F280Y AD', 'F280Y BBR', 'E276L AD', 'E276L BBR', 'E276F AD', 'E276F BBR', 'K279M AD', 'K279M BBR', 'K279W AD', 'K279W BBR', 'V287T AD', 'V287T BBR']
index = [30, 31, 0, 1, 2, 3, 4, 5, 6, 7, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29]
output = open('a3_a7_inters.txt', 'w')
for i in range(len(Label)):
    n = int(index[i])
    output.write(Label[i] + ': ' + str(a3_a7[n]) + '+/-' + str(a3_a7_err[n]) + '\n')
output.close()

output = open('a3_a6_inters.txt', 'w')
for i in range(len(Label)):
    n = int(index[i])
    output.write(Label[i] + ': ' + str(a3_a6[n]) + '+/-' + str(a3_a6_err[n]) + '\n')
output.close()

output = open('a6_a7_inters.txt', 'w')
for i in range(len(Label)):
    n = int(index[i])
    output.write(Label[i] + ': ' + str(a6_a7[n]) + '+/-' + str(a6_a7_err[n]) + '\n')
output.close()

#Plot comparison of all mutants
Label = ['Apo Open', 'Apo Closed', 'WT', 'F196A', 'L192A', 'L192F', 'L195A', 'L195F', 'L195N', 'S286A', 'F280Y', 'E276L', 'E276F', 'K279M', 'K279W', 'V287T']
index_AD = [30, 31, 0, 2, 4, 6, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28]
index_BBR = [30, 31, 1, 3, 5, 7, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29]

plot_mult_cmpr(Label, index_AD, 'AD', 'a3_a7_all', a3_a7, a3_a7_err, p_a3_a7_AD)
plot_mult_cmpr(Label, index_BBR, 'BBR', 'a3_a7_all', a3_a7, a3_a7_err, p_a3_a7_BBR)
plot_mult_cmpr(Label, index_AD, 'AD', 'a6_a7_all', a6_a7, a6_a7_err, p_a6_a7_AD)
plot_mult_cmpr(Label, index_BBR, 'BBR', 'a6_a7_all', a6_a7, a6_a7_err, p_a6_a7_BBR)
plot_mult_cmpr(Label, index_AD, 'AD', 'a3_a6_all', a3_a6, a3_a6_err, p_a3_a6_AD)
plot_mult_cmpr(Label, index_BBR, 'BBR', 'a3_a6_all', a3_a6, a3_a6_err, p_a3_a6_BBR)

#Plot comparison of mutants with similar BFE
Label_BFE, index_AD_BFE, index_BBR_BFE = [],[],[]
sel_index = [0, 1, 2, 4, 6]
for i in sel_index:
    Label_BFE.append(Label[i])
    index_AD_BFE.append(index_AD[i])
    index_BBR_BFE.append(index_BBR[i])


plot_mult_cmpr(Label_BFE, index_AD_BFE, 'AD', 'a3_a7_sim_BFE', a3_a7, a3_a7_err, p_a3_a7_AD)
plot_mult_cmpr(Label_BFE, index_BBR_BFE, 'BBR', 'a3_a7_sim_BFE', a3_a7, a3_a7_err, p_a3_a7_BBR)
plot_mult_cmpr(Label_BFE, index_AD_BFE, 'AD', 'a6_a7_sim_BFE', a6_a7, a6_a7_err, p_a6_a7_AD)
plot_mult_cmpr(Label_BFE, index_BBR_BFE, 'BBR', 'a6_a7_sim_BFE', a6_a7, a6_a7_err, p_a6_a7_BBR)
plot_mult_cmpr(Label_BFE, index_AD_BFE, 'AD', 'a3_a6_sim_BFE', a3_a6, a3_a6_err, p_a3_a6_AD)
plot_mult_cmpr(Label_BFE, index_BBR_BFE, 'BBR', 'a3_a6_sim_BFE', a3_a6, a3_a6_err, p_a3_a6_BBR)

#Set residue pairs
group_3 = [186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199, 200]
group_7 = [287, 288, 289, 290, 291, 292, 293, 294, 295, 296, 297, 298]
pair_a7_a3 = list(product(group_3, group_7))

a7_a3_inters_WT = np.zeros([len(pair_a7_a3), 2])

#Input Data for percent of presence of all a3 and a7 interactions
num = 0
for i in open("../../../WT/AD/analysis/a3_a7_inter_all.txt").readlines():
    a7_a3_inters_WT[num][0] = 100 - float(i)
    num += 1
num = 0
for i in open("../../../WT/BBR/analysis/a3_a7_inter_all.txt").readlines():
    a7_a3_inters_WT[num][1] = 100 - float(i)
    num += 1

corr_pairs = [(189, 291), (189, 292), (190, 291), (192, 291), (193, 291), (193, 292), (196, 291), (276, 291), (277, 291), (279, 291), (280, 291)]
mut_list = ['F196A', 'L192A', 'L192F', 'L195A', 'L195F', 'L195N', 'S286A', 'F280Y', 'E276L', 'E276F', 'K279M', 'K279W', 'V287T']
AD_corr = np.zeros([len(corr_pairs), len(mut_list)])
BBR_corr = np.zeros([len(corr_pairs), len(mut_list)])
for j in range(len(mut_list)):
    array = mut_inter('a3', 'a7', mut_list[j], pair_a7_a3)
    for i in range(len(pair_a7_a3)):
        if pair_a7_a3[i] in corr_pairs:
            for k in range(len(corr_pairs)):
                if pair_a7_a3[i] == corr_pairs[k]:
                    AD_corr[k][j], BBR_corr[k][j] = corr_inter(a7_a3_inters_WT, array, mut_list[j], i, pair_a7_a3)
        plot_inter(a7_a3_inters_WT, array, mut_list[j], 'a3', 'a7', pair_a7_a3)

#Set residue pairs
group_6 = [264, 265, 266, 267, 268, 269, 270, 271, 272, 273, 274, 275, 276, 277, 278, 278, 279, 280, 281]
pair_a7_a6 = list(product(group_6, group_7))

a6_a7_inters_WT = np.zeros([len(pair_a7_a6), 2])

#Input Data for percent of presence of all a3 and a7 interactions
num = 0
for i in open("../../../WT/AD/analysis/a6_a7_inter_all.txt").readlines():
    a6_a7_inters_WT[num][0] = float(i)
    num += 1
num = 0
for i in open("../../../WT/BBR/analysis/a6_a7_inter_all.txt").readlines():
    a6_a7_inters_WT[num][1] = float(i)
    num += 1

for j in range(len(mut_list)):
    array = mut_inter('a6', 'a7', mut_list[j], pair_a7_a6)
    for i in range(len(pair_a7_a6)):
        if pair_a7_a6[i] in corr_pairs:
            for k in range(len(corr_pairs)):
                if pair_a7_a6[i] == corr_pairs[k]:
                    AD_corr[k][j], BBR_corr[k][j] = corr_inter(a6_a7_inters_WT, array, mut_list[j], i, pair_a7_a6)
        plot_inter(a6_a7_inters_WT, array, mut_list[j], 'a6', 'a7', pair_a7_a6)

#Plot table comparing residue interactions to WT
ax = plt.figure(figsize=(10, 10), frameon=True) # no visible frame
ax = sns.heatmap(AD_corr, annot=False, cmap = 'bwr_r', cbar = True, cbar_kws={'label': 'Percentage Difference from WT'}, vmin = -200, vmax = 200, xticklabels = mut_list, yticklabels = corr_pairs)
ax.add_artist(lines.Line2D([0, 20], [7, 7], color = 'black', linestyle= '--', linewidth = 4))
plt.title('Bond Disruption Compared to WT')
plt.savefig('mutate_AD_corr_res.png')
plt.close()

ax = plt.figure(figsize=(10, 10), frameon=True) # no visible frame
ax = sns.heatmap(BBR_corr, annot=False, cmap = 'bwr_r', cbar = True, cbar_kws={'label': 'Percentage Difference from WT'}, vmin = -200, vmax = 200, xticklabels = mut_list, yticklabels = corr_pairs)
ax.add_artist(lines.Line2D([0, 20], [7, 7], color = 'black', linestyle= '--', linewidth = 4))
plt.title('Bond Disruption Compared to WT')
plt.savefig('mutate_BBR_corr_res.png')
plt.close()


