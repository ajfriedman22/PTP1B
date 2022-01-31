#Import Necessary Packages
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from itertools import product
def plot_func(dist, err, hel1, hel2, mut, p, p1, j, k):
    num = [5, 10, 15, 20, 25, 30]
    Method = ['Apo Open', 'Apo Closed', 'WT AD', 'WT BBR', mut + ' AD', mut + ' BBR']
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.set_title('Comparison of Distance b/w ' + hel1 + ' ' + hel2 + ' for ' + mut)    
    ax1.set_ylabel('Distance b/w Residues')
    ax1.bar(num, dist[[0, 1, 2, 3, j, k]], color = ['black', 'black', 'darkblue', 'darkred', 'blue', 'red'], width=4.5)
    plt.errorbar(num, dist[[0, 1, 2, 3, j, k]], yerr= err[[0, 1, 2, 3, j, k]], fmt='o', color='black')
    plt.xticks(num, Method, fontsize=8)
    if p < 0.05 and p > 0.01:
        x1, x2 = 15, 25 #Columns for WT and mutant AD
        y, h, col = (1.1*dist[[2, j]].max()), 1, 'b'
        plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
        plt.text((x1+x2)*0.5, y+h, "*" , ha='center', va='bottom', color=col)
    if p < 0.01 and p > 0.001:
        x1, x2 = 15, 25 #Columns for WT and mutant AD
        y, h, col = (1.1*dist[[2, j]].max()), 1, 'b'
        plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
        plt.text((x1+x2)*0.5, y+h, "**" , ha='center', va='bottom', color=col)
    if p < 0.001:
        x1, x2 = 15, 25 #Columns for WT and mutant AD
        y, h, col = (1.1*dist[[2, j]].max()), 1, 'b'
        plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
        plt.text((x1+x2)*0.5, y+h, "***" , ha='center', va='bottom', color=col)
    if p1 < 0.05 and p1 > 0.01:
        x1, x2 = 20, 30 #Columns for WT and mutant BBR
        y, h, col = (1.1*dist[[3, k]].max()), 1, 'r'
        plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
        plt.text((x1+x2)*0.5, y+h, "*" , ha='center', va='bottom', color=col)
    if p1 < 0.01 and p1 > 0.001:
        x1, x2 = 20, 30 #Columns for WT and mutant BBR
        y, h, col = (1.1*dist[[3, k]].max()), 1, 'r'
        plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
        plt.text((x1+x2)*0.5, y+h, "**" , ha='center', va='bottom', color=col)
    if p1 < 0.001:
        x1, x2 = 20, 30 #Columns for WT and mutant BBR
        y, h, col = (1.1*dist[[3, k]].max()), 1, 'r'
        plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
        plt.text((x1+x2)*0.5, y+h, "***" , ha='center', va='bottom', color=col)
    fig.savefig(hel1 + '_' + hel2 + '_' + mut + '_inter.png')
    plt.close(fig)


#Make open arrays for time and atomic distances
d200_282_WT_AD, d200_282_WT_BBR, d200_282_F196A_AD, d200_282_F196A_BBR = [],[],[],[]
d200_282_L192A_AD, d200_282_L192A_BBR, d200_282_L192F_AD, d200_282_L192F_BBR = [],[],[],[]
d200_282_L195A_AD, d200_282_L195A_BBR, d200_282_L195F_AD, d200_282_L195F_BBR, d200_282_L195N_AD, d200_282_L195N_BBR = [],[],[],[],[],[]
d200_282_Apo_open, d200_282_Apo_close = [],[]
d200_282_S286A_AD, d200_282_S286A_BBR, d200_282_F280Y_AD, d200_282_F280Y_BBR = [],[],[],[]
d200_282_E276L_AD, d200_282_E276L_BBR, d200_282_E276F_AD, d200_282_E276F_BBR = [],[],[],[]
d200_282_K279M_AD, d200_282_K279M_BBR, d200_282_K279W_AD, d200_282_K279W_BBR = [],[],[],[]
d200_282_V287T_AD, d200_282_V287T_BBR = [],[]

d200_287_WT_AD, d200_287_WT_BBR, d200_287_F196A_AD, d200_287_F196A_BBR = [],[],[],[]
d200_287_L192A_AD, d200_287_L192A_BBR, d200_287_L192F_AD, d200_287_L192F_BBR = [],[],[],[]
d200_287_L195A_AD, d200_287_L195A_BBR, d200_287_L195F_AD, d200_287_L195F_BBR, d200_287_L195N_AD, d200_287_L195N_BBR = [],[],[],[],[],[]
d200_287_Apo_open, d200_287_Apo_close = [],[]
d200_287_S286A_AD, d200_287_S286A_BBR, d200_287_F280Y_AD, d200_287_F280Y_BBR = [],[],[],[]
d200_287_E276L_AD, d200_287_E276L_BBR, d200_287_E276F_AD, d200_287_E276F_BBR = [],[],[],[]
d200_287_K279M_AD, d200_287_K279M_BBR, d200_287_K279W_AD, d200_287_K279W_BBR = [],[],[],[]
d200_287_V287T_AD, d200_287_V287T_BBR = [],[]

d276_292_WT_AD, d276_292_WT_BBR, d276_292_F196A_AD, d276_292_F196A_BBR = [],[],[],[]
d276_292_L192A_AD, d276_292_L192A_BBR, d276_292_L192F_AD, d276_292_L192F_BBR = [],[],[],[]
d276_292_L195A_AD, d276_292_L195A_BBR, d276_292_L195F_AD, d276_292_L195F_BBR, d276_292_L195N_AD, d276_292_L195N_BBR = [],[],[],[],[],[]
d276_292_Apo_open, d276_292_Apo_close = [],[]
d276_292_S286A_AD, d276_292_S286A_BBR, d276_292_F280Y_AD, d276_292_F280Y_BBR = [],[],[],[]
d276_292_E276L_AD, d276_292_E276L_BBR, d276_292_E276F_AD, d276_292_E276F_BBR = [],[],[],[]
d276_292_K279M_AD, d276_292_K279M_BBR, d276_292_K279W_AD, d276_292_K279W_BBR = [],[],[],[]
d276_292_V287T_AD, d276_292_V287T_BBR = [],[]

#Input Data for a3 and a6 residue dist
for i in open("../../../WT/AD/200_282_dist.txt").readlines():
    d200_282_WT_AD.append(float(i))
for i in open("../../../WT/BBR/200_282_dist.txt").readlines():
    d200_282_WT_BBR.append(float(i))
for i in open("../../../F196A/AD/200_282_dist.txt").readlines():
    d200_282_F196A_AD.append(float(i))
for i in open("../../../F196A/BBR/200_282_dist.txt").readlines():
    d200_282_F196A_BBR.append(float(i))
for i in open("../../../L192A/AD/200_282_dist.txt").readlines():
    d200_282_L192A_AD.append(float(i))
for i in open("../../../L192A/BBR/200_282_dist.txt").readlines():
    d200_282_L192A_BBR.append(float(i))
for i in open("../../../L192F/AD/200_282_dist.txt").readlines():
    d200_282_L192F_AD.append(float(i))
for i in open("../../../L192F/BBR/200_282_dist.txt").readlines():
    d200_282_L192F_BBR.append(float(i))
for i in open("../../../L195A/AD/200_282_dist.txt").readlines():
    d200_282_L195A_AD.append(float(i))
for i in open("../../../L195A/BBR/200_282_dist.txt").readlines():
    d200_282_L195A_BBR.append(float(i))
for i in open("../../../L195F/AD/200_282_dist.txt").readlines():
    d200_282_L195F_AD.append(float(i))
for i in open("../../../L195F/BBR/200_282_dist.txt").readlines():
    d200_282_L195F_BBR.append(float(i))
for i in open("../../../L195N/AD/200_282_dist.txt").readlines():
    d200_282_L195N_AD.append(float(i))
for i in open("../../../L195N/BBR/200_282_dist.txt").readlines():
    d200_282_L195N_BBR.append(float(i))
for i in open("../../../S286A/AD/200_282_dist.txt").readlines():
    d200_282_S286A_AD.append(float(i))
for i in open("../../../S286A/BBR/200_282_dist.txt").readlines():
    d200_282_S286A_BBR.append(float(i))
for i in open("../../../F280Y/AD/200_282_dist.txt").readlines():
    d200_282_F280Y_AD.append(float(i))
for i in open("../../../F280Y/BBR/200_282_dist.txt").readlines():
    d200_282_F280Y_BBR.append(float(i))
for i in open("../../../E276L/AD/200_282_dist.txt").readlines():
    d200_282_E276L_AD.append(float(i))
for i in open("../../../E276L/BBR/200_282_dist.txt").readlines():
    d200_282_E276L_BBR.append(float(i))
for i in open("../../../E276F/AD/200_282_dist.txt").readlines():
    d200_282_E276F_AD.append(float(i))
for i in open("../../../E276F/BBR/200_282_dist.txt").readlines():
    d200_282_E276F_BBR.append(float(i))
for i in open("../../../K279M/AD/200_282_dist.txt").readlines():
    d200_282_K279M_AD.append(float(i))
for i in open("../../../K279M/BBR/200_282_dist.txt").readlines():
    d200_282_K279M_BBR.append(float(i))
for i in open("../../../K279W/AD/200_282_dist.txt").readlines():
    d200_282_K279W_AD.append(float(i))
for i in open("../../../K279W/BBR/200_282_dist.txt").readlines():
    d200_282_K279W_BBR.append(float(i))
for i in open("../../../V287T/AD/200_282_dist.txt").readlines():
    d200_282_V287T_AD.append(float(i))
for i in open("../../../V287T/BBR/200_282_dist.txt").readlines():
    d200_282_V287T_BBR.append(float(i))
for i in open("../../../../rebuild_a7/200_282_dist.txt").readlines():
    d200_282_Apo_open.append(float(i))
for i in open("../../../../Apo_dis/config9/200_282_dist.txt").readlines():
    d200_282_Apo_open.append(float(i))
for i in open("../../../../Apo_dis/config11/200_282_dist.txt").readlines():
    d200_282_Apo_open.append(float(i))
for i in open("../../../../1sug/200_282_dist.txt").readlines():
    d200_282_Apo_close.append(float(i))
for i in open("../../../../1sug2/200_282_dist.txt").readlines():
    d200_282_Apo_close.append(float(i))
for i in open("../../../../1sug3/200_282_dist.txt").readlines():
    d200_282_Apo_open.append(float(i))

#Input Data for a3 and a7 residue dist
for i in open("../../../WT/AD/200_287_dist.txt").readlines():
    d200_287_WT_AD.append(float(i))
for i in open("../../../WT/BBR/200_287_dist.txt").readlines():
    d200_287_WT_BBR.append(float(i))
for i in open("../../../F196A/AD/200_287_dist.txt").readlines():
    d200_287_F196A_AD.append(float(i))
for i in open("../../../F196A/BBR/200_287_dist.txt").readlines():
    d200_287_F196A_BBR.append(float(i))
for i in open("../../../L192A/AD/200_287_dist.txt").readlines():
    d200_287_L192A_AD.append(float(i))
for i in open("../../../L192A/BBR/200_287_dist.txt").readlines():
    d200_287_L192A_BBR.append(float(i))
for i in open("../../../L192F/AD/200_287_dist.txt").readlines():
    d200_287_L192F_AD.append(float(i))
for i in open("../../../L192F/BBR/200_287_dist.txt").readlines():
    d200_287_L192F_BBR.append(float(i))
for i in open("../../../L195A/AD/200_287_dist.txt").readlines():
    d200_287_L195A_AD.append(float(i))
for i in open("../../../L195A/BBR/200_287_dist.txt").readlines():
    d200_287_L195A_BBR.append(float(i))
for i in open("../../../L195F/AD/200_287_dist.txt").readlines():
    d200_287_L195F_AD.append(float(i))
for i in open("../../../L195F/BBR/200_287_dist.txt").readlines():
    d200_287_L195F_BBR.append(float(i))
for i in open("../../../L195N/AD/200_287_dist.txt").readlines():
    d200_287_L195N_AD.append(float(i))
for i in open("../../../L195N/BBR/200_287_dist.txt").readlines():
    d200_287_L195N_BBR.append(float(i))
for i in open("../../../S286A/AD/200_287_dist.txt").readlines():
    d200_287_S286A_AD.append(float(i))
for i in open("../../../S286A/BBR/200_287_dist.txt").readlines():
    d200_287_S286A_BBR.append(float(i))
for i in open("../../../F280Y/AD/200_287_dist.txt").readlines():
    d200_287_F280Y_AD.append(float(i))
for i in open("../../../F280Y/BBR/200_287_dist.txt").readlines():
    d200_287_F280Y_BBR.append(float(i))
for i in open("../../../E276L/AD/200_287_dist.txt").readlines():
    d200_287_E276L_AD.append(float(i))
for i in open("../../../E276L/BBR/200_287_dist.txt").readlines():
    d200_287_E276L_BBR.append(float(i))
for i in open("../../../E276F/AD/200_287_dist.txt").readlines():
    d200_287_E276F_AD.append(float(i))
for i in open("../../../E276F/BBR/200_287_dist.txt").readlines():
    d200_287_E276F_BBR.append(float(i))
for i in open("../../../K279M/AD/200_287_dist.txt").readlines():
    d200_287_K279M_AD.append(float(i))
for i in open("../../../K279M/BBR/200_287_dist.txt").readlines():
    d200_287_K279M_BBR.append(float(i))
for i in open("../../../K279W/AD/200_287_dist.txt").readlines():
    d200_287_K279W_AD.append(float(i))
for i in open("../../../K279W/BBR/200_287_dist.txt").readlines():
    d200_287_K279W_BBR.append(float(i))
for i in open("../../../V287T/AD/200_287_dist.txt").readlines():
    d200_287_V287T_AD.append(float(i))
for i in open("../../../V287T/BBR/200_287_dist.txt").readlines():
    d200_287_V287T_BBR.append(float(i))
for i in open("../../../../rebuild_a7/200_287_dist.txt").readlines():
    d200_287_Apo_open.append(float(i))
for i in open("../../../../Apo_dis/config9/200_287_dist.txt").readlines():
    d200_287_Apo_open.append(float(i))
for i in open("../../../../Apo_dis/config11/200_287_dist.txt").readlines():
    d200_287_Apo_open.append(float(i))
for i in open("../../../../1sug/200_287_dist.txt").readlines():
    d200_287_Apo_close.append(float(i))
for i in open("../../../../1sug2/200_287_dist.txt").readlines():
    d200_287_Apo_close.append(float(i))
for i in open("../../../../1sug3/200_287_dist.txt").readlines():
    d200_287_Apo_open.append(float(i))

#Input Data for a6 and a7 residue dist
for i in open("../../../WT/AD/276_292_dist.txt").readlines():
    d276_292_WT_AD.append(float(i))
for i in open("../../../WT/BBR/276_292_dist.txt").readlines():
    d276_292_WT_BBR.append(float(i))
for i in open("../../../F196A/AD/276_292_dist.txt").readlines():
    d276_292_F196A_AD.append(float(i))
for i in open("../../../F196A/BBR/276_292_dist.txt").readlines():
    d276_292_F196A_BBR.append(float(i))
for i in open("../../../L192A/AD/276_292_dist.txt").readlines():
    d276_292_L192A_AD.append(float(i))
for i in open("../../../L192A/BBR/276_292_dist.txt").readlines():
    d276_292_L192A_BBR.append(float(i))
for i in open("../../../L192F/AD/276_292_dist.txt").readlines():
    d276_292_L192F_AD.append(float(i))
for i in open("../../../L192F/BBR/276_292_dist.txt").readlines():
    d276_292_L192F_BBR.append(float(i))
for i in open("../../../L195A/AD/276_292_dist.txt").readlines():
    d276_292_L195A_AD.append(float(i))
for i in open("../../../L195A/BBR/276_292_dist.txt").readlines():
    d276_292_L195A_BBR.append(float(i))
for i in open("../../../L195F/AD/276_292_dist.txt").readlines():
    d276_292_L195F_AD.append(float(i))
for i in open("../../../L195F/BBR/276_292_dist.txt").readlines():
    d276_292_L195F_BBR.append(float(i))
for i in open("../../../L195N/AD/276_292_dist.txt").readlines():
    d276_292_L195N_AD.append(float(i))
for i in open("../../../L195N/BBR/276_292_dist.txt").readlines():
    d276_292_L195N_BBR.append(float(i))
for i in open("../../../S286A/AD/276_292_dist.txt").readlines():
    d276_292_S286A_AD.append(float(i))
for i in open("../../../S286A/BBR/276_292_dist.txt").readlines():
    d276_292_S286A_BBR.append(float(i))
for i in open("../../../F280Y/AD/276_292_dist.txt").readlines():
    d276_292_F280Y_AD.append(float(i))
for i in open("../../../F280Y/BBR/276_292_dist.txt").readlines():
    d276_292_F280Y_BBR.append(float(i))
for i in open("../../../E276L/AD/276_292_dist.txt").readlines():
    d276_292_E276L_AD.append(float(i))
for i in open("../../../E276L/BBR/276_292_dist.txt").readlines():
    d276_292_E276L_BBR.append(float(i))
for i in open("../../../E276F/AD/276_292_dist.txt").readlines():
    d276_292_E276F_AD.append(float(i))
for i in open("../../../E276F/BBR/276_292_dist.txt").readlines():
    d276_292_E276F_BBR.append(float(i))
for i in open("../../../K279M/AD/276_292_dist.txt").readlines():
    d276_292_K279M_AD.append(float(i))
for i in open("../../../K279M/BBR/276_292_dist.txt").readlines():
    d276_292_K279M_BBR.append(float(i))
for i in open("../../../K279W/AD/276_292_dist.txt").readlines():
    d276_292_K279W_AD.append(float(i))
for i in open("../../../K279W/BBR/276_292_dist.txt").readlines():
    d276_292_K279W_BBR.append(float(i))
for i in open("../../../V287T/AD/276_292_dist.txt").readlines():
    d276_292_V287T_AD.append(float(i))
for i in open("../../../V287T/BBR/276_292_dist.txt").readlines():
    d276_292_V287T_BBR.append(float(i))
for i in open("../../../../rebuild_a7/276_292_dist.txt").readlines():
    d276_292_Apo_open.append(float(i))
for i in open("../../../../Apo_dis/config9/276_292_dist.txt").readlines():
    d276_292_Apo_open.append(float(i))
for i in open("../../../../Apo_dis/config11/276_292_dist.txt").readlines():
    d276_292_Apo_open.append(float(i))
for i in open("../../../../1sug/276_292_dist.txt").readlines():
    d276_292_Apo_close.append(float(i))
for i in open("../../../../1sug2/276_292_dist.txt").readlines():
    d276_292_Apo_close.append(float(i))
for i in open("../../../../1sug3/276_292_dist.txt").readlines():
    d276_292_Apo_open.append(float(i))

#Calculate mean and sem for interactions
d200_282 = np.array([np.mean(d200_282_Apo_open), np.mean(d200_282_Apo_close), np.mean(d200_282_WT_AD), 
    np.mean(d200_282_WT_BBR), np.mean(d200_282_F196A_AD), np.mean(d200_282_F196A_BBR), np.mean(d200_282_L192A_AD), 
    np.mean(d200_282_L192A_BBR), np.mean(d200_282_L192F_AD), np.mean(d200_282_L192F_BBR), 
    np.mean(d200_282_L195A_AD), np.mean(d200_282_L195A_BBR), np.mean(d200_282_L195F_AD), 
    np.mean(d200_282_L195F_BBR), np.mean(d200_282_L195N_AD), np.mean(d200_282_L195N_BBR), 
    np.mean(d200_282_S286A_AD),  np.mean(d200_282_S286A_BBR), np.mean(d200_282_F280Y_AD), 
    np.mean(d200_282_F280Y_BBR), np.mean(d200_282_E276L_AD), np.mean(d200_282_E276L_BBR),
    np.mean(d200_282_E276F_AD), np.mean(d200_282_E276F_BBR), np.mean(d200_282_K279M_AD),
    np.mean(d200_282_K279M_BBR), np.mean(d200_282_K279W_AD), np.mean(d200_282_K279W_BBR),
    np.mean(d200_282_V287T_AD), np.mean(d200_282_V287T_BBR)])

d200_287 = np.array([np.mean(d200_287_Apo_open), np.mean(d200_287_Apo_close), np.mean(d200_287_WT_AD), 
    np.mean(d200_287_WT_BBR), np.mean(d200_287_F196A_AD), np.mean(d200_287_F196A_BBR), 
    np.mean(d200_287_L192A_AD), np.mean(d200_287_L192A_BBR), np.mean(d200_287_L192F_AD), 
    np.mean(d200_287_L192F_BBR), np.mean(d200_287_L195A_AD), np.mean(d200_287_L195A_BBR), 
    np.mean(d200_287_L195F_AD), np.mean(d200_287_L195F_BBR), 
    np.mean(d200_287_L195N_AD), np.mean(d200_287_L195N_BBR),
    np.mean(d200_287_S286A_AD),  np.mean(d200_287_S286A_BBR), np.mean(d200_287_F280Y_AD), 
    np.mean(d200_287_F280Y_BBR), np.mean(d200_287_E276L_AD), np.mean(d200_287_E276L_BBR),
    np.mean(d200_287_E276F_AD), np.mean(d200_287_E276F_BBR), np.mean(d200_287_K279M_AD),
    np.mean(d200_287_K279M_BBR), np.mean(d200_287_K279W_AD), np.mean(d200_287_K279W_BBR),
    np.mean(d200_287_V287T_AD), np.mean(d200_287_V287T_BBR)])

d276_292 = np.array([np.mean(d276_292_Apo_open), np.mean(d276_292_Apo_close), np.mean(d276_292_WT_AD), 
    np.mean(d276_292_WT_BBR), np.mean(d276_292_F196A_AD), np.mean(d276_292_F196A_BBR), 
    np.mean(d276_292_L192A_AD), np.mean(d276_292_L192A_BBR), np.mean(d276_292_L192F_AD), 
    np.mean(d276_292_L192F_BBR), np.mean(d276_292_L195A_AD), np.mean(d276_292_L195A_BBR), 
    np.mean(d276_292_L195F_AD), np.mean(d276_292_L195F_BBR), 
    np.mean(d276_292_L195N_AD), np.mean(d276_292_L195N_BBR),
    np.mean(d276_292_S286A_AD),  np.mean(d276_292_S286A_BBR), np.mean(d276_292_F280Y_AD), 
    np.mean(d276_292_F280Y_BBR), np.mean(d276_292_E276L_AD), np.mean(d276_292_E276L_BBR),
    np.mean(d276_292_E276F_AD), np.mean(d276_292_E276F_BBR), np.mean(d276_292_K279M_AD),
    np.mean(d276_292_K279M_BBR), np.mean(d276_292_K279W_AD), np.mean(d276_292_K279W_BBR),
    np.mean(d276_292_V287T_AD), np.mean(d276_292_V287T_BBR)])

d200_282_err = np.array([stats.sem(d200_282_Apo_open), stats.sem(d200_282_Apo_close), stats.sem(d200_282_WT_AD),
    stats.sem(d200_282_WT_BBR), stats.sem(d200_282_F196A_AD), stats.sem(d200_282_F196A_BBR), 
    stats.sem(d200_282_L192A_AD), stats.sem(d200_282_L192A_BBR), stats.sem(d200_282_L192F_AD), 
    stats.sem(d200_282_L192F_BBR), stats.sem(d200_282_L195A_AD), stats.sem(d200_282_L195A_BBR), 
    stats.sem(d200_282_L195F_AD), stats.sem(d200_282_L195F_BBR), stats.sem(d200_282_L195N_AD), 
    stats.sem(d200_282_L195N_BBR),
    stats.sem(d200_282_S286A_AD),  stats.sem(d200_282_S286A_BBR), stats.sem(d200_282_F280Y_AD), 
    stats.sem(d200_282_F280Y_BBR), stats.sem(d200_282_E276L_AD), stats.sem(d200_282_E276L_BBR),
    stats.sem(d200_282_E276F_AD), stats.sem(d200_282_E276F_BBR), stats.sem(d200_282_K279M_AD),
    stats.sem(d200_282_K279M_BBR), stats.sem(d200_282_K279W_AD), stats.sem(d200_282_K279W_BBR),
    stats.sem(d200_282_V287T_AD), stats.sem(d200_282_V287T_BBR)])


d200_287_err = np.array([stats.sem(d200_287_Apo_open), stats.sem(d200_287_Apo_close), stats.sem(d200_287_WT_AD),
    stats.sem(d200_287_WT_BBR), stats.sem(d200_287_F196A_AD), stats.sem(d200_287_F196A_BBR), 
    stats.sem(d200_287_L192A_AD), stats.sem(d200_287_L192A_BBR), stats.sem(d200_287_L192F_AD), 
    stats.sem(d200_287_L192F_BBR), stats.sem(d200_287_L195A_AD), stats.sem(d200_287_L195A_BBR), 
    stats.sem(d200_287_L195F_AD), stats.sem(d200_287_L195F_BBR), stats.sem(d200_287_L195N_AD), 
    stats.sem(d200_287_L195N_BBR),
    stats.sem(d200_287_S286A_AD),  stats.sem(d200_287_S286A_BBR), stats.sem(d200_287_F280Y_AD), 
    stats.sem(d200_287_F280Y_BBR), stats.sem(d200_287_E276L_AD), stats.sem(d200_287_E276L_BBR),
    stats.sem(d200_287_E276F_AD), stats.sem(d200_287_E276F_BBR), stats.sem(d200_287_K279M_AD),
    stats.sem(d200_287_K279M_BBR), stats.sem(d200_287_K279W_AD), stats.sem(d200_287_K279W_BBR),
    stats.sem(d200_287_V287T_AD), stats.sem(d200_287_V287T_BBR)])


d276_292_err = np.array([stats.sem(d276_292_Apo_open), stats.sem(d276_292_Apo_close), stats.sem(d276_292_WT_AD),
    stats.sem(d276_292_WT_BBR), stats.sem(d276_292_F196A_AD), stats.sem(d276_292_F196A_BBR), 
    stats.sem(d276_292_L192A_AD), stats.sem(d276_292_L192A_BBR), stats.sem(d276_292_L192F_AD), 
    stats.sem(d276_292_L192F_BBR), stats.sem(d276_292_L195A_AD), stats.sem(d276_292_L195A_BBR), 
    stats.sem(d276_292_L195F_AD), stats.sem(d276_292_L195F_BBR), 
    stats.sem(d276_292_L195N_AD), stats.sem(d276_292_L195N_BBR),
    stats.sem(d276_292_S286A_AD),  stats.sem(d276_292_S286A_BBR), stats.sem(d276_292_F280Y_AD), 
    stats.sem(d276_292_F280Y_BBR), stats.sem(d276_292_E276L_AD), stats.sem(d276_292_E276L_BBR),
    stats.sem(d276_292_E276F_AD), stats.sem(d276_292_E276F_BBR), stats.sem(d276_292_K279M_AD),
    stats.sem(d276_292_K279M_BBR), stats.sem(d276_292_K279W_AD), stats.sem(d276_292_K279W_BBR),
    stats.sem(d276_292_V287T_AD), stats.sem(d276_292_V287T_BBR)])


#Run t-test between groups
st, p1 = stats.ttest_ind(d200_282_WT_AD, d200_282_F196A_AD, equal_var = False) #Welch's t-test b/w WT and F196A AD
st, p2 = stats.ttest_ind(d200_282_WT_AD, d200_282_L192A_AD, equal_var = False) #Welch's t-test b/w WT and L192A AD
st, p3 = stats.ttest_ind(d200_282_WT_AD, d200_282_L192F_AD, equal_var = False) #Welch's t-test b/w WT and L192F AD
st, p4 = stats.ttest_ind(d200_282_WT_AD, d200_282_L195A_AD, equal_var = False) #Welch's t-test b/w WT and L195A AD
st, p5 = stats.ttest_ind(d200_282_WT_AD, d200_282_L195F_AD, equal_var = False) #Welch's t-test b/w WT and L195F AD
st, p6 = stats.ttest_ind(d200_282_WT_AD, d200_282_L195N_AD, equal_var = False) #Welch's t-test b/w WT and L195N AD

st, p7 = stats.ttest_ind(d200_287_WT_AD, d200_287_F196A_AD, equal_var = False) #Welch's t-test b/w WT and F196A AD
st, p8 = stats.ttest_ind(d200_287_WT_AD, d200_287_L192A_AD, equal_var = False) #Welch's t-test b/w WT and L192A AD
st, p9 = stats.ttest_ind(d200_287_WT_AD, d200_287_L192F_AD, equal_var = False) #Welch's t-test b/w WT and L192F AD
st, p10 = stats.ttest_ind(d200_287_WT_AD, d200_287_L195A_AD, equal_var = False) #Welch's t-test b/w WT and L195A AD
st, p11 = stats.ttest_ind(d200_287_WT_AD, d200_287_L195F_AD, equal_var = False) #Welch's t-test b/w WT and L195F AD
st, p12 = stats.ttest_ind(d200_287_WT_AD, d200_287_L195N_AD, equal_var = False) #Welch's t-test b/w WT and L195N AD

st, p13 = stats.ttest_ind(d276_292_WT_AD, d276_292_F196A_AD, equal_var = False) #Welch's t-test b/w WT and F196A AD
st, p14 = stats.ttest_ind(d276_292_WT_AD, d276_292_L192A_AD, equal_var = False) #Welch's t-test b/w WT and L192A AD
st, p15 = stats.ttest_ind(d276_292_WT_AD, d276_292_L192F_AD, equal_var = False) #Welch's t-test b/w WT and L192F AD
st, p16 = stats.ttest_ind(d276_292_WT_AD, d276_292_L195A_AD, equal_var = False) #Welch's t-test b/w WT and L195A AD
st, p17 = stats.ttest_ind(d276_292_WT_AD, d276_292_L195F_AD, equal_var = False) #Welch's t-test b/w WT and L195F AD
st, p18 = stats.ttest_ind(d276_292_WT_AD, d276_292_L195N_AD, equal_var = False) #Welch's t-test b/w WT and L195N AD

st, p1B = stats.ttest_ind(d200_282_WT_BBR, d200_282_F196A_AD, equal_var = False) #Welch's t-test b/w WT and F196A AD
st, p2B = stats.ttest_ind(d200_282_WT_BBR, d200_282_L192A_AD, equal_var = False) #Welch's t-test b/w WT and L192A AD
st, p3B = stats.ttest_ind(d200_282_WT_BBR, d200_282_L192F_AD, equal_var = False) #Welch's t-test b/w WT and L192F AD
st, p4B = stats.ttest_ind(d200_282_WT_BBR, d200_282_L195A_AD, equal_var = False) #Welch's t-test b/w WT and L195A AD
st, p5B = stats.ttest_ind(d200_282_WT_BBR, d200_282_L195F_AD, equal_var = False) #Welch's t-test b/w WT and L195F AD
st, p6B = stats.ttest_ind(d200_282_WT_BBR, d200_282_L195N_AD, equal_var = False) #Welch's t-test b/w WT and L195N AD

st, p7B = stats.ttest_ind(d200_287_WT_BBR, d200_287_F196A_BBR, equal_var = False) #Welch's t-test b/w WT and F196A AD
st, p8B = stats.ttest_ind(d200_287_WT_BBR, d200_287_L192A_BBR, equal_var = False) #Welch's t-test b/w WT and L192A AD
st, p9B = stats.ttest_ind(d200_287_WT_BBR, d200_287_L192F_BBR, equal_var = False) #Welch's t-test b/w WT and L192F AD
st, p10B = stats.ttest_ind(d200_287_WT_BBR, d200_287_L195A_BBR, equal_var = False) #Welch's t-test b/w WT and L195A AD
st, p11B = stats.ttest_ind(d200_287_WT_BBR, d200_287_L195F_BBR, equal_var = False) #Welch's t-test b/w WT and L195F AD
st, p12B = stats.ttest_ind(d200_287_WT_BBR, d200_287_L195N_BBR, equal_var = False) #Welch's t-test b/w WT and L195N AD

st, p13B = stats.ttest_ind(d276_292_WT_BBR, d276_292_F196A_BBR, equal_var = False) #Welch's t-test b/w WT and F196A AD
st, p14B = stats.ttest_ind(d276_292_WT_BBR, d276_292_L192A_BBR, equal_var = False) #Welch's t-test b/w WT and L192A AD
st, p15B = stats.ttest_ind(d276_292_WT_BBR, d276_292_L192F_BBR, equal_var = False) #Welch's t-test b/w WT and L192F AD
st, p16B = stats.ttest_ind(d276_292_WT_BBR, d276_292_L195A_BBR, equal_var = False) #Welch's t-test b/w WT and L195A AD
st, p17B = stats.ttest_ind(d276_292_WT_BBR, d276_292_L195F_BBR, equal_var = False) #Welch's t-test b/w WT and L195F AD
st, p18B = stats.ttest_ind(d276_292_WT_BBR, d276_292_L195N_BBR, equal_var = False) #Welch's t-test b/w WT and L195N AD

st, p19 = stats.ttest_ind(d200_282_WT_AD, d200_282_S286A_AD, equal_var = False) #Welch's t-test b/w WT and F196A AD
st, p20 = stats.ttest_ind(d200_282_WT_AD, d200_282_F280Y_AD, equal_var = False) #Welch's t-test b/w WT and L192A AD
st, p21 = stats.ttest_ind(d200_282_WT_AD, d200_282_E276L_AD, equal_var = False) #Welch's t-test b/w WT and L192F AD
st, p22 = stats.ttest_ind(d200_282_WT_AD, d200_282_E276F_AD, equal_var = False) #Welch's t-test b/w WT and L195A AD
st, p23 = stats.ttest_ind(d200_282_WT_AD, d200_282_K279M_AD, equal_var = False) #Welch's t-test b/w WT and L195F AD
st, p24 = stats.ttest_ind(d200_282_WT_AD, d200_282_K279W_AD, equal_var = False) #Welch's t-test b/w WT and L195N AD
st, p25 = stats.ttest_ind(d200_282_WT_AD, d200_282_V287T_AD, equal_var = False) #Welch's t-test b/w WT and L195N AD

st, p26 = stats.ttest_ind(d200_287_WT_AD, d200_287_S286A_AD, equal_var = False) #Welch's t-test b/w WT and F196A AD
st, p27 = stats.ttest_ind(d200_287_WT_AD, d200_287_F280Y_AD, equal_var = False) #Welch's t-test b/w WT and L192A AD
st, p28 = stats.ttest_ind(d200_287_WT_AD, d200_287_E276L_AD, equal_var = False) #Welch's t-test b/w WT and L192F AD
st, p29 = stats.ttest_ind(d200_287_WT_AD, d200_287_E276F_AD, equal_var = False) #Welch's t-test b/w WT and L195A AD
st, p30 = stats.ttest_ind(d200_287_WT_AD, d200_287_K279M_AD, equal_var = False) #Welch's t-test b/w WT and L195F AD
st, p31 = stats.ttest_ind(d200_287_WT_AD, d200_287_K279W_AD, equal_var = False) #Welch's t-test b/w WT and L195N AD
st, p32 = stats.ttest_ind(d200_287_WT_AD, d200_287_V287T_AD, equal_var = False) #Welch's t-test b/w WT and L195N AD

st, p33 = stats.ttest_ind(d276_292_WT_AD, d276_292_S286A_AD, equal_var = False) #Welch's t-test b/w WT and F196A AD
st, p34 = stats.ttest_ind(d276_292_WT_AD, d276_292_F280Y_AD, equal_var = False) #Welch's t-test b/w WT and L192A AD
st, p35 = stats.ttest_ind(d276_292_WT_AD, d276_292_E276L_AD, equal_var = False) #Welch's t-test b/w WT and L192F AD
st, p36 = stats.ttest_ind(d276_292_WT_AD, d276_292_E276F_AD, equal_var = False) #Welch's t-test b/w WT and L195A AD
st, p37 = stats.ttest_ind(d276_292_WT_AD, d276_292_K279M_AD, equal_var = False) #Welch's t-test b/w WT and L195F AD
st, p38 = stats.ttest_ind(d276_292_WT_AD, d276_292_K279W_AD, equal_var = False) #Welch's t-test b/w WT and L195N AD
st, p39 = stats.ttest_ind(d276_292_WT_AD, d276_292_V287T_AD, equal_var = False) #Welch's t-test b/w WT and L195N AD

st, p19B = stats.ttest_ind(d200_282_WT_BBR, d200_282_S286A_BBR, equal_var = False) #Welch's t-test b/w WT and F196A AD
st, p20B = stats.ttest_ind(d200_282_WT_BBR, d200_282_F280Y_BBR, equal_var = False) #Welch's t-test b/w WT and L192A AD
st, p21B = stats.ttest_ind(d200_282_WT_BBR, d200_282_E276L_BBR, equal_var = False) #Welch's t-test b/w WT and L192F AD
st, p22B = stats.ttest_ind(d200_282_WT_BBR, d200_282_E276F_BBR, equal_var = False) #Welch's t-test b/w WT and L195A AD
st, p23B = stats.ttest_ind(d200_282_WT_BBR, d200_282_K279M_BBR, equal_var = False) #Welch's t-test b/w WT and L195F AD
st, p24B = stats.ttest_ind(d200_282_WT_BBR, d200_282_K279W_BBR, equal_var = False) #Welch's t-test b/w WT and L195N AD
st, p25B = stats.ttest_ind(d200_282_WT_BBR, d200_282_V287T_BBR, equal_var = False) #Welch's t-test b/w WT and L195N AD

st, p26B = stats.ttest_ind(d200_287_WT_BBR, d200_287_S286A_BBR, equal_var = False) #Welch's t-test b/w WT and F196A AD
st, p27B = stats.ttest_ind(d200_287_WT_BBR, d200_287_F280Y_BBR, equal_var = False) #Welch's t-test b/w WT and L192A AD
st, p28B = stats.ttest_ind(d200_287_WT_BBR, d200_287_E276L_BBR, equal_var = False) #Welch's t-test b/w WT and L192F AD
st, p29B = stats.ttest_ind(d200_287_WT_BBR, d200_287_E276F_BBR, equal_var = False) #Welch's t-test b/w WT and L195A AD
st, p30B = stats.ttest_ind(d200_287_WT_BBR, d200_287_K279M_BBR, equal_var = False) #Welch's t-test b/w WT and L195F AD
st, p31B = stats.ttest_ind(d200_287_WT_BBR, d200_287_K279W_BBR, equal_var = False) #Welch's t-test b/w WT and L195N AD
st, p32B = stats.ttest_ind(d200_287_WT_BBR, d200_287_V287T_BBR, equal_var = False) #Welch's t-test b/w WT and L195N AD

st, p33B = stats.ttest_ind(d276_292_WT_BBR, d276_292_S286A_BBR, equal_var = False) #Welch's t-test b/w WT and F196A AD
st, p34B = stats.ttest_ind(d276_292_WT_BBR, d276_292_F280Y_BBR, equal_var = False) #Welch's t-test b/w WT and L192A AD
st, p35B = stats.ttest_ind(d276_292_WT_BBR, d276_292_E276L_BBR, equal_var = False) #Welch's t-test b/w WT and L192F AD
st, p36B = stats.ttest_ind(d276_292_WT_BBR, d276_292_E276F_BBR, equal_var = False) #Welch's t-test b/w WT and L195A AD
st, p37B = stats.ttest_ind(d276_292_WT_BBR, d276_292_K279M_BBR, equal_var = False) #Welch's t-test b/w WT and L195F AD
st, p38B = stats.ttest_ind(d276_292_WT_BBR, d276_292_K279W_BBR, equal_var = False) #Welch's t-test b/w WT and L195N AD
st, p39B = stats.ttest_ind(d276_292_WT_BBR, d276_292_V287T_BBR, equal_var = False) #Welch's t-test b/w WT and L195N AD

plot_func(d200_282, d200_282_err, 'a3', 'a6', 'F196A', p1, p1B, 4, 5)
plot_func(d200_282, d200_282_err, 'a3', 'a6', 'L192A', p2, p2B, 6, 7)
plot_func(d200_282, d200_282_err, 'a3', 'a6', 'L192F', p3, p3B, 8, 9)
plot_func(d200_282, d200_282_err, 'a3', 'a6', 'L195A', p4, p4B, 10, 11)
plot_func(d200_282, d200_282_err, 'a3', 'a6', 'L195F', p5, p5B, 12, 13)
plot_func(d200_282, d200_282_err, 'a3', 'a6', 'L195N', p6, p6B, 14, 15)

plot_func(d200_287, d200_287_err, 'a3', 'a7', 'F196A', p7, p7B, 4, 5)
plot_func(d200_287, d200_287_err, 'a3', 'a7', 'L192A', p8, p8B, 6, 7)
plot_func(d200_287, d200_287_err, 'a3', 'a7', 'L192F', p9, p9B, 8, 9)
plot_func(d200_287, d200_287_err, 'a3', 'a7', 'L195A', p10, p10B, 10, 11)
plot_func(d200_287, d200_287_err, 'a3', 'a7', 'L195F', p11, p11B, 12, 13)
plot_func(d200_287, d200_287_err, 'a3', 'a7', 'L195N', p12, p12B, 14, 15)

plot_func(d276_292, d276_292_err, 'a6', 'a7', 'F196A', p13, p13B, 4, 5)
plot_func(d276_292, d276_292_err, 'a6', 'a7', 'L192A', p14, p14B, 6, 7)
plot_func(d276_292, d276_292_err, 'a6', 'a7', 'L192F', p15, p15B, 8, 9)
plot_func(d276_292, d276_292_err, 'a6', 'a7', 'L195A', p16, p16B, 10, 11)
plot_func(d276_292, d276_292_err, 'a6', 'a7', 'L195F', p17, p17B, 12, 13)
plot_func(d276_292, d276_292_err, 'a6', 'a7', 'L195N', p18, p18B, 14, 15)

plot_func(d200_282, d200_282_err, 'a3', 'a6', 'S286A', p19, p19B, 16, 17)
plot_func(d200_282, d200_282_err, 'a3', 'a6', 'F280Y', p20, p20B, 18, 19)
plot_func(d200_282, d200_282_err, 'a3', 'a6', 'E276L', p21, p21B, 20, 21)
plot_func(d200_282, d200_282_err, 'a3', 'a6', 'E276F', p22, p22B, 22, 23)
plot_func(d200_282, d200_282_err, 'a3', 'a6', 'K279M', p23, p23B, 24, 25)
plot_func(d200_282, d200_282_err, 'a3', 'a6', 'K279W', p24, p24B, 26, 27)
plot_func(d200_282, d200_282_err, 'a3', 'a6', 'V287T', p25, p25B, 28, 29)

plot_func(d200_287, d200_287_err, 'a3', 'a7', 'S286A', p26, p26B, 16, 17)
plot_func(d200_287, d200_287_err, 'a3', 'a7', 'F280Y', p27, p27B, 18, 19)
plot_func(d200_287, d200_287_err, 'a3', 'a7', 'E276L', p28, p28B, 20, 21)
plot_func(d200_287, d200_287_err, 'a3', 'a7', 'E276F', p29, p29B, 22, 23)
plot_func(d200_287, d200_287_err, 'a3', 'a7', 'K279M', p30, p30B, 24, 25)
plot_func(d200_287, d200_287_err, 'a3', 'a7', 'K279W', p31, p31B, 26, 27)
plot_func(d200_287, d200_287_err, 'a3', 'a7', 'V287T', p32, p32B, 28, 29)

plot_func(d276_292, d276_292_err, 'a6', 'a7', 'S286A', p33, p33B, 16, 17)
plot_func(d276_292, d276_292_err, 'a6', 'a7', 'F280Y', p34, p34B, 18, 19)
plot_func(d276_292, d276_292_err, 'a6', 'a7', 'E276L', p35, p35B, 20, 21)
plot_func(d276_292, d276_292_err, 'a6', 'a7', 'E276F', p36, p36B, 22, 23)
plot_func(d276_292, d276_292_err, 'a6', 'a7', 'K279M', p37, p37B, 24, 25)
plot_func(d276_292, d276_292_err, 'a6', 'a7', 'K279W', p38, p38B, 26, 27)
plot_func(d276_292, d276_292_err, 'a6', 'a7', 'V287T', p39, p39B, 28, 29)

output = open('a3_a6.txt', 'w')
output.write('Apo Open:' + str(d200_282[0]) + '+/-' + str(d200_282_err[0]) + '\n')
output.write('Apo Closed:' + str(d200_282[1]) + '+/-' + str(d200_282_err[1]) + '\n')
output.write('WT AD:' + str(d200_282[2]) + '+/-'  + str(d200_282_err[2]) + '\n')
output.write('WT BBR:' + str(d200_282[3]) + '+/-'  + str(d200_282_err[3]) + '\n')
output.write('F196A AD:' + str(d200_282[4]) + '+/-'  + str(d200_282_err[4]) + '\n')
output.write('F196A BBR:' + str(d200_282[5]) + '+/-'  + str(d200_282_err[5]) + '\n')
output.write('L192A AD:' + str(d200_282[6]) + '+/-'  + str(d200_282_err[6]) + '\n')
output.write('L192A BBR:' + str(d200_282[7]) + '+/-'  + str(d200_282_err[7]) + '\n')
output.write('L192F AD:' + str(d200_282[8]) + '+/-'  + str(d200_282_err[8]) + '\n')
output.write('L192F BBR:' + str(d200_282[9]) + '+/-'  + str(d200_282_err[9]) + '\n')
output.write('L195A AD:' + str(d200_282[10]) + '+/-'  + str(d200_282_err[10]) + '\n')
output.write('L195A BBR:' + str(d200_282[11]) + '+/-'  + str(d200_282_err[11]) + '\n')
output.write('L195F AD:' + str(d200_282[12]) + '+/-'  + str(d200_282_err[12]) + '\n')
output.write('L195F BBR:' + str(d200_282[13]) + '+/-'  + str(d200_282_err[13]) + '\n')
output.write('L195N AD:' + str(d200_282[14]) + '+/-'  + str(d200_282_err[14]) + '\n')
output.write('L195N BBR:' + str(d200_282[15]) + '+/-'  + str(d200_282_err[15]) + '\n')
output.write('S286A AD:' + str(d200_282[16]) + '+/-'  + str(d200_282_err[16]) + '\n')
output.write('S286A BBR:' + str(d200_282[17]) + '+/-'  + str(d200_282_err[17]) + '\n')
output.write('F280Y AD:' + str(d200_282[18]) + '+/-'  + str(d200_282_err[18]) + '\n')
output.write('F280Y BBR:' + str(d200_282[19]) + '+/-'  + str(d200_282_err[19]) + '\n')
output.write('E276L AD:' + str(d200_282[20]) + '+/-'  + str(d200_282_err[20]) + '\n')
output.write('E276L BBR:' + str(d200_282[21]) + '+/-'  + str(d200_282_err[21]) + '\n')
output.write('E276F AD:' + str(d200_282[22]) + '+/-'  + str(d200_282_err[22]) + '\n')
output.write('E276F BBR:' + str(d200_282[23]) + '+/-'  + str(d200_282_err[23]) + '\n')
output.write('K279M AD:' + str(d200_282[24]) + '+/-'  + str(d200_282_err[24]) + '\n')
output.write('K279M BBR:' + str(d200_282[25]) + '+/-'  + str(d200_282_err[25]) + '\n')
output.write('K279W AD:' + str(d200_282[26]) + '+/-'  + str(d200_282_err[26]) + '\n')
output.write('K279W BBR:' + str(d200_282[27]) + '+/-'  + str(d200_282_err[27]) + '\n')
output.write('V287T AD:' + str(d200_282[28]) + '+/-'  + str(d200_282_err[28]) + '\n')
output.write('V287T BBR:' + str(d200_282[29]) + '+/-'  + str(d200_282_err[29]) + '\n')

output2 = open('a3_a7.txt', 'w')
output2.write('Apo Open:' + str(d200_287[0]) + '+/-' + str(d200_287_err[0]) + '\n')
output2.write('Apo Closed:' + str(d200_287[1]) + '+/-' + str(d200_287_err[1]) + '\n')
output2.write('WT AD:' + str(d200_287[2]) + '+/-' + str(d200_287_err[2]) + '\n')
output2.write('WT BBR:' + str(d200_287[3]) + '+/-' + str(d200_287_err[3]) + '\n')
output2.write('F196A AD:' + str(d200_287[4]) + '+/-' + str(d200_287_err[4]) + '\n')
output2.write('F196A BBR:' + str(d200_287[5]) + '+/-' + str(d200_287_err[5]) + '\n')
output2.write('L192A AD:' + str(d200_287[6]) + '+/-' + str(d200_287_err[6]) + '\n')
output2.write('L192A BBR:' + str(d200_287[7]) + '+/-' + str(d200_287_err[7]) + '\n')
output2.write('L192F AD:' + str(d200_287[8]) + '+/-' + str(d200_287_err[8]) + '\n')
output2.write('L192F BBR:' + str(d200_287[9]) + '+/-' + str(d200_287_err[9]) + '\n')
output2.write('L195A AD:' + str(d200_287[10]) + '+/-' + str(d200_287_err[10]) + '\n')
output2.write('L195A BBR:' + str(d200_287[11]) + '+/-' + str(d200_287_err[11]) + '\n')
output2.write('L195F AD:' + str(d200_287[12]) + '+/-' + str(d200_287_err[12]) + '\n')
output2.write('L195F BBR:' + str(d200_287[13]) + '+/-' + str(d200_287_err[13]) + '\n')
output2.write('L195N AD:' + str(d200_287[14]) + '+/-' + str(d200_287_err[14]) + '\n')
output2.write('L195N BBR:' + str(d200_287[15]) + '+/-' + str(d200_287_err[15]))
output2.write('S286A AD:' + str(d200_287[16]) + '+/-'  + str(d200_287_err[16]) + '\n')
output2.write('S286A BBR:' + str(d200_287[17]) + '+/-'  + str(d200_287_err[17]) + '\n')
output2.write('F280Y AD:' + str(d200_287[18]) + '+/-'  + str(d200_287_err[18]) + '\n')
output2.write('F280Y BBR:' + str(d200_287[19]) + '+/-'  + str(d200_287_err[19]) + '\n')
output2.write('E276L AD:' + str(d200_287[20]) + '+/-'  + str(d200_287_err[20]) + '\n')
output2.write('E276L BBR:' + str(d200_287[21]) + '+/-'  + str(d200_287_err[21]) + '\n')
output2.write('E276F AD:' + str(d200_287[22]) + '+/-'  + str(d200_287_err[22]) + '\n')
output2.write('E276F BBR:' + str(d200_287[23]) + '+/-'  + str(d200_287_err[23]) + '\n')
output2.write('K279M AD:' + str(d200_287[24]) + '+/-'  + str(d200_287_err[24]) + '\n')
output2.write('K279M BBR:' + str(d200_287[25]) + '+/-'  + str(d200_287_err[25]) + '\n')
output2.write('K279W AD:' + str(d200_287[26]) + '+/-'  + str(d200_287_err[26]) + '\n')
output2.write('K279W BBR:' + str(d200_287[27]) + '+/-'  + str(d200_287_err[27]) + '\n')
output2.write('V287T AD:' + str(d200_287[28]) + '+/-'  + str(d200_287_err[28]) + '\n')
output2.write('V287T BBR:' + str(d200_287[29]) + '+/-'  + str(d200_287_err[29]) + '\n')

output3 = open('a6_a7.txt', 'w')
output3.write('Apo Open:' + str(d276_292[0]) + '+/-' + str(d276_292_err[0]) + '\n')
output3.write('Apo Closed:' + str(d276_292[1]) + '+/-' + str(d276_292_err[1]) + '\n')
output3.write('WT AD:' + str(d276_292[2]) + '+/-' + str(d276_292_err[2]) + '\n')
output3.write('WT BBR:' + str(d276_292[3]) + '+/-' + str(d276_292_err[3]) + '\n')
output3.write('F196A AD:' + str(d276_292[4]) + '+/-' + str(d276_292_err[4]) + '\n')
output3.write('F196A BBR:' + str(d276_292[5]) + '+/-' + str(d276_292_err[5]) + '\n')
output3.write('L192A AD:' + str(d276_292[6]) + '+/-' + str(d276_292_err[6]) + '\n')
output3.write('L192A BBR:' + str(d276_292[7]) + '+/-' + str(d276_292_err[7]) + '\n')
output3.write('L192F AD:' + str(d276_292[8]) + '+/-' + str(d276_292_err[8]) + '\n')
output3.write('L192F BBR:' + str(d276_292[9]) + '+/-' + str(d276_292_err[9]) + '\n')
output3.write('L195A AD:' + str(d276_292[10]) + '+/-' + str(d276_292_err[10]) + '\n')
output3.write('L195A BBR:' + str(d276_292[11]) + '+/-' + str(d276_292_err[11]) + '\n')
output3.write('L195F AD:' + str(d276_292[12]) + '+/-' + str(d276_292_err[12]) + '\n')
output3.write('L195F BBR:' + str(d276_292[13]) + '+/-' + str(d276_292_err[13]) + '\n')
output3.write('L195N AD:' + str(d276_292[14]) + '+/-' + str(d276_292_err[14]) + '\n')
output3.write('L195N BBR:' + str(d276_292[15]) + '+/-' + str(d276_292_err[15]))
output3.write('S286A AD:' + str(d276_292[16]) + '+/-'  + str(d276_292_err[16]) + '\n')
output3.write('S286A BBR:' + str(d276_292[17]) + '+/-'  + str(d276_292_err[17]) + '\n')
output3.write('F280Y AD:' + str(d276_292[18]) + '+/-'  + str(d276_292_err[18]) + '\n')
output3.write('F280Y BBR:' + str(d276_292[19]) + '+/-'  + str(d276_292_err[19]) + '\n')
output3.write('E276L AD:' + str(d276_292[20]) + '+/-'  + str(d276_292_err[20]) + '\n')
output3.write('E276L BBR:' + str(d276_292[21]) + '+/-'  + str(d276_292_err[21]) + '\n')
output3.write('E276F AD:' + str(d276_292[22]) + '+/-'  + str(d276_292_err[22]) + '\n')
output3.write('E276F BBR:' + str(d276_292[23]) + '+/-'  + str(d276_292_err[23]) + '\n')
output3.write('K279M AD:' + str(d276_292[24]) + '+/-'  + str(d276_292_err[24]) + '\n')
output3.write('K279M BBR:' + str(d276_292[25]) + '+/-'  + str(d276_292_err[25]) + '\n')
output3.write('K279W AD:' + str(d276_292[26]) + '+/-'  + str(d276_292_err[26]) + '\n')
output3.write('K279W BBR:' + str(d276_292[27]) + '+/-'  + str(d276_292_err[27]) + '\n')
output3.write('V287T AD:' + str(d276_292[28]) + '+/-'  + str(d276_292_err[28]) + '\n')
output3.write('V287T BBR:' + str(d276_292[29]) + '+/-'  + str(d276_292_err[29]) + '\n')

