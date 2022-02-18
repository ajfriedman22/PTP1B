#Import Necessary Packages
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from itertools import product
import seaborn as sns
#Load data, determine correlated samples and caculate mean and error
def load_data(Path):
    a3_a7_pt1, a3_a7_pt2, a6_a7_pt1, a6_a7_pt2, a6_a7_pt3, a7_L11 = [],[],[],[],[],[]
    for i in open('../../' + Path + '/a3_a7_pt1_tot_inter.txt'):
        a3_a7_pt1.append(float(i))
    for i in open('../../' + Path +'/a3_a7_pt2_tot_inter.txt'):
        a3_a7_pt2.append(float(i))
    for i in open('../../' + Path + '/a6_a7_pt1_tot_inter.txt'):
        a6_a7_pt1.append(float(i))
    for i in open('../../' + Path + '/a6_a7_pt2_tot_inter.txt'):
        a6_a7_pt2.append(float(i))
    for i in open('../../' + Path + '/a6_a7_pt3_tot_inter.txt'):
        a6_a7_pt3.append(float(i))
    for i in open('../../' + Path + '/a7_L11_inter.txt'):
        a7_L11.append(float(i))


    tot_data = [a3_a7_pt1, a3_a7_pt2, a6_a7_pt1, a6_a7_pt2, a6_a7_pt3, a7_L11]
    
    #Determine correlation time
    corr_time = []
    for i in range(6):
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
    uncorr_data = np.zeros((6, num_ucorr))
    n = 0
    for i in range(len(tot_data[0][:])):
        if i % mean_corr_time == 0:
            uncorr_data[0,n] = tot_data[0][i]
            uncorr_data[1,n] = tot_data[1][i]
            uncorr_data[2,n] = tot_data[2][i]
            uncorr_data[3,n] = tot_data[3][i]
            uncorr_data[4,n] = tot_data[4][i]
            uncorr_data[5,n] = tot_data[5][i]
            n += 1
    
    #Caculate Mean and standard error
    mean_data = [np.mean(uncorr_data[0,:]), np.mean(uncorr_data[1,:]), np.mean(uncorr_data[2,:]), np.mean(uncorr_data[3,:]), np.mean(uncorr_data[4,:]), np.mean(uncorr_data[5,:])]
    err_data = [stats.sem(uncorr_data[0,:]), stats.sem(uncorr_data[1,:]), stats.sem(uncorr_data[2,:]), stats.sem(uncorr_data[3,:]), stats.sem(uncorr_data[4,:]), stats.sem(uncorr_data[5,:])]

    #Seperate arrays for different measurements
    a3_a7_pt1 = uncorr_data[0,:]
    a3_a7_pt2 = uncorr_data[1,:]
    a6_a7_pt1 = uncorr_data[2,:]
    a6_a7_pt2 = uncorr_data[3,:]
    a6_a7_pt3 = uncorr_data[4,:]
    a7_L11 = uncorr_data[5,:]

    return mean_data, err_data, a3_a7_pt1, a3_a7_pt2, a6_a7_pt1, a6_a7_pt2, a6_a7_pt3, a7_L11

#Plot
def plot_indv_mut(inter_all, err_all, label, mut, P):
    for j in range(len(label)):
        for k in range(2, len(mut)):
            #Get data into proper format
            inter = [inter_all[0][j], inter_all[1][j], inter_all[k][j]]
            err = [err_all[0][j], err_all[1][j], err_all[k][j]]
            
            #Set Plot Colors
            Color = ['blue', 'blue']
            if P[k, j] < 0.01:
                if inter[2] - inter[0] < 0:
                    Color.append('green')
                if inter[2] - inter[0] > 0:
                    Color.append('red')
            else:
                Color.append('blue')

            #Plot Function
            num = [5, 10, 15]
            Method = ['Apo Open', 'Apo Closed', mut[k]]
            fig = plt.figure()
            ax1 = fig.add_subplot(111)
            ax1.set_title('Comparison of' + label[j] + 'Helix Interactions ' + mut[k])    
            ax1.set_ylabel('Mean Number of Interactions')
            ax1.bar(num, inter, color = Color, width=4.5)
            plt.errorbar(num, inter, yerr= err, fmt='o', color='black')
            plt.xticks(num, Method, fontsize=8)
            fig.savefig('./' + mut[k] + '/' + label[j] + '_' + mut[k] + '_inter.png')
            plt.close(fig)

def plot_mult_mut(inter, err, mut, label, P, which):
    for j in range(len(label)):
        #Process data
        num_mut = len(mut)
        num = np.linspace(0, num_mut*3, num = num_mut)
        Method = mut
        inter_j = inter[:,j]
        err_j = err[:,j]
        
        #Set Plot Colors
        Color = ['black', 'gray']
        for k in range(2, num_mut):
            if P[k, j] < 0.05:
                if inter_j[k] - inter_j[0] < 0:
                    Color.append('green')
                if inter_j[k] - inter_j[0] > 0:
                    Color.append('red')
            else:
                Color.append('lightblue')

        #Plot Function
        fig = plt.figure(figsize=(12, 8))
        ax1 = fig.add_subplot(111)
        ax1.set_title('Comparison of ' + label[j] + ' Helix Interactions for ' + which)
        ax1.set_ylabel('Mean Number of Interactions')
        ax1.bar(num, inter_j, color = Color, width=2.8)
        plt.errorbar(num, inter_j, yerr= err_j, fmt='o', color='black')
        plt.xticks(num, Method, fontsize=8)
        fig.savefig(label[j] + '_' + which + '_inter.png')
        plt.close(fig)
        
        if which == 'all':
            #Write to file
            output = open(label[j] + '_' + which + '_inter.txt', 'w')
            for k in range(2, len(mut)):
                output.write(str(mut[k]) + ': ' + str(inter[k,j]) + ' +/- ' + str(err[k,j]) + '\n')
                output.write('P-value:' + str(P[k,j]) + '\n')
            output.close()

#Load Data for all mutants and WT
mut = ['Apo Open', 'Apo Closed', 'F196A', 'L192F', 'L195F', 'F280Y', 'E276F', 'V287T']
inter = ['a3-a7-pt1', 'a3-a7-pt2', 'a6-a7-pt1', 'a6-a7-pt2', 'a6-a7-pt3', 'a7-L11']
all_mean = np.zeros((len(mut), len(inter)))
all_err = np.zeros((len(mut), len(inter)))

#File Paths
file_path = ['../Apo_dis/analysis', '../Apo_1SUG/analysis/1sug', 'F196A/Apo/analysis',  'L192F/Apo/analysis', 'L195F/Apo/analysis', 'F280Y/Apo/analysis', 'E276F/Apo/analysis', 'V287T/Apo/analysis']


#Load Data for Apo States
all_mean[0][:], all_err[0][:], a3_a7_pt1_ApoO, a3_a7_pt2_ApoO, a6_a7_pt1_ApoO, a6_a7_pt2_ApoO, a6_a7_pt3_ApoO, a7_L11_ApoO  = load_data(file_path[0])
all_mean[1][:], all_err[1][:], a3_a7_pt1_ApoC, a3_a7_pt2_ApoC, a6_a7_pt1_ApoC, a6_a7_pt2_ApoC, a6_a7_pt3_ApoC, a7_L11_ApoC  = load_data(file_path[1])

#empty array for p values
P = np.zeros((len(mut), len(inter)))

#Load data for mutants
for i in range(2, len(mut), 1):
    all_mean[i][:], all_err[i][:], a3_a7_pt1_mut, a3_a7_pt2_mut, a6_a7_pt1_mut, a6_a7_pt2_mut, a6_a7_pt3_mut, a7_L11_mut  = load_data(file_path[i])

    #Calculate all p values
    st, p1 = stats.ttest_ind(a3_a7_pt1_ApoO, a3_a7_pt1_mut, equal_var = False) #Welch's t-test b/w WT and mutant
    st, p2 = stats.ttest_ind(a3_a7_pt2_ApoO, a3_a7_pt2_mut, equal_var = False) #Welch's t-test b/w WT and mutant
    st, p3 = stats.ttest_ind(a6_a7_pt1_ApoO, a6_a7_pt1_mut, equal_var = False) #Welch's t-test b/w WT and mutant
    st, p4 = stats.ttest_ind(a6_a7_pt2_ApoO, a6_a7_pt2_mut, equal_var = False) #Welch's t-test b/w WT and mutant
    st, p5 = stats.ttest_ind(a6_a7_pt3_ApoO, a6_a7_pt3_mut, equal_var = False) #Welch's t-test b/w WT and mutant
    st, p6 = stats.ttest_ind(a7_L11_ApoO, a7_L11_mut, equal_var = False) #Welch's t-test b/w WT and mutant
  
    #Save to array
    P[i, :] = [p1, p2, p3, p4, p5, p6]

#All individual mutants
plot_indv_mut(all_mean, all_err, inter, mut, P)

#All mutants together
plot_mult_mut(all_mean, all_err, mut, inter, P, 'all')

#Determine percent difference from WT
mut_only = ['F196A', 'L192F', 'L195F', 'F280Y', 'E276F', 'V287T']
iters = [2, 3, 4, 5, 6, 7]
per_diff = np.zeros((len(inter), len(mut_only)))
for j in range(len(inter)):
    for i in range(len(mut_only)):
        n = iters[i]
        per_diff[j][i] = (all_mean[n][j] - all_mean[0][j])/((all_mean[n][j] + all_mean[0][j])/2) * 100
 
#Plot table comparing residue interactions to WT
ax = plt.figure(figsize=(12, 6), frameon=False) # no visible frame
ax = sns.heatmap(per_diff, annot=False, cmap = 'bwr', cbar = True, cbar_kws={'label': 'Percentage Difference from WT'}, vmin = -200, vmax = 200, xticklabels = mut_only, yticklabels = inter)
#ax.add_artist(lines.Line2D([0, 20], [7, 7], color = 'black', linestyle= '--', linewidth = 4))
plt.title('Mean Number of Interactions Compared to WT')
plt.savefig('mutate_Apo_helix_dist.png')
plt.close()

