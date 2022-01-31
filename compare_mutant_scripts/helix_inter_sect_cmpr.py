#Import Necessary Packages
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from itertools import product

#Load data, determine correlated samples and caculate mean and error
def load_data(mut, lig):
    a3_a7_pt1, a3_a7_pt2, a6_a7_pt1, a6_a7_pt2, a6_a7_pt3 = [],[],[],[],[]
    if lig == 'AD' or lig == 'BBR':
        for i in open('../../../' + mut + '/' + lig + '/a3_a7_pt1_tot_inter.txt'):
            a3_a7_pt1.append(float(i))
        for i in open('../../../' + mut + '/' + lig + '/a3_a7_pt2_tot_inter.txt'):
            a3_a7_pt2.append(float(i))
        for i in open('../../../' + mut + '/' + lig + '/a6_a7_pt1_tot_inter.txt'):
            a6_a7_pt1.append(float(i))
        for i in open('../../../' + mut + '/' + lig + '/a6_a7_pt2_tot_inter.txt'):
            a6_a7_pt2.append(float(i))
        for i in open('../../../' + mut + '/' + lig + '/a6_a7_pt3_tot_inter.txt'):
            a6_a7_pt3.append(float(i))

    elif lig == 'ApoO':
        for i in open('~/doc/PTP1B/Apo_dis/config11/a3_a7_pt1_tot_inter.txt'):
            a3_a7_pt1.append(float(i))
        for i in open('~/doc/PTP1B/Apo_dis/config11/a3_a7_pt2_tot_inter.txt'):
            a3_a7_pt2.append(float(i))
        for i in open('~/doc/PTP1B/Apo_dis/config11/a6_a7_pt1_tot_inter.txt'):
            a6_a7_pt1.append(float(i))
        for i in open('~/doc/PTP1B/Apo_dis/config11/a6_a7_pt2_tot_inter.txt'):
            a6_a7_pt2.append(float(i))
        for i in open('~/doc/PTP1B/Apo_dis/config11/a6_a7_pt3_tot_inter.txt'):
            a6_a7_pt3.append(float(i))

    elif lig == 'ApoC':
        for i in open('~/doc/PTP1B/1sug/a3_a7_pt1_tot_inter.txt'):
            a3_a7_pt1.append(float(i))
        for i in open('~/doc/PTP1B/1sug/a3_a7_pt2_tot_inter.txt'):
            a3_a7_pt2.append(float(i))
        for i in open('~/doc/PTP1B/1sug/a6_a7_pt1_tot_inter.txt'):
            a6_a7_pt1.append(float(i))
        for i in open('~/doc/PTP1B/1sug/a6_a7_pt2_tot_inter.txt'):
            a6_a7_pt2.append(float(i))
        for i in open('~/doc/PTP1B/1sug/a6_a7_pt3_tot_inter.txt'):
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
    
    #Caculate Mean and standard error
    mean_data = [np.mean(uncorr_data[0,:]), np.mean(uncorr_data[1,:]), np.mean(uncorr_data[2,:]), np.mean(uncorr_data[3,:]), np.mean(uncorr_data[4,:])]
    err_data = [stats.sem(uncorr_data[0,:]), stats.sem(uncorr_data[1,:]), stats.sem(uncorr_data[2,:]), stats.sem(uncorr_data[3,:]), stats.sem(uncorr_data[4,:])]

    #Seperate arrays for different measurements
    a3_a7_pt1 = uncorr_data[0,:]
    a3_a7_pt2 = uncorr_data[1,:]
    a6_a7_pt1 = uncorr_data[2,:]
    a6_a7_pt2 = uncorr_data[3,:]
    a6_a7_pt3 = uncorr_data[4,:]

    return mean_data, err_data, a3_a7_pt1, a3_a7_pt2, a6_a7_pt1, a6_a7_pt2, a6_a7_pt3

#Plot
def plot_indv_mut(inter_AD, inter_BBR, err_AD, err_BBR, label, mut, P_AD, P_BBR):
    for j in range(len(label)):
        for k in range(1, len(mut)):
            #Get data into proper format
            inter = [inter_AD[0][j], inter_BBR[0][j], inter_AD[k][j], inter_BBR[k][j]]
            err = [err_AD[0][j], err_BBR[0][j], err_AD[k][j], err_BBR[k][j]]
            
            #Set Plot Colors
            Color = ['blue', 'blue']
            if P_AD[k, j] < 0.01:
                if inter[2] - inter[0] < 0:
                    Color.append('green')
                if inter[2] - inter[0] > 0:
                    Color.append('red')
            else:
                Color.append('blue')
            if P_BBR[k, j] < 0.01:
                if inter[3] - inter[1] < 0:
                    Color.append('green')
                if inter[3] - inter[1] > 0:
                    Color.append('red')
            else:
                Color.append('blue')

            #Plot Function
            num = [5, 10, 15, 20]
            Method = ['WT AD', 'WT BBR', mut[k] + ' AD', mut[k] + ' BBR']
            fig = plt.figure()
            ax1 = fig.add_subplot(111)
            ax1.set_title('Comparison of' + label[j] + 'Helix Interactions ' + mut[k])    
            ax1.set_ylabel('Mean Number of Interactions')
            ax1.bar(num, inter, color = Color, width=4.5)
            plt.errorbar(num, inter, yerr= err, fmt='o', color='black')
            plt.xticks(num, Method, fontsize=8)
            fig.savefig('./' + mut[k] + '/' + label[j] + '_' + mut[k] + '_inter.png')
            plt.close(fig)

def plot_mult_mut(inter, err, mut, label, P, lig, which):
    for j in range(len(label)):
        #Process data
        num_mut = len(mut)
        num = np.linspace(0, num_mut*2, num = num_mut)
        Method = mut
        inter_j = inter[:,j]
        err_j = err[:,j]
        
        #Set Plot Colors
        Color = ['black', 'gray', 'blue']
        for k in range(3, num_mut):
            if P[k, j] < 0.01:
                if inter_j[k] - inter_j[0] < 0:
                    Color.append('green')
                if inter_j[k] - inter_j[0] > 0:
                    Color.append('red')
            else:
                Color.append('blue')

        #Plot Function
        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        ax1.set_title('Comparison of' + label[j] + 'Helix Interactions for ' + lig + ' ' + which)
        ax1.set_ylabel('Mean Number of Interactions')
        ax1.bar(num, inter_j, color = Color, width=1.8)
        plt.errorbar(num, inter_j, yerr= err_j, fmt='o', color='black')
        plt.xticks(num, Method, fontsize=8)
        fig.savefig(label[j] + '_' + lig + '_' + which + '_inter.png')
        plt.close(fig)
        
        if which == 'all':
            #Write to file
            output = open(label[j] + '_' + lig + '_' + which + '_inter.txt', 'w')
            for k in range(len(mut)):
                output.write(str(mut[k]) + ': ' + str(inter[k,j]) + ' +/- ' + str(err[k,j]) + '\n')
                output.write('P-value:' + str(P[k,j]) + '\n')
            output.close()

#Load Data for all mutants and WT
mut = ['Apo Open', 'Apo Closed', 'WT', 'F196A', 'L192A', 'L192F', 'L195A', 'L195F', 'L195N', 'S286A', 'F280Y', 'E276L', 'E276F', 'K279M', 'K279W', 'V287T']
inter = ['a3-a7-pt1', 'a3-a7-pt2', 'a6-a7-pt1', 'a6-a7-pt2', 'a6-a7-pt3']
all_mean_AD = np.zeros((len(mut), len(inter)))
all_mean_BBR = np.zeros((len(mut), len(inter)))
all_err_AD = np.zeros((len(mut), len(inter)))
all_err_BBR = np.zeros((len(mut), len(inter)))

#Load Data for Apo States
all_mean_AD[0][:], all_err_AD[0][:], a3_a7_pt1_ApoO, a3_a7_pt2_ApoO, a6_a7_pt1_ApoO, a6_a7_pt2_ApoO, a6_a7_pt3_ApoO  = load_data(mut[0], 'ApoO')
all_mean_AD[1][:], all_err_AD[1][:], a3_a7_pt1_ApoC, a3_a7_pt2_ApoC, a6_a7_pt1_ApoC, a6_a7_pt2_ApoC, a6_a7_pt3_ApoC  = load_data(mut[1], 'ApoC')
all_mean_BBR[0][:] = all_mean_AD[0][:]
all_mean_BBR[0][:] = all_mean_AD[0][:]

#Load data for WT
all_mean_AD[2][:], all_err_AD[2][:], a3_a7_pt1_ADWT, a3_a7_pt2_ADWT, a6_a7_pt1_ADWT, a6_a7_pt2_ADWT, a6_a7_pt3_ADWT  = load_data(mut[2], 'AD')
all_mean_BBR[2][:], all_err_BBR[2][:], a3_a7_pt1_BBRWT, a3_a7_pt2_BBRWT, a6_a7_pt1_BBRWT, a6_a7_pt2_BBRWT, a6_a7_pt3_BBR_WT = load_data(mut[2], 'BBR')

#empty array for p values
P_AD = np.zeros((len(mut), len(inter)))
P_BBR = np.zeros((len(mut), len(inter)))

#Load data for mutants
for i in range(3, len(mut), 1):
    all_mean_AD[i][:], all_err_AD[i][:], a3_a7_pt1_ADmut, a3_a7_pt2_ADmut, a6_a7_pt1_ADmut, a6_a7_pt2_ADmut, a6_a7_pt3_ADmut  = load_data(mut[i], 'AD')
    all_mean_BBR[i][:], all_err_BBR[i][:], a3_a7_pt1_BBRmut, a3_a7_pt2_BBRmut, a6_a7_pt1_BBRmut, a6_a7_pt2_BBRmut, a6_a7_pt3_BBRmut = load_data(mut[i], 'BBR')

    #Calculate all p values
    st, p1 = stats.ttest_ind(a3_a7_pt1_ADWT, a3_a7_pt1_ADmut, equal_var = False) #Welch's t-test b/w WT and mutant
    st, p2 = stats.ttest_ind(a3_a7_pt2_ADWT, a3_a7_pt2_ADmut, equal_var = False) #Welch's t-test b/w WT and mutant
    st, p3 = stats.ttest_ind(a6_a7_pt1_ADWT, a6_a7_pt1_ADmut, equal_var = False) #Welch's t-test b/w WT and mutant
    st, p4 = stats.ttest_ind(a6_a7_pt2_ADWT, a6_a7_pt2_ADmut, equal_var = False) #Welch's t-test b/w WT and mutant
    st, p5 = stats.ttest_ind(a6_a7_pt3_ADWT, a6_a7_pt3_ADmut, equal_var = False) #Welch's t-test b/w WT and mutant
   
    st, p6 = stats.ttest_ind(a3_a7_pt1_BBRWT, a3_a7_pt1_BBRmut, equal_var = False) #Welch's t-test b/w WT and mutant
    st, p7 = stats.ttest_ind(a3_a7_pt2_BBRWT, a3_a7_pt2_BBRmut, equal_var = False) #Welch's t-test b/w WT and mutant
    st, p8 = stats.ttest_ind(a6_a7_pt1_BBRWT, a6_a7_pt1_BBRmut, equal_var = False) #Welch's t-test b/w WT and mutant
    st, p9 = stats.ttest_ind(a6_a7_pt2_BBRWT, a6_a7_pt2_BBRmut, equal_var = False) #Welch's t-test b/w WT and mutant
    st, p10 = stats.ttest_ind(a6_a7_pt3_BBRWT, a6_a7_pt3_BBRmut, equal_var = False) #Welch's t-test b/w WT and mutant
  
    #Save to array
    P_AD[i, :] = [p1, p2, p3, p4, p5]
    P_BBR[i, :] = [p6, p7, p8, p9, p10]

#All individual mutants
plot_indv_mut(all_mean_AD, all_mean_BBR, all_err_AD, all_err_BBR, inter, mut, P_AD, P_BBR)

#All mutants together
plot_mult_mut(all_mean_AD, all_err_AD, mut, inter, P_AD, 'AD', 'all')
plot_mult_mut(all_mean_BBR, all_err_BBR, mut, inter, P_BBR, 'BBR', 'all')

#Only mutants with similar binding affinity
sim_AD = [0, 1, 2, 4]
sim_mean_AD = np.zeros((len(sim_AD), len(inter)))
sim_err_AD = np.zeros((len(sim_AD), len(inter)))
sim_P_AD = np.zeros((len(sim_AD), len(inter)))
sim_mut_AD = []

sim_BBR = [0, 1, 2, 4]
sim_mean_BBR = np.zeros((len(sim_BBR), len(inter)))
sim_err_BBR = np.zeros((len(sim_BBR), len(inter)))
sim_P_BBR = np.zeros((len(sim_BBR), len(inter)))
sim_mut_BBR = []

n = 0
m = 0
for i in range(len(mut)):
    if i in sim_AD:
        sim_mut_AD.append(mut[i])
        sim_mean_AD[n][:] = all_mean_AD[i][:]
        sim_err_AD[n][:] = all_err_AD[i][:]
        sim_P_AD[n][:] = P_AD[i][:]
        n += 1
    if i in sim_BBR:
        sim_mut_BBR.append(mut[i])
        sim_mean_BBR[m][:] = all_mean_BBR[i][:]
        sim_err_BBR[m][:] = all_err_BBR[i][:]
        sim_P_BBR[m][:] = P_BBR[i][:]
        m += 1

plot_mult_mut(sim_mean_AD, sim_err_AD, sim_mut_AD, inter, sim_P_AD, 'AD', 'sim_bind')
plot_mult_mut(sim_mean_BBR, sim_err_BBR, sim_mut_BBR, inter, sim_P_BBR, 'BBR', 'sim_bind')

#Determine percent difference from WT
mut_only = ['F196A', 'L192A', 'L192F', 'L195A', 'L195F', 'L195N', 'S286A', 'F280Y', 'E276L', 'E276F', 'K279M', 'K279W', 'V287T']
per_diff_AD = np.zeros(len(mut_only))
per_diff_BBR = np.zeros(len(mut_only))
for i in range(len(mut_only)):
    per_diff_AD[i] = (all_mean_AD[i+3] - all_mean_AD[2])/((all_mean_AD[i+3] + all_mean_AD[2])/2) * 100
    per_diff_BBR[i] = (all_mean_BBR[i+3] - all_mean_BBR[2])/((all_mean_BBR[i+3] + all_mean_BBR[2])/2) * 100
   
#Plot table comparing residue interactions to WT
ax = plt.figure(figsize=(8, 10), frameon=False) # no visible frame
ax = sns.heatmap(per_diff_AD, annot=False, cmap = 'RdYlBu', xticklabels = mut_only, yticklabels = inter)
#ax.add_artist(lines.Line2D([0, 20], [7, 7], color = 'black', linestyle= '--', linewidth = 4))
plt.title('Helical Distance Compared to WT for AD')
plt.savefig('mutate_AD_helix_dist.png')
plt.close()

ax = plt.figure(figsize=(8, 10), frameon=False) # no visible frame
ax = sns.heatmap(per_diff_BBR, annot=False, cmap = 'RdYlBu', xticklabels = mut_only, yticklabels = inter)
#ax.add_artist(lines.Line2D([0, 20], [7, 7], color = 'black', linestyle= '--', linewidth = 4))
plt.title('Helical Distance Compared to WT for BBR')
plt.savefig('mutate_AD_helix_dist.png')
plt.close()

