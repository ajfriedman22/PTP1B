#Determine the difference in interactions b/w helicex in Apo as well as ligand bound states
from matplotlib import pyplot as plt
import numpy as np
from scipy import stats
from itertools import product
import seaborn as sns
import pandas as pd
def load_file(file_dir, inters, num_bonds):
    inter_per = np.zeros(num_bonds)
    n = 0
    for i in inters:
        data = open('../../../' + file_dir + '/' + i + '_inter_all.txt', 'r').readlines()
        for j in range(len(data)):
            inter_per[n] = float(data[j])
            n += 1

    return inter_per

def box_plot(Apo_open, Apo_close, AD, BBR, pair, output_dir, p, p1):
    Apo_open_df = pd.DataFrame({'Apo Open': Apo_open})
    Apo_close_df = pd.DataFrame({'Apo Closed': Apo_close})
    AD_df = pd.DataFrame({'AD': AD})
    BBR_df = pd.DataFrame({'BBR': BBR})
    mean = np.array([np.mean(Apo_open), np.mean(Apo_close), np.mean(AD), np.mean(BBR)])

    df = pd.concat([Apo_open_df, Apo_close_df, AD_df, BBR_df])

    ax = sns.stripplot(data = df, dodge=True, alpha=0.25, zorder=1, palette='bright')
    ax = sns.pointplot(data = df, join=False, scale=0.75, palette='dark')
    
    if p < 0.05 and p > 0.01:
        x1, x2 = 0, 2 #Columns for Apo and AD
        y, h, col = (1.1*mean[[0, 2]].max()), 1, 'b'
        plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
        plt.text((x1+x2)*0.5, y+h, "*" , ha='center', va='bottom', color=col)
    if p < 0.01 and p > 0.001:
        x1, x2 = 0, 2 #Columns for Apo and AD
        y, h, col = (1.1*mean[[0, 2]].max()), 1, 'b'
        plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
        plt.text((x1+x2)*0.5, y+h, "**" , ha='center', va='bottom', color=col)
    if p < 0.001:
        x1, x2 = 0, 2 #Columns for Apo and AD
        y, h, col = (1.1*mean[[0, 2]].max()), 1, 'b'
        plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
        plt.text((x1+x2)*0.5, y+h, "***" , ha='center', va='bottom', color=col)
    if p1 < 0.05 and p1 > 0.01:
        x1, x2 = 0, 3 #Columns for Apo and BBR
        y, h, col = (1.1*mean[[0, 3]].max()), 1, 'r'
        plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
        plt.text((x1+x2)*0.5, y+h, "*" , ha='center', va='bottom', color=col)
    if p1 < 0.01 and p1 > 0.001:
        x1, x2 = 0, 3 #Columns for Apo and BBR
        y, h, col = (1.1*mean[[0, 3]].max()), 1, 'r'
        plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
        plt.text((x1+x2)*0.5, y+h, "**" , ha='center', va='bottom', color=col)
    if p1 < 0.001:
        x1, x2 = 0, 3 #Columns for Apo and BBR
        y, h, col = (1.1*mean[[0, 3]].max()), 1, 'r'
        plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
        plt.text((x1+x2)*0.5, y+h, "***" , ha='center', va='bottom', color=col)

    plt.title('Residue Interactions b/w ' + str(pair))
    plt.savefig(output_dir + '/Hel_inter_' + str(pair) + '_box.png')
    plt.close()

#Set interaction pairs
group_3 = np.linspace(185, 199, num = 15) #residues in the a3 helix
group_6 = np.linspace(263, 280, num = 18) #residues in the a6 helix
group_7 = np.linspace(286, 297, num = 12) #residues in the a7 helix
group_L11 = np.linspace(149, 152, num = 4) #residues in the L11 loop

pair_a7_a3 = list(product(group_3, group_7))
pair_a7_a6 = list(product(group_6, group_7))

#Array with all interactions
pairs_all = pair_a7_a3
pairs_all.extend(pair_a7_a6)

#Open input files or contacts over time
Apo_open_dir_list = ['rebuild_a7_high/config9/analysis', 'rebuild_a7_high/config11/analysis', 'Apo_dis/analysis']
Apo_close_dir_list = ['Apo_1SUG/analysis/1sug', 'Apo_1SUG/analysis/1sug2']
AD_dir_list = ['1sug_dis_AD/analysis/config11', '1sug_dis_AD/analysis/config_alt', '1sug_dis_AD/analysis/config_alt2', 'mutate/WT/AD/analysis']
BBR_dir_list = ['BBR_a7/analysis', 'BBR_dis/analysis/config9', 'mutate/WT/BBR/analysis'] 

#Load all percentages into arrays
inters = ['a3_a7', 'a6_a7']
num_bonds = len(pairs_all)
Apo_open_per_inter = np.zeros([len(Apo_open_dir_list), num_bonds])
Apo_close_per_inter = np.zeros([len(Apo_close_dir_list), num_bonds])
AD_per_inter = np.zeros([len(AD_dir_list), num_bonds])
BBR_per_inter = np.zeros([len(BBR_dir_list), num_bonds])

for i in range(len(Apo_open_dir_list)):
    Apo_open_per_inter[i,:] = load_file(Apo_open_dir_list[i], inters, num_bonds)
    
for i in range(len(Apo_close_dir_list)):
    Apo_close_per_inter[i,:] = load_file(Apo_close_dir_list[i], inters, num_bonds)

for i in range(len(AD_dir_list)):
    AD_per_inter[i,:] = load_file(AD_dir_list[i], inters, num_bonds)

for i in range(len(BBR_dir_list)):
    BBR_per_inter[i,:] = load_file(BBR_dir_list[i], inters, num_bonds)

#Output File for significant bonds
output_lig_disrupt = open('Sign_inter_lig_disrupt.txt', 'w')
output_lig_inc = open('Sign_inter_lig_inc.txt', 'w')
output_Apo_disrupt = open('Sign_inter_Apo_disrupt.txt', 'w')

#Determine the mean and sem for each helix interactions
mean_inter = np.zeros([4, num_bonds])
sem_inter = np.zeros([4, num_bonds])

for n in range(num_bonds):
    Apo_open = Apo_open_per_inter[:,n]
    Apo_close = Apo_close_per_inter[:,n]
    AD = AD_per_inter[:,n]
    BBR = BBR_per_inter[:,n]

    mean_inter[0][n] = np.mean(Apo_open) #Apo Open
    mean_inter[1][n] = np.mean(Apo_close) #Apo Closed
    mean_inter[2][n] = np.mean(AD) #AD
    mean_inter[3][n] = np.mean(BBR) #BBR

    sem_inter[0][n] = stats.sem(Apo_open) #Apo Open
    sem_inter[1][n] = stats.sem(Apo_close) #Apo Closed
    sem_inter[2][n] = stats.sem(AD) #AD
    sem_inter[3][n] = stats.sem(BBR) #BBR
    
    #Run Welch's test on BBR vs AD helix interactions for a3-a7
    st1, p1 = stats.ttest_ind(AD, BBR, equal_var = False) #Welch's t-test b/w interactions with AD open + BBR open
    st2, p2 = stats.ttest_ind(AD, Apo_open, equal_var = False) #Welch's t-test b/w interactions with AD + Apo Open
    st3, p3 = stats.ttest_ind(BBR, Apo_open, equal_var = False) #Welch's t-test b/w interactions with BBR + Apo Open
    st4, p4 = stats.ttest_ind(Apo_close, Apo_open, equal_var = False) #Welch's t-test b/w interactions with Apo Open + closed
   
    diff_AD = mean_inter[0][n] - mean_inter[2][n]
    diff_BBR = mean_inter[0][n] - mean_inter[3][n]
    
    #Both AD amd BBR disrupt the interactions
    if (p2 < 0.1 or p3 < 0.05) and (diff_AD > 10 or diff_BBR > 10):
        output_lig_disrupt.write(str(pairs_all[n]) + ':\n'+ 'AD: ' + str(mean_inter[2][n]) + ' +/- ' + str(sem_inter[2][n]) + '\np to Apo open: ' + str(p2) + '\nBBR: ' + str(mean_inter[3][n]) + ' +/- ' + str(sem_inter[3][n])+ '\np to Apo open: ' + str(p2) + '\n')
        output_lig_disrupt.write('Apo Open: ' + str(mean_inter[0][n]) + '+/-' + str(sem_inter[0][n]) + '\n')
        output_lig_disrupt.write('Apo Closed: ' + str(mean_inter[1][n]) + '+/-' + str(sem_inter[1][n]) + '\n')

        box_plot(Apo_open, Apo_close, AD, BBR, pairs_all[n], 'Ligand_disrupt/', p2, p3)

    if (p2 < 0.1 or p3 < 0.05) and diff_AD < -10 and diff_BBR < -10:
        output_lig_inc.write(str(pairs_all[n]) + ':\n'+ 'AD: ' + str(mean_inter[2][n]) + ' +/- ' + str(sem_inter[2][n]) + '\np to Apo open: ' + str(p2) + '\nBBR: ' + str(mean_inter[3][n]) + ' +/- ' + str(sem_inter[3][n])+ '\np to Apo open: ' + str(p2) + '\n')
        output_lig_inc.write('Apo Open: ' + str(mean_inter[0][n]) + '+/-' + str(sem_inter[0][n]) + '\n')
        output_lig_inc.write('Apo Closed: ' + str(mean_inter[1][n]) + '+/-' + str(sem_inter[1][n]) + '\n')

        box_plot(Apo_open, Apo_close, AD, BBR, pairs_all[n], 'Ligand_inc/', p2, p3)

    diff_Apo = mean_inter[1][n] - mean_inter[0][n]
    if p4 < 0.05 and diff_Apo > 0:
        output_Apo_disrupt.write(str(pairs_all[n]) + ':\n'+ 'AD: ' + str(mean_inter[2][n]) + ' +/- ' + str(sem_inter[2][n]) + '\np to Apo open: ' + str(p2) + '\nBBR: ' + str(mean_inter[3][n]) + ' +/- ' + str(sem_inter[3][n])+ '\np to Apo open: ' + str(p2) + '\n')
        output_Apo_disrupt.write('Apo Open: ' + str(mean_inter[0][n]) + '+/-' + str(sem_inter[0][n]) + '\n')
        output_Apo_disrupt.write('Apo Closed: ' + str(mean_inter[1][n]) + '+/-' + str(sem_inter[1][n]) + '\n')

        box_plot(Apo_open, Apo_close, AD, BBR, pairs_all[n], 'Apo_disrupt/', p2, p3)

