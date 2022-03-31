#Determine the difference in helical interactions b/w AD and BBR
from matplotlib import pyplot as plt
import numpy as np
from scipy import stats

def load_file(file_dir):
    num_inter = []
    for i in open('../../../' + file_dir + '/all_iter_frac.txt', 'r').readlines():
        num_inter.append(float(i.split(' ')[1]))
    bond = 0
    num_frame = max(num_inter)
    all_inter = np.zeros(len(num_inter))
    for i in num_inter:
        all_inter[bond] = (i / num_frame) *100
        bond += 1

    return all_inter
#Bond names
group_l = [300]
group_p  = [186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199, 200, 264, 265, 266, 267, 268, 269, 270, 271, 272, 273, 274, 275, 276, 277, 278, 279, 280, 281, 282, 283, 284, 285, 286, 287, 288, 289, 290, 291, 292, 293, 294, 295, 296, 297, 298]

#Empty array for each set of percentages
num_bonds = len(group_p)

#Open input files or contacts over time
dir_list = ['1sug_dis_AD/analysis/config11', '1sug_dis_AD/analysis/config_alt', '1sug_dis_AD/analysis/config_alt2', 'AD_dis/analysis/config11', 'mutate/WT/AD/analysis', 'BBR_a7/analysis', 'BBR_1sug/analysis', 'BBR_dis/analysis/config9', 'BBR_1sug_dis/analysis/config7', 'mutate/WT/BBR/analysis'] 

#Empty array for each set of percentages
num_bonds = len(group_p)
all_inter = np.zeros([len(dir_list), num_bonds])

for i in range(len(dir_list)):
    all_inter[i,:] = load_file(dir_list[i])

#Output File for significant bonds
output_AD_freq = open('Sign_inter_AD_freq.txt', 'w')
output_AD_infreq = open('Sign_inter_AD_infreq.txt', 'w')
output_BBR_freq = open('Sign_inter_BBR_freq.txt', 'w')
output_BBR_infreq = open('Sign_inter_BBR_infreq.txt', 'w')
output_cmpr = open('Sign_inter_cmpr.txt', 'w')

#Determine the mean and sem for each helix interactions
mean_inter = np.zeros([4, num_bonds])
sem_inter = np.zeros([4, num_bonds])

for n in range(num_bonds):
    AD_open = [all_inter[0][n], all_inter[1][n], all_inter[2][n], all_inter[4][n]]
    AD_close = all_inter[3][n]
    BBR_open = [all_inter[5][n], all_inter[6][n], all_inter[9][n]]
    BBR_close = [all_inter[7][n], all_inter[8][n]]

    mean_inter[0][n] = np.mean(AD_open) #AD Open
    mean_inter[1][n] = AD_close #AD Closed
    mean_inter[2][n] = np.mean(BBR_open) #BBR Open
    mean_inter[3][n] = np.mean(BBR_close) #BBR Closed

    sem_inter[0][n] = stats.sem(AD_open) #AD Open
    sem_inter[1][n] = 0 #AD Closed
    sem_inter[2][n] = stats.sem(BBR_open) #BBR Open
    sem_inter[3][n] = stats.sem(BBR_close) #BBR Closed
    
    #Run Welch's test on BBR vs AD helix interactions for a3-a7
    st1, p1 = stats.ttest_ind(AD_open, BBR_open, equal_var = False) #Welch's t-test b/w interactions with AD open + BBR open
    #st2, p2 = stats.ttest_ind(AD_open, AD_close, equal_var = False) #Welch's t-test b/w interactions with AD open + AD closed
    st3, p3 = stats.ttest_ind(BBR_open, BBR_close, equal_var = False) #Welch's t-test b/w interactions with BBR open + closed
    
    diff = abs(mean_inter[0][n] - mean_inter[1][n])
    
    if mean_inter[0][n] >= 50:
        output_AD_freq.write(str(group_p[n]) + ' + AD:\n'+ 'Open: ' + str(mean_inter[0][n]) + ' +/- ' + str(sem_inter[0][n]) + '\nClosed: ' + str(mean_inter[1][n]) + '\n')
    if mean_inter[0][n] < 50:
        output_AD_infreq.write(str(group_p[n]) + ' + AD:\n'+ 'Open: ' + str(mean_inter[0][n]) + ' +/- ' + str(sem_inter[0][n]) + '\nClosed: ' + str(mean_inter[1][n]) + '\n')
    
    if mean_inter[2][n] >= 50:
        output_BBR_freq.write(str(group_p[n]) + ' + BBR:\n'+ 'Open: ' + str(mean_inter[2][n]) + ' +/- ' + str(sem_inter[2][n]) + '\nClosed: ' + str(mean_inter[3][n]) + ' +/- ' + str(sem_inter[3][n]) + '\np = ' + str(p3) + '\n')
    if mean_inter[2][n] < 50:
        output_BBR_infreq.write(str(group_p[n]) + ' + BBR:\n'+ 'Open: ' + str(mean_inter[2][n]) + ' +/- ' + str(sem_inter[2][n]) + '\nClosed: ' + str(mean_inter[3][n]) + ' +/- ' + str(sem_inter[3][n]) + '\np = ' + str(p3) + '\n')
   
    diff = abs(mean_inter[0][n] - mean_inter[3][n])

    if p1 < 0.05:
        output_cmpr.write(str(group_p[n]) + ' + AD/BBR:\n'+ 'AD: ' + str(mean_inter[0][n]) + ' +/- ' + str(sem_inter[0][n]) + '\nBBR: ' + str(mean_inter[3][n]) + ' +/- ' + str(sem_inter[3][n]) + '\np = ' + str(p1) + '\n')
        #Bar plot
        num = [1, 2]
        label = ['AD', 'BBR']
        per = np.array([mean_inter[0][n], mean_inter[2][n]])
        sem = np.array([sem_inter[0][n], sem_inter[2][n]])
        fig = plt.figure()
        plt.title('Interactions b/w Ligand and res' + str(group_p[n]))
        plt.ylabel('Mean Percent of Time Interaction is Formed')
        plt.bar(num, per, width=0.8)
        plt.errorbar(num, per, yerr=sem, fmt='o', color='black')
        plt.xticks(num, label, fontsize=8)
        if p1 < 0.05 and p1 > 0.01:
            x1, x2 = 1, 2 #Columns for AD crystal and BBR
            y, h, col = (1.1*per.max()) + 2, 2, 'k'
            plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
            plt.text((x1+x2)*0.5, y+h, "*" , ha='center', va='bottom', color=col)
        if p1 < 0.01 and p1 > 0.001:
            x1, x2 = 1, 2 #Columns for AD crystal and BBR
            y, h, col = (1.1*per.max()) + 2, 2, 'k'
            plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
            plt.text((x1+x2)*0.5, y+h, "**" , ha='center', va='bottom', color=col)
        if p1 < 0.001:
            x1, x2 = 1, 2 #Columns for AD crystal and BBR
            y, h, col = (1.1*per.max()) + 2, 2, 'k'
            plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
            plt.text((x1+x2)*0.5, y+h, "***" , ha='center', va='bottom', color=col)
        fig.savefig('AD_BBR_' + str(group_p[n]) +'_cmpr.png')
        plt.close(fig)


