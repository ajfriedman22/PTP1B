from matplotlib import pyplot as plt
import numpy as np
from scipy import stats
import sys

#Import custom modules
sys.path.insert(1,'/ocean/projects/cts160011p/afriedma/code/PTP1B/util')
import plot

def load_file(dir_path, inters):
    a3, a4, a5, a6, a7 = [],[],[],[],[]
    for i in range(len(dir_path)):
        for j in range(len(inters)):
            data = open('../../../' + dir_path[i] + '/' + inters[j] + '_inter.txt').readlines()
            if j == 0:
                for k in data:
                    a3.append(float(k))
            if j == 1:
                for k in data:
                    a4.append(float(k))
            if j == 2:
                for k in data:
                    a5.append(float(k))
            if j == 3:
                for k in data:
                    a6.append(float(k))
            if j == 4:
                for k in data:
                    a7.append(float(k))
    return a3, a4, a5, a6, a7

def plot_helix_inter(helix, err, hel_num, p):
    num = [1, 2]
    label = ['AD', 'BBR']

    fig = plt.figure(figsize=(8,6))
    plt.title(r'Interactions with $\alpha$' + str(hel_num) + ' Helix', fontsize = 18)
    plt.ylabel('Mean Simultaneous Interactions', fontsize = 16)
    plt.bar(num, helix, color = ['blue', 'purple'], width=0.8)
    plt.errorbar(num, helix, yerr=err, fmt='o', color='black')
    plt.xticks(num, label, fontsize=16)
    plt.yticks(fontsize = 14)
    plot.error_bar(1, 2, helix[0], helix[1], p, 0.5, 'k')
    fig.savefig('a' + str(hel_num) + '_helix_inter_cmpr.png')
    plt.close(fig)

#List of all directory paths for each group
dir_path_AD_crys = ['mutate/WT/AD/analysis', '1sug_dis_AD/analysis/config_alt2']
dir_path_AD_alt = ['1sug_dis_AD/analysis/config_alt']
dir_path_AD_alt2 = ['AD_rebuild_a7/analysis', 'AD/analysis', '1sug_AD/analysis']
dir_path_BBR = ['mutate/WT/BBR/analysis', 'BBR_a7/analysis', 'BBR_dis/analysis/config11']

#List interactions of interest
inters = ['a3', 'a4', 'a5', 'a6', 'a7']

#Load all data
a3_AD_crys, a4_AD_crys, a5_AD_crys, a6_AD_crys, a7_AD_crys = load_file(dir_path_AD_crys, inters)
a3_AD_alt, a4_AD_alt, a5_AD_alt, a6_AD_alt, a7_AD_alt = load_file(dir_path_AD_alt, inters)
#a3_AD_alt2, a4_AD_alt2, a5_AD_alt2, a6_AD_alt2, a7_AD_alt2 = load_file(dir_path_AD_alt2, inters)
a3_BBR, a4_BBR, a5_BBR, a6_BBR, a7_BBR = load_file(dir_path_BBR, inters)


#Determine the mean and sem for each helix interactions
a3_AD = np.append(a3_AD_alt, a3_AD_crys)
a4_AD = np.append(a4_AD_alt, a4_AD_crys)
a5_AD = np.append(a5_AD_alt, a5_AD_crys)
a6_AD = np.append(a6_AD_alt, a6_AD_crys)
a7_AD = np.append(a7_AD_alt, a7_AD_crys)

#a3_AD_all = np.append(a3_AD, a3_AD_alt2)

helix_a3 = np.array([np.mean(a3_AD), np.mean(a3_BBR)])
helix_a4 = np.array([np.mean(a4_AD), np.mean(a4_BBR)])
helix_a5 = np.array([np.mean(a5_AD), np.mean(a5_BBR)])
helix_a6 = np.array([np.mean(a6_AD), np.mean(a6_BBR)])
helix_a7 = np.array([np.mean(a7_AD), np.mean(a7_BBR)])

err_a3 = np.array([stats.sem(a3_AD), stats.sem(a3_BBR)])
err_a4 = np.array([stats.sem(a4_AD), stats.sem(a4_BBR)])
err_a5 = np.array([stats.sem(a5_AD), stats.sem(a5_BBR)])
err_a6 = np.array([stats.sem(a6_AD), stats.sem(a6_BBR)])
err_a7 = np.array([stats.sem(a7_AD), stats.sem(a7_BBR)])

#Run Welch's test on BBR vs AD helix interactions for a3-a7
st_a3, p_a3 = stats.ttest_ind(a3_AD, a3_BBR, equal_var = False) #Welch's t-test b/w a3 interactions with AD binding and BBR binding
st_a4, p_a4 = stats.ttest_ind(a4_AD, a4_BBR, equal_var = False) #Welch's t-test b/w a4 interactions with AD binding and BBR binding
st_a5, p_a5 = stats.ttest_ind(a5_AD, a5_BBR, equal_var = False) #Welch's t-test b/w a4 interactions with AD binding and BBR binding
st_a6, p_a6 = stats.ttest_ind(a6_AD, a6_BBR, equal_var = False) #Welch's t-test b/w a6 interactions with AD binding and BBR binding
st_a7, p_a7 = stats.ttest_ind(a7_AD, a7_BBR, equal_var = False) #Welch's t-test b/w a7 interactions with AD binding and BBR binding

#Print output to file
output = open('Helix_inter.txt', 'w')
output.write('a3 Helix Interactions\n')
output.write('Full for AD: ' + str(helix_a3[0]) + ' +/- ' + str(err_a3[0]) + '\n')
output.write('Full for BBR: ' + str(helix_a3[1]) + ' +/- ' + str(err_a3[1]) + '\n')
output.write('p-value b/w AD and BBR: ' + str(p_a3) + '\n')

output.write('a4 Helix Interactions\n')
output.write('Full for AD: ' + str(helix_a4[0]) + ' +/- ' + str(err_a4[0]) + '\n')
output.write('Full for BBR: ' + str(helix_a4[1]) + ' +/- ' + str(err_a4[1]) + '\n')
output.write('p-value b/w AD and BBR: ' + str(p_a4) + '\n')

output.write('a5 Helix Interactions\n')
output.write('AD crys: ' + str(helix_a5[0]) + ' +/- ' + str(err_a5[0]) + '\n')
output.write('BBR: ' + str(helix_a5[1]) + ' +/- ' + str(err_a5[1]) + '\n')
output.write('p-value b/w AD and BBR: ' + str(p_a5) + '\n')

output.write('a6 Helix Interactions\n')
output.write('AD crys: ' + str(helix_a6[0]) + ' +/- ' + str(err_a6[0]) + '\n')
output.write('BBR: ' + str(helix_a6[1]) + ' +/- ' + str(err_a6[1]) + '\n')
output.write('p-value b/w AD and BBR: ' + str(p_a6) + '\n')

output.write('a7 Helix Interactions\n')
output.write('AD crys: ' + str(helix_a7[0]) + ' +/- ' + str(err_a7[0]) + '\n')
output.write('BBR: ' + str(helix_a7[1]) + ' +/- ' + str(err_a7[1]) + '\n')
output.write('p-value b/w AD and BBR: ' + str(p_a7) + '\n')

#Plot comparison of helix interactions for a3-a7
num = [1, 2]
label = ['AD', 'BBR']
#Plot Bar graph comparing averages for a3 helix
plot_helix_inter(helix_a3, err_a3, 3, p_a3)

#Plot Bar graph comparing averages for a4 helix
plot_helix_inter(helix_a4, err_a4, 4, p_a4)

#Plot Bar graph comparing averages for a5 helix
plot_helix_inter(helix_a5, err_a5, 5, p_a5)

#Plot Bar graph comparing averages for a6 helix
plot_helix_inter(helix_a6, err_a6, 6, p_a6)

#Plot Bar graph comparing averages for a7 helix
plot_helix_inter(helix_a7, err_a7, 7, p_a7)

#plot.plot_gen_box(a3_AD, a3_BBR, 'AD', 'BBR', p_a3, 12, 'Ligand', 'Mean Simultaneous Interactison', '', 'a3_helix_box_cmpr.png', 8, 8)
plot.plot_gen_box(a3_AD_all, a3_BBR, 'AD', 'BBR', p_a3, 12, 'Ligand', 'Mean Simultaneous Interactison', '', 'a3_helix_box_all_cmpr.png', 8, 8)

