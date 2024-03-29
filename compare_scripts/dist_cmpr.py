#Import Necessary Packages
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from itertools import product

#Import custom modules
sys.path.insert(1, '/ocean/projects/cts160011p/afriedma/PTP1B/util/')
import plot

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

    plot.error_bar(5, 15, dist[0], dist[2], p, 1, 'b'):
    plot.error_bar(5, 20, dist[0], dist[3], p1, 1, 'r'):
    plot.error_bar(5, 10, dist[0], dist[1], p2, 1, 'k'):

    fig.savefig(label + '_inter.png')
    plt.close(fig)

#Make open arrays for time and atomic distances
pairs = ['151_191_', '152_297_', '178_150_', '179_191_', '185_191_', '189_295_', '200_282_', '200_287_', '264_185_', '276_292_', '280_287_']
label = ['L11_a3', 'L11_a7', 'WPD_L11', 'WPD_a3_top', 'WPD_a3', 'a3_a7_top', 'a3_a6', 'a3_a7', 'a3_a6_top', 'a6_a7_top', 'a6_a7']
eq_time = [5, 5, 50, 5, 75, 5, 5, 5, 70, 5, 60, 30]
tot_time = [200, 200, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300]

#Load file path
file_path = open('dist_file_path.txt', 'r').readlines()

#Load distance for all pairs
for i in range(len(pairs)):
    d_AD, d_BBR = [],[]
    d_Apo_open, d_Apo_close = [],[]
    d_WT_AD, d_WT_BBR, d_a7, d_dis9, d_dis11, d_1sug, d_1sug2, d_1sug3 = [],[],[],[],[],[],[],[]
    d_BBR_a7, d_AD_dis11, d_AD_alt, d_AD_alt2 = [],[],[],[]

    #Input Data for a3 and a6 residue dist
    for j in open(file_path[0] + pairs[i] + 'dist.txt').readlines():
        d_WT_AD.append(float(j))
    for j in open(file_path[1] + pairs[i] + 'dist.txt').readlines():
        d_WT_BBR.append(float(j))
    for j in open(file_path[2] + pairs[i] + 'dist.txt').readlines():
        d_a7.append(float(j))
    for j in open(file_path[3] + pairs[i] + 'dist.txt').readlines():
        d_dis9.append(float(j))
    for j in open(file_path[4] + pairs[i] + 'dist.txt').readlines():
        d_dis11.append(float(j))
    for j in open(file_path[5] + pairs[i] + 'dist.txt').readlines():
        d_1sug.append(float(j))
    for j in open(file_path[6] + pairs[i] + 'dist.txt').readlines():
        d_1sug2.append(float(j))
    for j in open(file_path[7] + pairs[i] + 'dist.txt').readlines():
        d_1sug3.append(float(j))
    for j in open(file_path[8] + pairs[i] + 'dist.txt').readlines():
        d_BBR_a7.append(float(j))
    for j in open(file_path[9] + pairs[i] + 'dist.txt').readlines():
        d_AD_dis11.append(float(j))
    for j in open(file_path[10] + pairs[i] + 'dist.txt').readlines():
        d_AD_alt.append(float(j))
    for j in open(file_path[11] + pairs[i] + 'dist.txt').readlines():
        d_AD_alt2.append(float(j))
    
    #Seperate into categories
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

