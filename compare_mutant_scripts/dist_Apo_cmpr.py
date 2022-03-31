#Import Necessary Packages
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from itertools import product
def load_rmsf(path, pair, WT_data):
    raw = open('../../../' + path + '/analysis/' + pair + '_dist.txt').readlines()
    num = int(len(raw)/250 + 1)
    r = np.zeros(num)
    n = 0
    for i in range(len(raw)):
        if i%250 == 0:
            r[n] = float(raw[i])
            n += 1
    mean = np.mean(r)
    err = stats.sem(r)
    st, p = stats.ttest_ind(WT_data, r, equal_var = False) #Welch's t-test
    return mean, err, p

def plot_compare(mean, err, Label, sect, n):
    mean_n = mean[:,n]
    err_n = err[:,n]

    section = sect[n]
    bar_color = ['gray', 'gray', 'pink', 'blue', 'pink', 'red', 'red']

    num = np.linspace(1, len(Label)+1, num = len(Label))
    fig = plt.figure(figsize=(18,10))
    ax1 = fig.add_subplot(111)
    ax1.set_title("Comparison of Residue distance for " + section)
    ax1.set_ylabel('Distance(nm)')
    ax1.bar(num, mean_n, color = bar_color)
    plt.xticks(num, Label, fontsize=14)
    plt.errorbar(num, mean_n, yerr= err_n, fmt='o', color='black')
    fig.savefig('dist_compare_' + section + '.png') 

#File paths for all input files
file_path = ['../Apo_dis', 'L192F/Apo', 'E276F/Apo', 'F280Y/Apo', 'L195F/Apo', 'F196A/Apo', 'V287T/Apo'] #Indices to rank in order of closest activity to WT to Furthest
pairs = ['81_199', '117_182', '117_217', '192_225', '182_220', '187_269', '178_190', '152_177', '151_191', '152_297', '178_150', '179_191', '182_220', '185_191', '187_269', '189_295', '192_225', '200_282', '200_287', '264_185', '276_292', '280_287']
label = ['WT', 'L192F', 'E276F', 'F280Y', 'L195F', 'F196A', 'V287T']
mut = ['L192F', 'E276F', 'F280Y', 'L195F', 'F196A', 'V287T']

#Load WT
raw = open('../../../' + file_path[0] + '/analysis/' + pairs[0] + '_dist.txt').readlines()
num = int(len(raw)/250 + 1)
WT_data = np.zeros(num)
n = 0
for i in range(len(raw)):
    if i%250 == 0:
        WT_data[n] = float(raw[i])
        n += 1

#Empty Vector for means
Dist_mean = np.zeros((len(file_path), len(pairs)))
Dist_err = np.zeros((len(file_path), len(pairs)))

Dist_mean[0,0] = np.mean(WT_data)
Dist_err[0,0] = stats.sem(WT_data)

p = np.zeros((len(mut), len(pairs)))

#Load the rest of the data
for i in range(1, len(file_path)):
    n = i-1
    for j in range(len(pairs)):
        if j >0:
            raw = open('../../../' + file_path[0] + '/analysis/' + pairs[j] + '_dist.txt').readlines()
            WT_data = np.zeros(len(raw))
            for k in range(len(raw)):
                WT_data[k] = float(raw[k])
            Dist_mean[0,j] = np.mean(WT_data)
            Dist_err[0,j] = stats.sem(WT_data)


        Dist_mean[i][j], Dist_err[i][j], p[n][j] = load_rmsf(file_path[i], pairs[j], WT_data)

#Plot compare
for i in range(len(pairs)):
    plot_compare(Dist_mean, Dist_err, label, pairs, i)

