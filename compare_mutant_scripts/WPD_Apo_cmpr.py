#Import Necessary Packages
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from itertools import product
def load_rmsf(path, WT_data, name):
    raw = []
    with open('../../' + path + '/analysis/WPD_' + name + '.xvg') as f:
        for _ in range(17):
            next(f)
        for line in f:
            cols = line.split()
            if len(cols) == 2:
                raw.append(float(cols[1]))
    num = int(len(raw)/250) + 1
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

def plot_compare(mean, err, Label):
    bar_color = ['gray', 'gray', 'pink', 'blue', 'pink', 'red', 'red']

    num = np.linspace(1, len(Label)+1, num = len(Label))
    fig = plt.figure(figsize=(18,10))
    ax1 = fig.add_subplot(111)
    ax1.set_title("Comparison of WPD Loop Conformation")
    ax1.set_ylabel('Distance(nm)')
    ax1.bar(num, mean, color = bar_color)
    plt.xticks(num, Label, fontsize=14)
    plt.errorbar(num, mean, yerr= err, fmt='o', color='black')
    fig.savefig('WPD_compare_Apo.png') 

#File paths for all input files
file_path = ['../Apo_dis', 'L192F/Apo', 'E276F/Apo', 'F280Y/Apo', 'L195F/Apo', 'F196A/Apo', 'V287T/Apo'] #Indices to rank in order of closest activity to WT to Furthest
mut = ['Apo_alt', 'L192F_Apo', 'E276F_Apo', 'F280Y_Apo', 'L195F_Apo', 'F196A_Apo', 'V287T_Apo']

#Load WT
raw = []
with open('../../' + file_path[0] + '/analysis/WPD_' + mut[0] + '.xvg') as f:
    for _ in range(17):
        next(f)
    for line in f:
        cols = line.split()
        if len(cols) == 2:
            raw.append(float(cols[1]))
num = int(len(raw)/250) + 1
WT_data = np.zeros(num)
n = 0
for i in range(len(raw)):
    if i%250 == 0:
        WT_data[n] = float(raw[i])
        n += 1

#Empty Vector for means
Dist_mean = np.zeros(len(file_path))
Dist_err = np.zeros(len(file_path))

Dist_mean[0] = np.mean(WT_data)
Dist_err[0] = stats.sem(WT_data)

p = np.zeros(len(mut))

#Load the rest of the data
for i in range(1, len(file_path)):
    n = i-1
    Dist_mean[i], Dist_err[i], p[n] = load_rmsf(file_path[i], WT_data, mut[i])

label = ['WT', 'L192F', 'E276F', 'F280Y', 'L195F', 'F196A', 'V287T'] 
#Plot compare
plot_compare(Dist_mean, Dist_err, label)

