#Import Necessary Packages
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from matplotlib_venn import venn2, venn3, venn2_circles

def load_data(path):
    raw = open('../../../' + path + '/Hbonds_bw_a3_sect_and_not_a3.txt').readlines()
    num = int(len(raw)/10 + 1)
    r = np.zeros(num)
    n=0
    for i in range(len(raw)):
        if i%10 == 0:
            r[n] = float(raw[i])
            n+=1
    return r

def load_names(path):
    names = open('../../../' + path + '/Hbonds_bw_a3_sect_and_not_a3_name.txt').readlines()
    return names

def process_data(data, WT_data):
    mean = np.mean(data)
    err = stats.sem(data)
    st, p = stats.ttest_ind(WT_data, data, equal_var = False) #Welch's t-test

    return mean, err, p

#File paths for all input files
file_path = ['../Apo_dis/analysis', '../Apo_1SUG/analysis/1sug', 'F196A/Apo/analysis', 'L192F/Apo/analysis', 'L195F/Apo/analysis', 'F280Y/Apo/analysis', 'E276F/Apo/analysis', 'V287T/Apo/analysis']

#Load WT
WT_data = load_data(file_path[0])

#Array of means, err, and p-values
mean = np.zeros(len(file_path))
err = np.zeros(len(file_path))
P = np.zeros(len(file_path))

#Set WT Values
mean[0] = np.mean(WT_data)
err[0] = stats.sem(WT_data)

#Load the rest of the data and caculate p-values
for i in range(1, len(file_path)):
    #Load Data
    dist = load_data(file_path[i])
    mean[i], err[i], P[i] = process_data(dist, WT_data)


#Name Labels
Label = ['WT', 'Apo Closed', 'F196A', 'L192F', 'L195F', 'F280Y', 'E276F', 'V287T']

bar_color = ['gray', 'gray', 'red', 'gray', 'pink', 'blue', 'pink', 'red']

num = np.linspace(1, len(Label)+1, num = len(Label))
fig = plt.figure(figsize=(18,10))
ax1 = fig.add_subplot(111)
ax1.set_title("Comparison of the Number of H-bonds with the bottom of the a3")
ax1.set_ylabel('Mean # H-bonds')
ax1.bar(num, mean, color = bar_color)
plt.xticks(num, Label, fontsize=14)
plt.errorbar(num, mean, yerr= err, fmt='o', color='black')
fig.savefig('a3_hbond_compare_.png') 

#Load names of H-bonds
names_WT = load_names(file_path[0])
names_ApoC = load_names(file_path[1])
names_F196A = load_names(file_path[2])
names_L192F = load_names(file_path[3])
names_L195F = load_names(file_path[4])
names_F280Y = load_names(file_path[5])
names_E276F = load_names(file_path[6])
names_V287T = load_names(file_path[7])

file_names = open('Hbond_names.txt', 'w')
file_names.write('WT\n')
for i in names_WT:
    file_names.write(i)
file_names.write('Apo Closed\n')
for i in names_ApoC:
    file_names.write(i)
file_names.write('F196A\n')
for i in names_F196A:
    file_names.write(i)
file_names.write('L192F\n')
for i in names_L192F:
    file_names.write(i)
file_names.write('L195F\n')
for i in names_L195F:
    file_names.write(i)
file_names.write('F280Y\n')
for i in names_F280Y:
    file_names.write(i)
file_names.write('E276F\n')
for i in names_E276F:
    file_names.write(i)
file_names.write('V287T\n')
for i in names_V287T:
    file_names.write(i)


