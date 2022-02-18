#Import Necessary Packages
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

def load_rmsf(path):
    raw = open('../../' + path + '/analysis/rmsf_ref.txt').readlines()
    r = np.zeros(len(raw))
    for i in range(len(raw)):
        r[i] = float(raw[i])
    return r

def plot_two(RMSF_all, ind, Label):
    #Seperate Data for Plotting
    rmsf_0 = RMSF_all[ind[0],:]
    rmsf_1 = RMSF_all[ind[1],:]
    label_0 = Label[ind[0]]
    label_1 = Label[ind[1]]
    
    #Establish atom numbers
    res = np.linspace(1, 299, num=len(rmsf_0))
    
    #Plot
    fig = plt.figure(figsize=(18,10))
    ax1 = fig.add_subplot(111)
    ax1.set_title("Comparison of RMSF between " + str(label_0) + " and " + str(label_1))
    ax1.set_ylabel('RMSF(nm)')
    ax1.set_xlabel('Residue Number')
    ax1.plot(res, rmsf_0, label=label_0)
    ax1.plot(res, rmsf_1, label=label_1)
    ax1.legend(loc='upper right')
    fig.savefig('RMSF_compare_' + label_0 + '_' + label_1 + '.png') 

def plot_three(RMSF_all, ind, Label):
    #Seperate Data for Plotting
    rmsf_0 = RMSF_all[ind[0],:]
    rmsf_1 = RMSF_all[ind[1],:]
    rmsf_2 = RMSF_all[ind[2],:]
    label_0 = Label[ind[0]]
    label_1 = Label[ind[1]]
    label_2 = Label[ind[2]]

    #Establish atom numbers
    res = np.linspace(1, 299, num=len(rmsf_0))
    
    #Plot
    fig = plt.figure(figsize=(18,10))
    ax1 = fig.add_subplot(111)
    ax1.set_title("Comparison of RMSF between " + str(label_0) + " and " + str(label_1) + " and " + str(label_2))
    ax1.set_ylabel('RMSF(nm)')
    ax1.set_xlabel('Residue Number')
    ax1.plot(res, rmsf_0, label=label_0)
    ax1.plot(res, rmsf_1, label=label_1)
    ax1.plot(res, rmsf_2, label=label_2)
    ax1.legend(loc='upper right')
    fig.savefig('RMSF_compare_' + label_0 + '_' + label_1 + '_' + label_2 + '.png') 

def plot_three_zoom(RMSF_all, ind, Label, res_range):
    #Seperate Data for Plotting
    rmsf_0all = RMSF_all[ind[0],:]
    rmsf_1all = RMSF_all[ind[1],:]
    rmsf_2all = RMSF_all[ind[2],:]
    label_0 = Label[ind[0]]
    label_1 = Label[ind[1]]
    label_2 = Label[ind[2]]

    #Establish atom numbers
    res_all = np.linspace(1, 299, num=len(rmsf_0all))
    
    #Limit to residue range
    res, rmsf_0, rmsf_1, rmsf_2 = [],[],[],[]
    for i in range(len(res_all)):
        if res_all[i] >= res_range[0] and res_all[i] <= res_range[1]:
            rmsf_0.append(rmsf_0all[i])
            rmsf_1.append(rmsf_1all[i])
            rmsf_2.append(rmsf_2all[i])
            res.append(res_all[i])
    #Plot
    fig = plt.figure(figsize=(18,10))
    ax1 = fig.add_subplot(111)
    ax1.set_title("Comparison of RMSF between " + str(label_0) + " and " + str(label_1) + " and " + str(label_2))
    ax1.set_ylabel('RMSF(nm)')
    ax1.set_xlabel('Residue Number')
    ax1.plot(res, rmsf_0, label=label_0)
    ax1.plot(res, rmsf_1, label=label_1)
    ax1.plot(res, rmsf_2, label=label_2)
    ax1.legend(loc='upper right')
    fig.savefig('RMSF_compare_' + label_0 + '_' + label_1 + '_' + label_2 + '_' + str(res_range[0]) + '_' + str(res_range[1]) + '.png') 

#File paths for all input files
file_path = ['../Apo_dis', 'F196A/Apo', 'L192F/Apo', 'L195F/Apo', 'F280Y/Apo', 'E276F/Apo', 'V287T/Apo']

#Load WT
rmsf = load_rmsf(file_path[0])

#Set empty array sizes
RMSF = np.zeros((len(file_path), len(rmsf)))
RMSF[0][:] = rmsf

for i in range(1, len(file_path)):
    #Load Data
    RMSF[i][:] = load_rmsf(file_path[i])
print(RMSF[0,:])

#Name Labels
Label = ['WT', 'F196A', 'L192F', 'L195F', 'F280Y', 'E276F', 'V287T']

#Plot all compared to WT
for i in range(1, len(Label)):
    plot_two(RMSF, [0,i], Label)

plot_three(RMSF, [0, 6, 4], Label)
plot_three(RMSF, [0, 2, 6], Label)
plot_three(RMSF, [0, 1, 6], Label)

#Plot only the substrate binding loop
plot_three_zoom(RMSF, [0, 6, 4], Label, [113, 120])
plot_three_zoom(RMSF, [0, 2, 6], Label, [113, 120])
plot_three_zoom(RMSF, [0, 1, 6], Label, [113, 120])

