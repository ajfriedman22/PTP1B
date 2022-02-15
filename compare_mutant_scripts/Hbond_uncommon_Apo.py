
#import required packages
import mdtraj as md
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt

def import_data(file_path):
    file_bonds = open(file_path + '/analysis/Hbonds_per_Apo_uncommon.txt','r').readlines() #Read file with the percentage of each bond
    per = np.zeros(len(file_bonds)) #Declare empty array for the converted percentages
    for i in range(len(file_bonds)): #convert percentages to floats
        per[i] = float(file_bonds[i])
    return per

def plot_bond(save_path, hbond_name, mean, err, mean_all):
    num = [0, 1, 2, 3, 4, 5, 6]
    Label = ['WT', 'F196A', 'L192F', 'L195F', 'E276F', 'F280Y', 'V287T']
    fig = plt.figure(figsize=(14,8))
    ax1 = fig.add_subplot(111)
    ax1.set_title('Comparison of Presence of bond b/w ' + str(hbond_name)) 
    ax1.set_ylabel('Percent time h-bond was formed')
    ax1.bar(num, mean_all, width=0.9)
    plt.xticks(num, Label, fontsize=8)
    fig.savefig(save_path + '/' + hbond_name + '_all.png')
    plt.close(fig)         

    num = [1, 2, 3, 4]
    Label = ['~WT', '>WT', '>>WT', '<<WT']
    fig = plt.figure(figsize=(14,8))
    ax1 = fig.add_subplot(111)
    ax1.set_title('Comparison of Presence of bond b/w ' + str(hbond_name)) 
    ax1.set_ylabel('Percent time h-bond was formed')
    ax1.bar(num, mean, width=0.9)
    plt.errorbar(num, mean, yerr= err, fmt='o', color='black')
    plt.xticks(num, Label, fontsize=8)
    fig.savefig(save_path + '/' + hbond_name + '.png')
    plt.close(fig)         


#Load h-bond names file
hbond_names = open('Hbonds_uncommon_Apo.txt', 'r').readlines()

#Load all data for percentages of h-bonds
file_paths = ['../../../Apo_dis', '../../F196A/Apo', '../../L192F/Apo', '../../L195F/Apo', '../../E276F/Apo', '../../F280Y/Apo', '../../V287T/Apo']
per_hbonds = np.zeros((len(file_paths), len(hbond_names)))
for i in range(len(file_paths)):
    per_hbonds[i][:] = import_data(file_paths[i])

#Average for all bonds into 4 groups close to WT activity, slightly above WT activity, far above WT activity, and below WT activity
group_per = np.zeros((4, len(hbond_names))) #Mean of groups
group_per_err = np.zeros((4, len(hbond_names))) #SEM of groups

for i in range(len(hbond_names)):
    group_per[0][i] = np.mean([per_hbonds[0][i], per_hbonds[2][i]])#similar to WT
    group_per[1][i] = np.mean([per_hbonds[3][i], per_hbonds[4][i]])#slightly above WT
    group_per[2][i] = np.mean([per_hbonds[1][i], per_hbonds[6][i]])#significanty above WT
    group_per[3][i] = per_hbonds[5][i] #below WT

    group_per_err[0][i] = stats.sem([per_hbonds[0][i], per_hbonds[2][i]])
    group_per_err[1][i] = stats.sem([per_hbonds[3][i], per_hbonds[4][i]])
    group_per_err[2][i] = stats.sem([per_hbonds[1][i], per_hbonds[6][i]])
    group_per_err[3][i] = 0 #Only one value

#Determine bonds correlated to Apo activity level
for i in range(len(hbond_names)):
    mean = group_per[:,i]
    err = group_per_err[:,i]
    mean_all = per_hbonds[:,i]
    #Present more with increased activity
    diff = mean[2] - mean[0]
    diff1 = mean[1] - mean[0]
#    if diff > 30:
#        plot_bond('more_inc', hbond_names, i, group_per, group_per_err)
    #Present less with increased activity
#    if diff < 30:
#        plot_bond('less_inc', hbond_names, i, group_per, group_per_err)

    #Present more with decreased activity
    diff2 = mean[3] - mean[0]
#    if diff2 > 30:
#        plot_bond('more_dec', hbond_names, i, group_per, group_per_err)
    #Present less with increased activity
#    if diff2 < 30:
#        plot_bond('less_dec', hbond_names, i, group_per, group_per_err)
    
    #If bond is present more with increased activity and less with decreased activity
    if diff > 15 and diff2 < -15 and (diff1 > 0 or abs(diff1) < 5):
        plot_bond('less_dec_more_inc', hbond_names[i], mean, err, mean_all)

    #vice versa
    elif diff < -15 and diff2 > 15 and (diff1 < 0 or abs(diff1) < 5):
        plot_bond('more_dec_less_inc', hbond_names[i], mean, err, mean_all)
    
    #If bond is just more present with increased activity
    elif diff > 15 and (diff1 > 0 or abs(diff1) < 5):
        plot_bond('more_inc', hbond_names[i], mean, err, mean_all)
    
    #If bond is just less present with decreased activity
    elif diff2 < -15:
        plot_bond('less_dec', hbond_names[i], mean, err, mean_all)
    
    #If bond is just less present with increased activity
    elif diff < -15 and (diff1 < 0 or abs(diff1) < 5):
        plot_bond('less_inc', hbond_names[i], mean, err, mean_all)

    #If bond is just more present with decreased activity
    elif diff2 > 15:
        plot_bond('more_dec', hbond_names[i], mean, err, mean_all)
   
    
    #Specific requests
    if hbond_names[i] == 'THR177-N -- TYR152-O\n':
        plot_bond('.', hbond_names[i], mean, err, mean_all)



