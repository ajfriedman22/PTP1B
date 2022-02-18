#Import Necessary Packages
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import seaborn as sns

def load_rmsf(path, sect):
    raw = open('../../' + path + '/analysis/rmsd_' + sect + '_ref.txt').readlines()
    r = np.zeros(len(raw))
    for i in range(len(raw)):
        r[i] = float(raw[i])

    return r

def plot_compare(RMSD_mean, RMSD_err, Label, sect, n):
    rmsd_n = RMSD_mean[:,n]
    rmsd_err_n = RMSD_err[:,n]

    section = sect[n]
    bar_color = ['gray', 'red', 'gray', 'pink', 'blue', 'pink', 'red']

    num = np.linspace(1, len(Label)+1, num = len(Label))
    fig = plt.figure(figsize=(18,10))
    ax1 = fig.add_subplot(111)
    ax1.set_title("Comparison of RMSD for " + section)
    ax1.set_ylabel('RMSD(nm)')
    ax1.bar(num, rmsd_n, color = bar_color)
    plt.xticks(num, Label, fontsize=14)
    plt.errorbar(num, rmsd_n, yerr= rmsd_err_n, fmt='o', color='black')
    fig.savefig('RMSD_compare_' + section + '.png') 


#File paths for all input files
file_path = ['../Apo_dis', 'L192F/Apo', 'E276F/Apo', 'F280Y/Apo', 'L195F/Apo', 'F196A/Apo', 'V287T/Apo'] #Indices to rank in order of closest activity to WT to Furthest
sections = ['WPD', 'P', 'a3', 'a4', 'a5', 'a6', 'a7', 'Q']

#open all files
RMSD_mean = np.zeros((len(file_path), len(sections)))
RMSD_err = np.zeros((len(file_path), len(sections)))


for i in range(len(file_path)):
    for j in range(len(sections)):
        #Load Data
        rmsd = load_rmsf(file_path[i], sections[j])
        #Mean and SEM for each trajectory
        RMSD_mean[i][j] = np.mean(rmsd)
        RMSD_err[i][j] = stats.sem(rmsd)
print(np.shape(RMSD_mean))

#Name Labels
Label = ['WT', 'L192F', 'E276F', 'F280Y', 'L195F', 'F196A', 'V287T']
Labels_mut = ['L192F', 'E276F', 'F280Y', 'L195F', 'F196A', 'V287T']

#Plot all compared to WT
for i in range(len(sections)):
    plot_compare(RMSD_mean, RMSD_err, Label, sections, i)

#Determine % difference from WT
RMSD_diff = np.zeros((len(sections), len(Labels_mut)))
for i in range(1, len(Label)):
    n = i-1
    for j in range(len(sections)):
        WT = RMSD_mean[0][j]
        Mut = RMSD_mean[i][j]
        RMSD_diff[j][n] = ((Mut-WT)/((Mut+WT)/2)) * 100

#Plot table comparing residue interactions to WT
ax = plt.figure(figsize=(12, 6), frameon=False) # no visible frame
ax = sns.heatmap(RMSD_diff, annot=False, cmap = 'jet', cbar = True, cbar_kws={'label': 'Percentage Difference from WT'}, vmin = 0, vmax = 150, xticklabels = Labels_mut, yticklabels = sections)
#ax.add_artist(lines.Line2D([0, 20], [7, 7], color = 'black', linestyle= '--', linewidth = 4))
plt.title('Section RMSD Compared to WT')
plt.savefig('mutate_RMSD_Apo.png')
plt.close()

