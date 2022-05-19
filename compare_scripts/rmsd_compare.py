#Import Necessary Packages
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import seaborn as sns
import ruptures as rpt
from statistics import stdev
import pandas as pd

def load_rmsf(path, sect, ref):
    raw_string = open('../../' + path + '/rmsd_' + sect + '_ref_' + ref + '.txt').readlines()
    #Convert data to fload
    raw = np.zeros(len(raw_string))
    for i in range(len(raw_string)):
        raw[i] = float(raw_string[i])*10
    return raw

def plot_compare(RMSD_mean, RMSD_err, Label, sect, n, ref):
    rmsd_n = RMSD_mean[:,n]
    rmsd_err_n = RMSD_err[:,n]

    section = sect[n]
#    bar_color = ['gray', 'gray', 'pink', 'blue', 'pink', 'red', 'red']

    num = np.linspace(1, len(Label)+1, num = len(Label))
    fig = plt.figure(figsize=(18,10))
    ax1 = fig.add_subplot(111)
    ax1.set_title("Comparison of RMSD for " + section + ' to ' + ref)
    ax1.set_ylabel(r'RMSD($\AA$)')
    ax1.bar(num, rmsd_n)
    plt.xticks(num, Label, fontsize=14)
    plt.errorbar(num, rmsd_n, yerr= rmsd_err_n, fmt='o', color='black')
    fig.savefig('RMSD_compare_' + section + '_' + ref + '.png') 
    plt.close(fig)

def plot_kernel_cmpr_lig(apo_df, AD_df, BBR_df, mut, sect, n):
    df = pd.concat([apo_df, AD_df, BBR_df])

    sns.kdeplot(data = df, fill=True, alpha=0.5, common_norm = True, common_grid = True)
    plt.axvline(x = RMSD_mean_close[0,n], color = 'r')
    plt.axvline(x = RMSD_mean_close[1,n], color = 'b')
    plt.xlabel(r'RMSD($\AA$)')
    plt.ylabel(r'Normalized Density')
    plt.title(sect + ' RMSD Compared to WT Closed')
    plt.savefig('mutate_RMSD_' + sect + '_' + mut + '.png')
    plt.close()

def sort_data(file_path, sections, ref, RMSD_mean, RMSD_err):
    RMSD = []
    for i in range(len(file_path)):
        for j in range(len(sections)):
            #Load Data for reference open
            rmsd = load_rmsf(file_path[i], sections[j], ref[0])

            #Mean and SEM for each trajectory
            RMSD_mean[i][j] = np.mean(rmsd)
            RMSD_err[i][j] = stats.sem(rmsd)
        
            #if section is for the a3_top
            if sections[j] == 'a3_top':
                for k in rmsd:
                    RMSD.append(k)
    return RMSD, RMSD_mean, RMSD_err


#File paths for all input files
file_path_Apo_close = ['/Apo_1SUG/analysis/1sug', 'Apo_1SUG/analysis/1sug3']
file_path_Apo_open = ['Apo_dis/analysis', 'rebuild_a7_high/config11/analysis']
file_path_AD = ['1sug_dis_AD/analysis/config11', '1sug_dis_AD/analysis/config_alt', 'mutate/WT/AD/analysis']
file_path_BBR = ['BBR_a7/analysis', 'BBR_dis/analysis/config7', 'BBR_dis/analysis/config9', 'mutate/WT/BBR/analysis']

sections = ['WPD', 'WPD_a3', 'SBL', 'beg', 'P', 'CYS', 'a3', 'a3_top', 'a4', 'a4_P', 'a5', 'a6', 'a6_bot', 'a7', 'Q']
ref = ['closed']

#open all files
RMSD_mean_Apo_open = np.zeros((len(file_path_Apo_open), len(sections))) #Mean for reference open
RMSD_err_Apo_open = np.zeros((len(file_path_Apo_open), len(sections))) #SEM for reference open
RMSD_mean_Apo_close = np.zeros((len(file_path_Apo_close), len(sections))) #Mean for reference closed
RMSD_err_Apo_close = np.zeros((len(file_path_Apo_close), len(sections))) #SEM for reference closed
RMSD_mean_AD = np.zeros((len(file_path_AD), len(sections))) #Mean for reference closed
RMSD_err_AD = np.zeros((len(file_path_AD), len(sections))) #SEM for reference closed
RMSD_mean_BBR = np.zeros((len(file_path_BBR), len(sections))) #Mean for reference closed
RMSD_err_BBR = np.zeros((len(file_path_BBR), len(sections))) #SEM for reference closed

#Save all rmsd values for a3_top, a4_top, and a6 helix
RMSD_Apo_open, RMSD_mean_Apo_open, RMSD_err_Apo_open = sort_data(file_path_Apo_open, sections, ref, RMSD_mean_Apo_open, RMSD_err_Apo_open)
RMSD_Apo_close, RMSD_mean_Apo_close, RMSD_err_Apo_close = sort_data(file_path_Apo_close, sections, ref, RMSD_mean_Apo_close, RMSD_err_Apo_close)
RMSD_AD, RMSD_mean_AD, RMSD_err_AD = sort_data(file_path_AD, sections, ref, RMSD_mean_AD, RMSD_err_AD)
RMSD_BBR, RMSD_mean_BBR, RMSD_err_BBR = sort_data(file_path_BBR, sections, ref, RMSD_mean_BBR, RMSD_err_BBR)

#Compare a3_top
a3_Apo_open_df = pd.DataFrame({'Apo Open': RMSD_Apo_open})
a3_Apo_close_df = pd.DataFrame({'Apo Closed': RMSD_Apo_close})
a3_AD_df = pd.DataFrame({'AD': RMSD_AD})
a3_BBR_df = pd.DataFrame({'BBR': RMSD_BBR})

df = pd.concat([a3_Apo_open_df, a3_Apo_close_df, a3_AD_df, a3_BBR_df])

ax = plt.figure(figsize=(12, 6), frameon=False) # no visible frame
sns.kdeplot(data = df, fill=True, alpha=0.5, common_norm = True, common_grid = True)
plt.xlabel(r'RMSD($\AA$)')
plt.xlim(0, 2)
plt.ylabel(r'Normalized Density')
plt.title(r'Top of $\alpha$-3 RMSD Compared to WT Closed')
plt.savefig('mutate_RMSD_a3_top_all.png')
plt.close()

#Determine p-values for each of the sections of focus
file_p = open('p_values_mut.txt', 'w')
p = np.zeros((1, 3))
st, p[0,0] = stats.ttest_ind(RMSD_Apo_open, RMSD_AD, equal_var = False) #Welch's t-test
st, p[0,1] = stats.ttest_ind(RMSD_Apo_open, RMSD_BBR, equal_var = False) #Welch's t-test
st, p[0,2] = stats.ttest_ind(RMSD_Apo_open, RMSD_Apo_close, equal_var = False) #Welch's t-test

sections_mini = ['a3']
Labels = ['AD', 'BBR', 'Apo Closed']
file_p.write('P values of RMSD with Apo closed reference structure Relative to Apo Open')
for i in range(len(sections_mini)):
    file_p.write(str(sections_mini[i]) + '\n') 
    for j in range(len(Labels)):
        file_p.write(Labels[j] + ': ' + str(p[i,j]) + '\n')

