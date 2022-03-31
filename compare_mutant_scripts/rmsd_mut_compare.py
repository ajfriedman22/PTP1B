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



#File paths for all input files
file_path = ['../Apo_dis/analysis', 'L192F/Apo/analysis', 'E276F/Apo/analysis', 'F280Y/Apo/analysis', 'L195F/Apo/analysis', 'F196A/Apo/analysis', 'V287T/Apo/analysis'] #Indices to rank in order of closest activity to WT to Furthest

file_path_close = ['../Apo_1SUG/analysis/1sug', '../Apo_dis/analysis', 'L192F/Apo/analysis', 'E276F/Apo/analysis', 'F280Y/Apo/analysis', 'L195F/Apo/analysis', 'F196A/Apo/analysis', 'V287T/Apo/analysis']
file_path_close_AD = ['L192F/AD/analysis', 'E276F/AD/analysis', 'F280Y/AD/analysis', 'L195F/AD/analysis', 'F196A/AD/analysis', 'V287T/AD/analysis']
file_path_close_BBR = ['L192F/BBR/analysis', 'E276F/BBR/analysis', 'F280Y/BBR/analysis', 'L195F/BBR/analysis', 'F196A/BBR/analysis', 'V287T/BBR/analysis']

sections = ['WPD', 'WPD_a3', 'SBL', 'beg', 'P', 'CYS', 'a3', 'a3_top', 'a4', 'a4_P', 'a5', 'a6', 'a6_bot', 'a7', 'Q']
ref = ['open', 'closed', 'self', 'F196A', 'V287T']

#open all files
RMSD_mean = np.zeros((len(file_path), len(sections))) #Mean for reference open
RMSD_err = np.zeros((len(file_path), len(sections))) #SEM for reference open
RMSD_mean_close = np.zeros((len(file_path_close), len(sections))) #Mean for reference closed
RMSD_err_close = np.zeros((len(file_path_close), len(sections))) #SEM for reference closed
RMSD_mean_close_AD = np.zeros((len(file_path_close_AD), len(sections))) #Mean for reference closed
RMSD_err_close_AD = np.zeros((len(file_path_close_AD), len(sections))) #SEM for reference closed
RMSD_mean_close_BBR = np.zeros((len(file_path_close_BBR), len(sections))) #Mean for reference closed
RMSD_err_close_BBR = np.zeros((len(file_path_close_BBR), len(sections))) #SEM for reference closed

#Save all rmsd values for a3_top, a4_top, and a6 helix

for i in range(len(file_path_close)):
    for j in range(len(sections)):
        if i < len(file_path):
            #Load Data for reference open
            rmsd = load_rmsf(file_path[i], sections[j], ref[0])
            #Mean and SEM for each trajectory
            RMSD_mean[i][j] = np.mean(rmsd)
            RMSD_err[i][j] = stats.sem(rmsd)

        #Load Data for reference open
        rmsd_Apo = load_rmsf(file_path_close[i], sections[j], ref[1])
        rmsd_Apo_rapo = load_rmsf(file_path_close[i], sections[j], ref[3])

        #Mean and SEM for each trajectory
        RMSD_mean_close[i][j] = np.mean(rmsd)
        RMSD_err_close[i][j] = stats.sem(rmsd)
        if i < len(file_path_close_AD):#Mean and RMSD for AD and BBR bound mutants
            #Load Data for reference open
            rmsd_AD_rapo = load_rmsf(file_path_close_AD[i], sections[j], ref[3])
            rmsd_AD = load_rmsf(file_path_close_AD[i], sections[j], ref[1])
            RMSD_mean_close_AD[i][j] = np.mean(rmsd_AD)
            RMSD_err_close_AD[i][j] = stats.sem(rmsd_AD)

            rmsd_BBR_rapo = load_rmsf(file_path_close_BBR[i], sections[j], ref[3])
            rmsd_BBR = load_rmsf(file_path_close_BBR[i], sections[j], ref[1])
            RMSD_mean_close_BBR[i][j] = np.mean(rmsd_BBR)
            RMSD_err_close_BBR[i][j] = stats.sem(rmsd_BBR)

        if sections[j] == 'a3_top':
            if i == 0:
                rmsd_a3_top_1sug = rmsd_Apo
                rmsd_a3_top_L192F_AD = rmsd_AD
                rmsd_a3_top_L192F_BBR = rmsd_BBR
            if i == 1:
                rmsd_a3_top_apo = rmsd_Apo
                rmsd_a3_top_E276F_AD = rmsd_AD
                rmsd_a3_top_E276F_BBR = rmsd_BBR
            if i == 2:
                rmsd_a3_top_L192F = rmsd_Apo
                rmsd_a3_top_F280Y_AD = rmsd_AD
                rmsd_a3_top_F280Y_BBR = rmsd_BBR
            if i == 3:
                rmsd_a3_top_E276F = rmsd_Apo
                rmsd_a3_top_L195F_AD = rmsd_AD
                rmsd_a3_top_L195F_BBR = rmsd_BBR
            if i == 4:
                rmsd_a3_top_F280Y = rmsd_Apo
                rmsd_a3_top_F196A_AD = rmsd_AD
                rmsd_a3_top_F196A_BBR = rmsd_BBR
            if i == 5:
                rmsd_a3_top_L195F = rmsd_Apo
                rmsd_a3_top_V287T_AD = rmsd_AD
                rmsd_a3_top_V287T_BBR = rmsd_BBR
            if i == 6:
                rmsd_a3_top_F196A = rmsd_Apo
            if i == 7:
                rmsd_a3_top_V287T = rmsd_Apo

        if sections[j] == 'a4':
            if i == 0:
                rmsd_a4_1sug = rmsd_Apo
                rmsd_a4_1sug_rapo = rmsd_Apo_rapo
                rmsd_a4_L192F_AD = rmsd_AD
                rmsd_a4_L192F_BBR = rmsd_BBR
            if i == 1:
                rmsd_a4_apo = rmsd_Apo
                rmsd_a4_apo_rapo = rmsd_Apo_rapo
                rmsd_a4_E276F_AD = rmsd_AD
                rmsd_a4_E276F_BBR = rmsd_BBR
            if i == 2:
                rmsd_a4_L192F = rmsd_Apo
                rmsd_a4_F280Y_AD = rmsd_AD
                rmsd_a4_F280Y_BBR = rmsd_BBR
            if i == 3:
                rmsd_a4_E276F = rmsd_Apo
                rmsd_a4_L195F_AD = rmsd_AD
                rmsd_a4_L195F_BBR = rmsd_BBR
            if i == 4:
                rmsd_a4_F280Y = rmsd_Apo
                rmsd_a4_F196A_AD = rmsd_AD
                rmsd_a4_F196A_BBR = rmsd_BBR
                rmsd_a4_F196A_AD_rapo = rmsd_AD_rapo
                rmsd_a4_F196A_BBR_rapo = rmsd_BBR_rapo

            if i == 5:
                rmsd_a4_L195F = rmsd_Apo
                rmsd_a4_V287T_AD = rmsd_AD
                rmsd_a4_V287T_BBR = rmsd_BBR
            if i == 6:
                rmsd_a4_F196A = rmsd_Apo
                rmsd_a4_F196A_rapo = rmsd_Apo_rapo
            if i == 7:
                rmsd_a4_V287T = rmsd_Apo

        if sections[j] == 'a6_bot':
            if i == 0:
                rmsd_a6_bot_1sug = rmsd
            if i == 1:
                rmsd_a6_bot_apo = rmsd
            if i == 2:
                rmsd_a6_bot_L192F = rmsd
            if i == 3:
                rmsd_a6_bot_E276F = rmsd
            if i == 4:
                rmsd_a6_bot_F280Y = rmsd
            if i == 5:
                rmsd_a6_bot_L195F = rmsd
            if i == 6:
                rmsd_a6_bot_F196A = rmsd
            if i == 7:
                rmsd_a6_bot_V287T = rmsd

        if sections[j] == 'a6':
            if i == 0:
                rmsd_a6_1sug = rmsd
                rmsd_a6_L192F_AD = rmsd_AD
                rmsd_a6_L192F_BBR = rmsd_BBR
            if i == 1:
                rmsd_a6_apo = rmsd
                rmsd_a6_E276F_AD = rmsd_AD
                rmsd_a6_E276F_BBR = rmsd_BBR
            if i == 2:
                rmsd_a6_L192F = rmsd
                rmsd_a6_F280Y_AD = rmsd_AD
                rmsd_a6_F280Y_BBR = rmsd_BBR
            if i == 3:
                rmsd_a6_E276F = rmsd
                rmsd_a6_L195F_AD = rmsd_AD
                rmsd_a6_L195F_BBR = rmsd_BBR
            if i == 4:
                rmsd_a6_F280Y = rmsd
                rmsd_a6_F196A_AD = rmsd_AD
                rmsd_a6_F196A_BBR = rmsd_BBR
            if i == 5:
                rmsd_a6_L195F = rmsd
                rmsd_a6_V287T_AD = rmsd_AD
                rmsd_a6_V287T_BBR = rmsd_BBR
            if i == 6:
                rmsd_a6_F196A = rmsd
            if i == 7:
                rmsd_a6_V287T = rmsd

        if sections[j] == 'CYS':
            if i == 0:
                rmsd_cys_1sug = rmsd
            if i == 1:
                rmsd_cys_apo = rmsd
            if i == 2:
                rmsd_cys_L192F = rmsd
            if i == 3:
                rmsd_cys_E276F = rmsd
            if i == 4:
                rmsd_cys_F280Y = rmsd
            if i == 5:
                rmsd_cys_L195F = rmsd
            if i == 6:
                rmsd_cys_F196A = rmsd
            if i == 7:
                rmsd_cys_V287T = rmsd

#Name Labels
Label = ['WT', 'L192F', 'E276F', 'F280Y', 'L195F', 'F196A', 'V287T']
Label_close = ['WT Close', 'WT Open', 'L192F', 'E276F', 'F280Y', 'L195F', 'F196A', 'V287T']
Labels_mut = ['L192F', 'E276F', 'F280Y', 'L195F', 'F196A', 'V287T']

#Plot all compared to WT Open
for i in range(len(sections)):
    plot_compare(RMSD_mean, RMSD_err, Label, sections, i, ref[0])
    plot_compare(RMSD_mean_close, RMSD_err_close, Label_close, sections, i, ref[1])

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

RMSD_mean_mut = np.zeros((len(Labels_mut), len(sections))) #Mean for reference open
RMSD_err_mut = np.zeros((len(Labels_mut), len(sections))) #SEM for reference open

#Plot self and two references
for i in range(2, 5):
    for j in range(len(sections)):
        #Load Data
        for k in range(len(Labels_mut)):
            n = k + 1
            rmsd = load_rmsf(file_path[n], sections[j], ref[i])
            RMSD_mean_mut[k][j] = np.mean(rmsd)
            RMSD_err_mut[k][j] = stats.sem(rmsd)
        plot_compare(RMSD_mean_mut, RMSD_err_mut, Labels_mut, sections, j, ref[i])

#Plot Kernel DEnsity Estimate Plot
#Compare a3_top for L192F, E276F, L195F, V287T
a3_top_Apo_open_df = pd.DataFrame({'Apo Open':rmsd_a3_top_apo})
a3_top_Apo_close_df = pd.DataFrame({'Apo Closed': rmsd_a3_top_1sug})
a3_top_L192F_df = pd.DataFrame({'L192F': rmsd_a3_top_L192F})
a3_top_L195F_df = pd.DataFrame({'L195F': rmsd_a3_top_L195F})
a3_top_F280Y_df = pd.DataFrame({'F280Y': rmsd_a3_top_F280Y})
a3_top_E276F_df = pd.DataFrame({'E276F': rmsd_a3_top_E276F})
a3_top_F196A_df = pd.DataFrame({'F196A': rmsd_a3_top_F196A})
a3_top_V287T_df = pd.DataFrame({'V287T': rmsd_a3_top_V287T})
a3_top_L192F_AD_df = pd.DataFrame({'L192F AD': rmsd_a3_top_L192F_AD})
a3_top_L195F_AD_df = pd.DataFrame({'L195F AD': rmsd_a3_top_L195F_AD})
a3_top_F280Y_AD_df = pd.DataFrame({'F280Y AD': rmsd_a3_top_F280Y_AD})
a3_top_E276F_AD_df = pd.DataFrame({'E276F AD': rmsd_a3_top_E276F_AD})
a3_top_F196A_AD_df = pd.DataFrame({'F196A AD': rmsd_a3_top_F196A_AD})
a3_top_V287T_AD_df = pd.DataFrame({'V287T AD': rmsd_a3_top_V287T_AD})
a3_top_L192F_BBR_df = pd.DataFrame({'L192F BBR': rmsd_a3_top_L192F_BBR})
a3_top_L195F_BBR_df = pd.DataFrame({'L195F BBR': rmsd_a3_top_L195F_BBR})
a3_top_F280Y_BBR_df = pd.DataFrame({'F280Y BBR': rmsd_a3_top_F280Y_BBR})
a3_top_E276F_BBR_df = pd.DataFrame({'E276F BBR': rmsd_a3_top_E276F_BBR})
a3_top_F196A_BBR_df = pd.DataFrame({'F196A BBR': rmsd_a3_top_F196A_BBR})
a3_top_V287T_BBR_df = pd.DataFrame({'V287T BBR': rmsd_a3_top_V287T_BBR})

df = pd.concat([a3_top_Apo_open_df, a3_top_Apo_close_df, a3_top_L192F_df, a3_top_E276F_df, a3_top_V287T_df, a3_top_F196A_df, a3_top_F280Y_df, a3_top_L195F_df])

ax = plt.figure(figsize=(12, 6), frameon=False) # no visible frame
sns.kdeplot(data = df, fill=True, alpha=0.5, common_norm = True, common_grid = True)
plt.xlabel(r'RMSD($\AA$)')
plt.xlim(0, 2)
plt.ylabel(r'Normalized Density')
plt.title('Top of the a3 RMSD Compared to WT Closed')
plt.savefig('mutate_RMSD_a3_top_all.png')
plt.close()

df = pd.concat([a3_top_L192F_df, a3_top_E276F_df, a3_top_V287T_df, a3_top_F196A_df, a3_top_F280Y_df, a3_top_L195F_df])

sns.kdeplot(data = df, fill=True, alpha=0.5, common_norm = True, common_grid = True)
plt.axvline(x = RMSD_mean_close[0,7], color = 'r')
plt.axvline(x = RMSD_mean_close[1,7], color = 'b')
plt.xlabel(r'RMSD($\AA$)')
plt.xlim(0, 2)
plt.ylabel(r'Normalized Density')
plt.title('Top of the a3 RMSD Compared to WT Closed')
plt.savefig('mutate_RMSD_a3_top_mut_all.png')
plt.close()

df = pd.concat([a3_top_Apo_open_df, a3_top_Apo_close_df, a3_top_V287T_df, a3_top_F280Y_df])

sns.kdeplot(data = df, fill=True, alpha=0.5, common_norm = True, common_grid = True)
plt.xlabel(r'RMSD($\AA$)')
plt.xlim(0, 2)
plt.ylabel(r'Normalized Density')
plt.title('Top of the a3 RMSD Compared to WT Closed')
plt.savefig('mutate_RMSD_a3_top_mut_extr.png')
plt.close()

plot_kernel_cmpr_lig(a3_top_L192F_df, a3_top_L192F_AD_df, a3_top_L192F_BBR_df, 'L192F', sections[7], 7)
plot_kernel_cmpr_lig(a3_top_L195F_df, a3_top_L195F_AD_df, a3_top_L195F_BBR_df, 'L195F', sections[7], 7)
plot_kernel_cmpr_lig(a3_top_E276F_df, a3_top_E276F_AD_df, a3_top_E276F_BBR_df, 'E276F', sections[7], 7)
plot_kernel_cmpr_lig(a3_top_V287T_df, a3_top_V287T_AD_df, a3_top_V287T_BBR_df, 'V287T', sections[7], 7)

#Compare a4 for L192F, E276F, L195F, V287T
a4_Apo_open_df = pd.DataFrame({'Apo Open':rmsd_a4_apo})
a4_Apo_close_df = pd.DataFrame({'Apo Closed': rmsd_a4_1sug})
a4_L192F_df = pd.DataFrame({'L192F': rmsd_a4_L192F})
a4_L195F_df = pd.DataFrame({'L195F': rmsd_a4_L195F})
a4_F280Y_df = pd.DataFrame({'F280Y': rmsd_a4_F280Y})
a4_E276F_df = pd.DataFrame({'E276F': rmsd_a4_E276F})
a4_F196A_df = pd.DataFrame({'F196A': rmsd_a4_F196A})
a4_V287T_df = pd.DataFrame({'V287T': rmsd_a4_V287T})
a4_L192F_AD_df = pd.DataFrame({'L192F AD': rmsd_a4_L192F_AD})
a4_L195F_AD_df = pd.DataFrame({'L195F AD': rmsd_a4_L195F_AD})
a4_F280Y_AD_df = pd.DataFrame({'F280Y AD': rmsd_a4_F280Y_AD})
a4_E276F_AD_df = pd.DataFrame({'E276F AD': rmsd_a4_E276F_AD})
a4_F196A_AD_df = pd.DataFrame({'F196A AD': rmsd_a4_F196A_AD})
a4_V287T_AD_df = pd.DataFrame({'V287T AD': rmsd_a4_V287T_AD})
a4_L192F_BBR_df = pd.DataFrame({'L192F BBR': rmsd_a4_L192F_BBR})
a4_L195F_BBR_df = pd.DataFrame({'L195F BBR': rmsd_a4_L195F_BBR})
a4_F280Y_BBR_df = pd.DataFrame({'F280Y BBR': rmsd_a4_F280Y_BBR})
a4_E276F_BBR_df = pd.DataFrame({'E276F BBR': rmsd_a4_E276F_BBR})
a4_F196A_BBR_df = pd.DataFrame({'F196A BBR': rmsd_a4_F196A_BBR})
a4_V287T_BBR_df = pd.DataFrame({'V287T BBR': rmsd_a4_V287T_BBR})

df = pd.concat([a4_Apo_open_df, a4_Apo_close_df, a4_L192F_df, a4_E276F_df, a4_V287T_df, a4_F196A_df, a4_F280Y_df, a4_L195F_df])

ax = plt.figure(figsize=(12, 6), frameon=False) # no visible frame
sns.kdeplot(data = df, fill=True, alpha=0.5, common_norm = True, common_grid = True)
plt.xlabel(r'RMSD($\AA$)')
plt.xlim(0, 2)
plt.ylabel(r'Normalized Density')
plt.title(r'$\alpha$-4 RMSD Compared to WT Closed')
plt.savefig('mutate_RMSD_a4_all.png')
plt.close()

df = pd.concat([a4_Apo_open_df, a4_Apo_close_df, a4_V287T_df, a4_F196A_df, a4_F280Y_df])

ax = plt.figure(figsize=(12, 6), frameon=False) # no visible frame
sns.kdeplot(data = df, fill=True, alpha=0.5, common_norm = True, common_grid = True)
plt.xlabel(r'RMSD($\AA$)')
plt.xlim(0, 2)
plt.ylabel(r'Normalized Density')
plt.title(r'$\alpha$-4 RMSD Compared to WT Closed')
plt.savefig('mutate_RMSD_a4.png')
plt.close()

plot_kernel_cmpr_lig(a4_F196A_df, a4_F196A_AD_df, a4_F196A_BBR_df, 'F196A', sections[8], 8)
plot_kernel_cmpr_lig(a4_F280Y_df, a4_F280Y_AD_df, a4_F280Y_BBR_df, 'F280Y', sections[8], 8)

a4_Apo_open_df = pd.DataFrame({'Apo Open':rmsd_a4_apo_rapo})
a4_Apo_close_df = pd.DataFrame({'Apo Closed': rmsd_a4_1sug_rapo})
a4_F196A_df = pd.DataFrame({'F196A': rmsd_a4_F196A_rapo})
a4_F196A_AD_df = pd.DataFrame({'F196A AD': rmsd_a4_F196A_AD_rapo})
a4_F196A_BBR_df = pd.DataFrame({'F196A BBR': rmsd_a4_F196A_BBR_rapo})

df = pd.concat([a4_Apo_open_df, a4_Apo_close_df, a4_F196A_df, a4_F196A_AD_df, a4_F196A_BBR_df])

ax = plt.figure(figsize=(12, 6), frameon=False) # no visible frame
sns.kdeplot(data = df, fill=True, alpha=0.5, common_norm = True, common_grid = True)
plt.xlabel(r'RMSD($\AA$)')
plt.ylabel(r'Normalized Density')
plt.title(r'$\alpha$-4 RMSD Compared to Apo F196A')
plt.savefig('mutate_RMSD_a4_ref_F196A.png')
plt.close()

#Just CYS215
rmsd_cys = [rmsd_cys_1sug, rmsd_cys_apo, rmsd_cys_L192F, rmsd_cys_E276F, rmsd_cys_V287T, rmsd_cys_F196A, rmsd_cys_F280Y, rmsd_cys_L195F]
ax = plt.figure(figsize=(12, 6), frameon=False) # no visible frame
sns.kdeplot(data = rmsd_cys, fill=True, alpha=0.5)
plt.title('CYS215 RMSD Compared to WT Closed')
plt.savefig('mutate_RMSD_cys_all.png')
plt.close()

rmsd_cys = [rmsd_cys_1sug, rmsd_cys_apo, rmsd_cys_F196A]
ax = plt.figure(figsize=(12, 6), frameon=False) # no visible frame
sns.kdeplot(data = rmsd_cys, fill=True, alpha=0.5)
plt.title('CYS215 RMSD Compared to WT Closed')
plt.savefig('mutate_RMSD_cys_F196A.png')
plt.close()

#Determine p-values for each of the sections of focus
file_p = open('p_values_mut.txt', 'w')
p = np.zeros((4, 6))
st, p[0,0] = stats.ttest_ind(rmsd_a3_top_apo, rmsd_a3_top_L192F, equal_var = False) #Welch's t-test
st, p[0,1] = stats.ttest_ind(rmsd_a3_top_apo, rmsd_a3_top_E276F, equal_var = False) #Welch's t-test
st, p[0,2] = stats.ttest_ind(rmsd_a3_top_apo, rmsd_a3_top_F280Y, equal_var = False) #Welch's t-test
st, p[0,3] = stats.ttest_ind(rmsd_a3_top_apo, rmsd_a3_top_L195F, equal_var = False) #Welch's t-test
st, p[0,4] = stats.ttest_ind(rmsd_a3_top_apo, rmsd_a3_top_F196A, equal_var = False) #Welch's t-test
st, p[0,5] = stats.ttest_ind(rmsd_a3_top_apo, rmsd_a3_top_V287T, equal_var = False) #Welch's t-test

st, p[1,0] = stats.ttest_ind(rmsd_a4_apo, rmsd_a4_L192F, equal_var = False) #Welch's t-test
st, p[1,1] = stats.ttest_ind(rmsd_a4_apo, rmsd_a4_E276F, equal_var = False) #Welch's t-test
st, p[1,2] = stats.ttest_ind(rmsd_a4_apo, rmsd_a4_F280Y, equal_var = False) #Welch's t-test
st, p[1,3] = stats.ttest_ind(rmsd_a4_apo, rmsd_a4_L195F, equal_var = False) #Welch's t-test
st, p[1,4] = stats.ttest_ind(rmsd_a4_apo, rmsd_a4_F196A, equal_var = False) #Welch's t-test
st, p[1,5] = stats.ttest_ind(rmsd_a4_apo, rmsd_a4_V287T, equal_var = False) #Welch's t-test

st, p[2,0] = stats.ttest_ind(rmsd_a6_bot_apo, rmsd_a6_bot_L192F, equal_var = False) #Welch's t-test
st, p[2,1] = stats.ttest_ind(rmsd_a6_bot_apo, rmsd_a6_bot_E276F, equal_var = False) #Welch's t-test
st, p[2,2] = stats.ttest_ind(rmsd_a6_bot_apo, rmsd_a6_bot_F280Y, equal_var = False) #Welch's t-test
st, p[2,3] = stats.ttest_ind(rmsd_a6_bot_apo, rmsd_a6_bot_L195F, equal_var = False) #Welch's t-test
st, p[2,4] = stats.ttest_ind(rmsd_a6_bot_apo, rmsd_a6_bot_F196A, equal_var = False) #Welch's t-test
st, p[2,5] = stats.ttest_ind(rmsd_a6_bot_apo, rmsd_a6_bot_V287T, equal_var = False) #Welch's t-test

st, p[3,0] = stats.ttest_ind(rmsd_cys_apo, rmsd_cys_L192F, equal_var = False) #Welch's t-test
st, p[3,1] = stats.ttest_ind(rmsd_cys_apo, rmsd_cys_E276F, equal_var = False) #Welch's t-test
st, p[3,2] = stats.ttest_ind(rmsd_cys_apo, rmsd_cys_F280Y, equal_var = False) #Welch's t-test
st, p[3,3] = stats.ttest_ind(rmsd_cys_apo, rmsd_cys_L195F, equal_var = False) #Welch's t-test
st, p[3,4] = stats.ttest_ind(rmsd_cys_apo, rmsd_cys_F196A, equal_var = False) #Welch's t-test
st, p[3,5] = stats.ttest_ind(rmsd_cys_apo, rmsd_cys_V287T, equal_var = False) #Welch's t-test

sections_mini = ['a3_top', 'a4', 'a6_bot', 'CYS215']
for i in range(len(sections_mini)):
    file_p.write(str(sections_mini[i]) + '\n')
    for j in range(len(Labels_mut)):
        file_p.write(Labels_mut[j] + ': ' + str(p[i,j]) + '\n')

