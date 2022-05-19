#Import Necessary Packages
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import seaborn as sns
import ruptures as rpt
from statistics import stdev
import pandas as pd

def load_rms(path, sect, ref):
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

def rmsd_sect(sect, file_path_close, file_path_close_AD, file_path_close_BBR, ref, n):
    rmsd_1sug = load_rms(file_path_close[0], sect, ref[n])
    rmsd_apo = load_rms(file_path_close[1], sect, ref[n])
    rmsd_L192F = load_rms(file_path_close[2], sect, ref[n])
    rmsd_E276F = load_rms(file_path_close[3], sect, ref[n])
    rmsd_F280Y = load_rms(file_path_close[4], sect, ref[n])
    rmsd_L195F = load_rms(file_path_close[5], sect, ref[n])
    rmsd_F196A = load_rms(file_path_close[6], sect, ref[n])
    rmsd_V287T = load_rms(file_path_close[7], sect, ref[n])

    rmsd_L192F_AD = load_rms(file_path_close_AD[0], sect, ref[n])
    rmsd_L192F_BBR = load_rms(file_path_close_BBR[0], sect, ref[n])
    rmsd_E276F_AD = load_rms(file_path_close_AD[1], sect, ref[n])
    rmsd_E276F_BBR = load_rms(file_path_close_BBR[1], sect, ref[n])
    rmsd_F280Y_AD = load_rms(file_path_close_AD[2], sect, ref[n])
    rmsd_F280Y_BBR = load_rms(file_path_close_BBR[2], sect, ref[n])
    rmsd_L195F_AD = load_rms(file_path_close_AD[3], sect, ref[n])
    rmsd_L195F_BBR = load_rms(file_path_close_BBR[3], sect, ref[n])
    rmsd_F196A_AD = load_rms(file_path_close_AD[4], sect, ref[n])
    rmsd_F196A_BBR = load_rms(file_path_close_BBR[4], sect, ref[n])
    rmsd_V287T_AD = load_rms(file_path_close_AD[4], sect, ref[n])
    rmsd_V287T_BBR = load_rms(file_path_close_BBR[4], sect, ref[n])
    
    return rmsd_1sug, rmsd_apo, rmsd_L192F, rmsd_E276F, rmsd_F280Y, rmsd_L195F, rmsd_F196A, rmsd_V287T, rmsd_L192F_AD, rmsd_E276F_AD, rmsd_F280Y_AD, rmsd_L195F_AD, rmsd_F196A_AD, rmsd_V287T_AD, rmsd_L192F_BBR, rmsd_E276F_BBR, rmsd_F280Y_BBR, rmsd_L195F_BBR, rmsd_F196A_BBR, rmsd_V287T_BBR

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
rmsd_a3_1sug, rmsd_a3_apo, rmsd_a3_L192F, rmsd_a3_E276F, rmsd_a3_F280Y, rmsd_a3_L195F, rmsd_a3_F196A, rmsd_a3_V287T, rmsd_a3_L192F_AD, rmsd_a3_E276F_AD, rmsd_a3_F280Y_AD, rmsd_a3_L195F_AD, rmsd_a3_F196A_AD, rmsd_a3_V287T_AD, rmsd_a3_L192F_BBR, rmsd_a3_E276F_BBR, rmsd_a3_F280Y_BBR, rmsd_a3_L195F_BBR, rmsd_a3_F196A_BBR, rmsd_a3_V287T_BBR = rmsd_sect('a3', file_path_close, file_path_close_AD, file_path_close_BBR, ref, 1)

rmsd_a3_top_1sug, rmsd_a3_top_apo, rmsd_a3_top_L192F, rmsd_a3_top_E276F, rmsd_a3_top_F280Y, rmsd_a3_top_L195F, rmsd_a3_top_F196A, rmsd_a3_top_V287T, rmsd_a3_top_L192F_AD, rmsd_a3_top_E276F_AD, rmsd_a3_top_F280Y_AD, rmsd_a3_top_L195F_AD, rmsd_a3_top_F196A_AD, rmsd_a3_top_V287T_AD, rmsd_a3_top_L192F_BBR, rmsd_a3_top_E276F_BBR, rmsd_a3_top_F280Y_BBR, rmsd_a3_top_L195F_BBR, rmsd_a3_top_F196A_BBR, rmsd_a3_top_V287T_BBR = rmsd_sect('a3_top', file_path_close, file_path_close_AD, file_path_close_BBR, ref, 1)

rmsd_a4_1sug, rmsd_a4_apo, rmsd_a4_L192F, rmsd_a4_E276F, rmsd_a4_F280Y, rmsd_a4_L195F, rmsd_a4_F196A, rmsd_a4_V287T, rmsd_a4_L192F_AD, rmsd_a4_E276F_AD, rmsd_a4_F280Y_AD, rmsd_a4_L195F_AD, rmsd_a4_F196A_AD, rmsd_a4_V287T_AD, rmsd_a4_L192F_BBR, rmsd_a4_E276F_BBR, rmsd_a4_F280Y_BBR, rmsd_a4_L195F_BBR, rmsd_a4_F196A_BBR, rmsd_a4_V287T_BBR = rmsd_sect('a4', file_path_close, file_path_close_AD, file_path_close_BBR, ref, 1)

rmsd_a6_1sug, rmsd_a6_apo, rmsd_a6_L192F, rmsd_a6_E276F, rmsd_a6_F280Y, rmsd_a6_L195F, rmsd_a6_F196A, rmsd_a6_V287T, rmsd_a6_L192F_AD, rmsd_a6_E276F_AD, rmsd_a6_F280Y_AD, rmsd_a6_L195F_AD, rmsd_a6_F196A_AD, rmsd_a6_V287T_AD, rmsd_a6_L192F_BBR, rmsd_a6_E276F_BBR, rmsd_a6_F280Y_BBR, rmsd_a6_L195F_BBR, rmsd_a6_F196A_BBR, rmsd_a6_V287T_BBR = rmsd_sect('a6', file_path_close, file_path_close_AD, file_path_close_BBR, ref, 1)

rmsd_a6_bot_1sug, rmsd_a6_bot_apo, rmsd_a6_bot_L192F, rmsd_a6_bot_E276F, rmsd_a6_bot_F280Y, rmsd_a6_bot_L195F, rmsd_a6_bot_F196A, rmsd_a6_bot_V287T, rmsd_a6_bot_L192F_AD, rmsd_a6_bot_E276F_AD, rmsd_a6_bot_F280Y_AD, rmsd_a6_bot_L195F_AD, rmsd_a6_bot_F196A_AD, rmsd_a6_bot_V287T_AD, rmsd_a6_bot_L192F_BBR, rmsd_a6_bot_E276F_BBR, rmsd_a6_bot_F280Y_BBR, rmsd_a6_bot_L195F_BBR, rmsd_a6_bot_F196A_BBR, rmsd_a6_bot_V287T_BBR = rmsd_sect('a6_bot', file_path_close, file_path_close_AD, file_path_close_BBR, ref, 1)

rmsd_CYS_1sug, rmsd_CYS_apo, rmsd_CYS_L192F, rmsd_CYS_E276F, rmsd_CYS_F280Y, rmsd_CYS_L195F, rmsd_CYS_F196A, rmsd_CYS_V287T, rmsd_CYS_L192F_AD, rmsd_CYS_E276F_AD, rmsd_CYS_F280Y_AD, rmsd_CYS_L195F_AD, rmsd_CYS_F196A_AD, rmsd_CYS_V287T_AD, rmsd_CYS_L192F_BBR, rmsd_CYS_E276F_BBR, rmsd_CYS_F280Y_BBR, rmsd_CYS_L195F_BBR, rmsd_CYS_F196A_BBR, rmsd_CYS_V287T_BBR = rmsd_sect('CYS', file_path_close, file_path_close_AD, file_path_close_BBR, ref, 1)

rmsd_beg_1sug, rmsd_beg_apo, rmsd_beg_L192F, rmsd_beg_E276F, rmsd_beg_F280Y, rmsd_beg_L195F, rmsd_beg_F196A, rmsd_beg_V287T, rmsd_beg_L192F_AD, rmsd_beg_E276F_AD, rmsd_beg_F280Y_AD, rmsd_beg_L195F_AD, rmsd_beg_F196A_AD, rmsd_beg_V287T_AD, rmsd_beg_L192F_BBR, rmsd_beg_E276F_BBR, rmsd_beg_F280Y_BBR, rmsd_beg_L195F_BBR, rmsd_beg_F196A_BBR, rmsd_beg_V287T_BBR = rmsd_sect('beg', file_path_close, file_path_close_AD, file_path_close_BBR, ref, 1)

for i in range(len(file_path_close)):
    for j in range(len(sections)):
        #Load Data for reference open
        rmsd_Apo = load_rms(file_path_close[i], sections[j], ref[1])

        #Mean and SEM for each trajectory
        RMSD_mean_close[i][j] = np.mean(rmsd_Apo)
        RMSD_err_close[i][j] = stats.sem(rmsd_Apo)

for i in range(len(file_path)):
    #Load Data for reference open
    rmsd = load_rms(file_path[i], sections[j], ref[0])
    #Mean and SEM for each trajectory
    RMSD_mean[i][j] = np.mean(rmsd)
    RMSD_err[i][j] = stats.sem(rmsd)

for i in range(len(file_path_close_AD)):
    #Load Data for reference open
    rmsd_AD = load_rms(file_path_close_AD[i], sections[j], ref[1])
    RMSD_mean_close_AD[i][j] = np.mean(rmsd_AD)
    RMSD_err_close_AD[i][j] = stats.sem(rmsd_AD)

    rmsd_BBR = load_rms(file_path_close_BBR[i], sections[j], ref[1])
    RMSD_mean_close_BBR[i][j] = np.mean(rmsd_BBR)
    RMSD_err_close_BBR[i][j] = stats.sem(rmsd_BBR)

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

RMSD_mean_mut = np.zeros((len(Label_close), len(sections))) #Mean for reference open
RMSD_err_mut = np.zeros((len(Label_close), len(sections))) #SEM for reference open

#Plot self and two references
for i in [0, 2, 3, 4]:
    for j in range(len(sections)):
        #Load Data
        RMSD_mean_mut[0][j] = RMSD_mean_close[0][j]
        RMSD_err_mut[0][j] = RMSD_err_close[0][j]
        for k in range(1, len(Label_close)):
            rmsd = load_rms(file_path_close[k], sections[j], ref[i])
            RMSD_mean_mut[k][j] = np.mean(rmsd)
            RMSD_err_mut[k][j] = stats.sem(rmsd)
        plot_compare(RMSD_mean_mut, RMSD_err_mut, Label_close, sections, j, ref[i])

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

sns.kdeplot(data = df, fill=True, alpha=0.5, common_grid = True)
plt.xlabel(r'RMSD($\AA$)', fontsize = 12)
plt.xlim(0, 2)
plt.xticks(fontsize = 10)
plt.yticks(fontsize = 10)
plt.ylabel(r'Normalized Density', fontsize = 12)
plt.title(r'Top of the $\alpha$3 RMSD Compared to WT Closed', fontsize = 14)
plt.savefig('mutate_RMSD_a3_top_all.png')
plt.close()

df = pd.concat([a3_top_L192F_df, a3_top_E276F_df, a3_top_V287T_df, a3_top_F196A_df, a3_top_F280Y_df, a3_top_L195F_df])

sns.kdeplot(data = df, fill=True, alpha=0.5, common_norm = True, common_grid = True)
plt.xlabel(r'RMSD($\AA$)', fontsize = 12)
plt.xlim(0, 2)
plt.xticks(fontsize = 10)
plt.yticks(fontsize = 10)
plt.ylabel(r'Normalized Density', fontsize = 12)
plt.title(r'Top of the $\alpha$3 RMSD Compared to WT Closed', fontsize = 14)
plt.savefig('mutate_RMSD_a3_top_mut_all.png')
plt.close()

df = pd.concat([a3_top_Apo_open_df, a3_top_Apo_close_df, a3_top_V287T_df, a3_top_F280Y_df])

sns.kdeplot(data = df, fill=True, alpha=0.5, common_norm = True, common_grid = True)
plt.xlabel(r'RMSD($\AA$)', fontsize = 12)
plt.xlim(0, 2)
plt.xticks(fontsize = 10)
plt.yticks(fontsize = 10)
plt.ylabel(r'Normalized Density', fontsize = 12)
plt.title(r'Top of the $\alpha$3 RMSD Compared to WT Closed', fontsize = 14)
plt.savefig('mutate_RMSD_a3_top_mut_extr.png')
plt.close()

plot_kernel_cmpr_lig(a3_top_L192F_df, a3_top_L192F_AD_df, a3_top_L192F_BBR_df, 'L192F', sections[7], 7)
plot_kernel_cmpr_lig(a3_top_L195F_df, a3_top_L195F_AD_df, a3_top_L195F_BBR_df, 'L195F', sections[7], 7)
plot_kernel_cmpr_lig(a3_top_E276F_df, a3_top_E276F_AD_df, a3_top_E276F_BBR_df, 'E276F', sections[7], 7)
plot_kernel_cmpr_lig(a3_top_V287T_df, a3_top_V287T_AD_df, a3_top_V287T_BBR_df, 'V287T', sections[7], 7)

#Compare a3_top for L192F, E276F, L195F, V287T
a3_Apo_open_df = pd.DataFrame({'Apo Open':rmsd_a3_apo})
a3_Apo_close_df = pd.DataFrame({'Apo Closed': rmsd_a3_1sug})
a3_L192F_df = pd.DataFrame({'L192F': rmsd_a3_L192F})
a3_L195F_df = pd.DataFrame({'L195F': rmsd_a3_L195F})
a3_F280Y_df = pd.DataFrame({'F280Y': rmsd_a3_F280Y})
a3_E276F_df = pd.DataFrame({'E276F': rmsd_a3_E276F})
a3_F196A_df = pd.DataFrame({'F196A': rmsd_a3_F196A})
a3_V287T_df = pd.DataFrame({'V287T': rmsd_a3_V287T})
a3_L192F_AD_df = pd.DataFrame({'L192F AD': rmsd_a3_L192F_AD})
a3_L195F_AD_df = pd.DataFrame({'L195F AD': rmsd_a3_L195F_AD})
a3_F280Y_AD_df = pd.DataFrame({'F280Y AD': rmsd_a3_F280Y_AD})
a3_E276F_AD_df = pd.DataFrame({'E276F AD': rmsd_a3_E276F_AD})
a3_F196A_AD_df = pd.DataFrame({'F196A AD': rmsd_a3_F196A_AD})
a3_V287T_AD_df = pd.DataFrame({'V287T AD': rmsd_a3_V287T_AD})
a3_L192F_BBR_df = pd.DataFrame({'L192F BBR': rmsd_a3_L192F_BBR})
a3_L195F_BBR_df = pd.DataFrame({'L195F BBR': rmsd_a3_L195F_BBR})
a3_F280Y_BBR_df = pd.DataFrame({'F280Y BBR': rmsd_a3_F280Y_BBR})
a3_E276F_BBR_df = pd.DataFrame({'E276F BBR': rmsd_a3_E276F_BBR})
a3_F196A_BBR_df = pd.DataFrame({'F196A BBR': rmsd_a3_F196A_BBR})
a3_V287T_BBR_df = pd.DataFrame({'V287T BBR': rmsd_a3_V287T_BBR})

df = pd.concat([a3_Apo_open_df, a3_Apo_close_df, a3_L192F_df, a3_E276F_df, a3_V287T_df, a3_F196A_df, a3_F280Y_df, a3_L195F_df])

ax = plt.figure(figsize=(12, 6), frameon=False) # no visible frame
sns.kdeplot(data = df, fill=True, alpha=0.5, common_grid = True)
plt.xlabel(r'RMSD($\AA$)')
plt.xlim(0, 2)
plt.ylabel(r'Normalized Density')
plt.title(r'$\alpha$-3 RMSD Compared to WT Closed')
plt.savefig('mutate_RMSD_a3_all.png')
plt.close()

df = pd.concat([a3_L192F_df, a3_E276F_df, a3_V287T_df, a3_F196A_df, a3_F280Y_df, a3_L195F_df])

sns.kdeplot(data = df, fill=True, alpha=0.5, common_norm = True, common_grid = True)
plt.xlabel(r'RMSD($\AA$)', fontsize = 12)
plt.xlim(0, 2)
plt.xticks(fontsize = 10)
plt.yticks(fontsize = 10)
plt.ylabel(r'Normalized Density', fontsize = 12)
plt.title(r'$\alpha$-3 RMSD Compared to WT Closed', fontsize = 14)
plt.savefig('mutate_RMSD_a3_mut_all.png')
plt.close()

df = pd.concat([a3_Apo_open_df, a3_Apo_close_df, a3_V287T_df, a3_F280Y_df])

sns.kdeplot(data = df, fill=True, alpha=0.5, common_norm = True, common_grid = True)
plt.xlabel(r'RMSD($\AA$)')
plt.xlim(0, 2)
plt.ylabel(r'Normalized Density')
plt.title(r'$\alpha$-3 RMSD Compared to WT Closed')
plt.savefig('mutate_RMSD_a3_mut_extr.png')
plt.close()

plot_kernel_cmpr_lig(a3_L192F_df, a3_L192F_AD_df, a3_L192F_BBR_df, 'L192F', sections[6], 6)
plot_kernel_cmpr_lig(a3_L195F_df, a3_L195F_AD_df, a3_L195F_BBR_df, 'L195F', sections[6], 6)
plot_kernel_cmpr_lig(a3_E276F_df, a3_E276F_AD_df, a3_E276F_BBR_df, 'E276F', sections[6], 6)
plot_kernel_cmpr_lig(a3_V287T_df, a3_V287T_AD_df, a3_V287T_BBR_df, 'V287T', sections[6], 6)

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

sns.kdeplot(data = df, fill=True, alpha=0.5, common_norm = True, common_grid = True)
plt.xlabel(r'RMSD($\AA$)', fontsize = 12)
plt.xticks(fontsize = 10)
plt.yticks(fontsize = 10)
plt.xlim(0, 1.5)
plt.ylabel(r'Normalized Density', fontsize = 12)
plt.title(r'$\alpha$-4 RMSD Compared to WT Closed', fontsize = 14)
plt.savefig('mutate_RMSD_a4_all.png')
plt.close()

df = pd.concat([a4_L192F_df, a4_E276F_df, a4_V287T_df, a4_F196A_df, a4_F280Y_df, a4_L195F_df])

sns.kdeplot(data = df, fill=True, alpha=0.5, common_norm = True, common_grid = True)
plt.xlabel(r'RMSD($\AA$)', fontsize = 12)
plt.xticks(fontsize = 10)
plt.yticks(fontsize = 10)
plt.xlim(0, 1.5)
plt.ylabel(r'Normalized Density', fontsize = 12)
plt.title(r'$\alpha$-4 RMSD Compared to WT Closed', fontsize = 14)
plt.savefig('mutate_RMSD_a4_mut_all.png')
plt.close()

df = pd.concat([a4_Apo_open_df, a4_Apo_close_df, a4_V287T_df, a4_F196A_df, a4_F280Y_df])

ax = plt.figure(figsize=(12, 6), frameon=False) # no visible frame
sns.kdeplot(data = df, fill=True, alpha=0.5, common_norm = True, common_grid = True)
plt.xlabel(r'RMSD($\AA$)')
plt.xlim(0, 1.5)
plt.ylabel(r'Normalized Density')
plt.title(r'$\alpha$-4 RMSD Compared to WT Closed')
plt.savefig('mutate_RMSD_a4.png')
plt.close()

plot_kernel_cmpr_lig(a4_F196A_df, a4_F196A_AD_df, a4_F196A_BBR_df, 'F196A', sections[8], 8)
plot_kernel_cmpr_lig(a4_F280Y_df, a4_F280Y_AD_df, a4_F280Y_BBR_df, 'F280Y', sections[8], 8)

#a4_Apo_open_df = pd.DataFrame({'Apo Open':rmsd_a4_apo_rapo})
#a4_Apo_close_df = pd.DataFrame({'Apo Closed': rmsd_a4_1sug_rapo})
#a4_F196A_df = pd.DataFrame({'F196A': rmsd_a4_F196A_rapo})
#a4_F196A_AD_df = pd.DataFrame({'F196A AD': rmsd_a4_F196A_AD_rapo})
#a4_F196A_BBR_df = pd.DataFrame({'F196A BBR': rmsd_a4_F196A_BBR_rapo})

#df = pd.concat([a4_Apo_open_df, a4_Apo_close_df, a4_F196A_df, a4_F196A_AD_df, a4_F196A_BBR_df])

#ax = plt.figure(figsize=(12, 6), frameon=False) # no visible frame
#sns.kdeplot(data = df, fill=True, alpha=0.5, common_norm = True, common_grid = True)
#plt.xlabel(r'RMSD($\AA$)')
#plt.ylabel(r'Normalized Density')
#plt.title(r'$\alpha$-4 RMSD Compared to Apo F196A')
#plt.savefig('mutate_RMSD_a4_ref_F196A.png')
#plt.close()

#a6 comparison
a6_Apo_open_df = pd.DataFrame({'Apo Open':rmsd_a6_apo})
a6_Apo_close_df = pd.DataFrame({'Apo Closed': rmsd_a6_1sug})
a6_L192F_df = pd.DataFrame({'L192F': rmsd_a6_L192F})
a6_L195F_df = pd.DataFrame({'L195F': rmsd_a6_L195F})
a6_F280Y_df = pd.DataFrame({'F280Y': rmsd_a6_F280Y})
a6_E276F_df = pd.DataFrame({'E276F': rmsd_a6_E276F})
a6_F196A_df = pd.DataFrame({'F196A': rmsd_a6_F196A})
a6_V287T_df = pd.DataFrame({'V287T': rmsd_a6_V287T})
a6_L192F_AD_df = pd.DataFrame({'L192F AD': rmsd_a6_L192F_AD})
a6_L195F_AD_df = pd.DataFrame({'L195F AD': rmsd_a6_L195F_AD})
a6_F280Y_AD_df = pd.DataFrame({'F280Y AD': rmsd_a6_F280Y_AD})
a6_E276F_AD_df = pd.DataFrame({'E276F AD': rmsd_a6_E276F_AD})
a6_F196A_AD_df = pd.DataFrame({'F196A AD': rmsd_a6_F196A_AD})
a6_V287T_AD_df = pd.DataFrame({'V287T AD': rmsd_a6_V287T_AD})
a6_L192F_BBR_df = pd.DataFrame({'L192F BBR': rmsd_a6_L192F_BBR})
a6_L195F_BBR_df = pd.DataFrame({'L195F BBR': rmsd_a6_L195F_BBR})
a6_F280Y_BBR_df = pd.DataFrame({'F280Y BBR': rmsd_a6_F280Y_BBR})
a6_E276F_BBR_df = pd.DataFrame({'E276F BBR': rmsd_a6_E276F_BBR})
a6_F196A_BBR_df = pd.DataFrame({'F196A BBR': rmsd_a6_F196A_BBR})
a6_V287T_BBR_df = pd.DataFrame({'V287T BBR': rmsd_a6_V287T_BBR})

df = pd.concat([a6_Apo_open_df, a6_Apo_close_df, a6_L192F_df, a6_E276F_df, a6_V287T_df, a6_F196A_df, a6_F280Y_df, a6_L195F_df])

sns.kdeplot(data = df, fill=True, alpha=0.5, common_grid = True)
plt.xlabel(r'RMSD($\AA$)', fontsize = 12)
plt.xlim(0, 2)
plt.xticks(fontsize = 10)
plt.yticks(fontsize = 10)
plt.ylabel(r'Normalized Density', fontsize = 12)
plt.title(r'$\alpha$-6 RMSD Compared to WT Closed', fontsize = 14)
plt.savefig('mutate_RMSD_a6_all.png')
plt.close()

df = pd.concat([a6_L192F_df, a6_E276F_df, a6_V287T_df, a6_F196A_df, a6_F280Y_df, a6_L195F_df])

sns.kdeplot(data = df, fill=True, alpha=0.5, common_norm = True, common_grid = True)
plt.xlabel(r'RMSD($\AA$)', fontsize = 12)
plt.xlim(0, 2)
plt.xticks(fontsize = 10)
plt.yticks(fontsize = 10)
plt.ylabel(r'Normalized Density', fontsize = 12)
plt.title(r'$\alpha$-6 RMSD Compared to WT Closed', fontsize = 14)
plt.savefig('mutate_RMSD_a6_mut_all.png')
plt.close()

plot_kernel_cmpr_lig(a6_L192F_df, a6_L192F_AD_df, a6_L192F_BBR_df, 'L192F', sections[11], 11)
plot_kernel_cmpr_lig(a6_L195F_df, a6_L195F_AD_df, a6_L195F_BBR_df, 'L195F', sections[11], 11)
plot_kernel_cmpr_lig(a6_E276F_df, a6_E276F_AD_df, a6_E276F_BBR_df, 'E276F', sections[11], 11)
plot_kernel_cmpr_lig(a6_V287T_df, a6_V287T_AD_df, a6_V287T_BBR_df, 'V287T', sections[11], 11)

#Just CYS215
cys_Apo_open_df = pd.DataFrame({'Apo Open':rmsd_CYS_apo})
cys_Apo_close_df = pd.DataFrame({'Apo Closed': rmsd_CYS_1sug})
cys_L192F_df = pd.DataFrame({'L192F': rmsd_CYS_L192F})
cys_L195F_df = pd.DataFrame({'L195F': rmsd_CYS_L195F})
cys_F280Y_df = pd.DataFrame({'F280Y': rmsd_CYS_F280Y})
cys_E276F_df = pd.DataFrame({'E276F': rmsd_CYS_E276F})
cys_F196A_df = pd.DataFrame({'F196A': rmsd_CYS_F196A})
cys_V287T_df = pd.DataFrame({'V287T': rmsd_CYS_V287T})

df = pd.concat([cys_Apo_open_df, cys_Apo_close_df, cys_L192F_df, cys_E276F_df, cys_V287T_df, cys_F196A_df, cys_F280Y_df, cys_L195F_df])

sns.kdeplot(data = df, fill=True, alpha=0.5, common_grid = True)
plt.xlabel(r'RMSD($\AA$)', fontsize = 12)
plt.xlim(0, 1)
plt.xticks(fontsize = 10)
plt.yticks(fontsize = 10)
plt.ylabel(r'Normalized Density', fontsize = 12)
plt.title(r'CYS215 RMSD Compared to WT Closed', fontsize = 14)
plt.savefig('mutate_RMSD_cys_all.png')
plt.close()

df = pd.concat([cys_L192F_df, cys_E276F_df, cys_V287T_df, cys_F196A_df, cys_F280Y_df, cys_L195F_df])

sns.kdeplot(data = df, fill=True, alpha=0.5, common_grid = True)
plt.xlabel(r'RMSD($\AA$)', fontsize = 12)
plt.xlim(0, 1)
plt.xticks(fontsize = 10)
plt.yticks(fontsize = 10)
plt.ylabel(r'Normalized Density', fontsize = 12)
plt.title(r'CYS215 RMSD Compared to WT Closed', fontsize = 14)
plt.savefig('mutate_RMSD_cys_mut_all.png')
plt.close()

rmsd_cys = [rmsd_CYS_1sug, rmsd_CYS_apo, rmsd_CYS_F196A]
ax = plt.figure(figsize=(12, 6), frameon=False) # no visible frame
sns.kdeplot(data = rmsd_cys, fill=True, alpha=0.5)
plt.title('CYS215 RMSD Compared to WT Closed')
plt.savefig('mutate_RMSD_cys_F196A.png')
plt.close()

#BEG loop (L1)
beg_Apo_open_df = pd.DataFrame({'Apo Open':rmsd_beg_apo})
beg_Apo_close_df = pd.DataFrame({'Apo Closed': rmsd_beg_1sug})
beg_L192F_df = pd.DataFrame({'L192F': rmsd_beg_L192F})
beg_L195F_df = pd.DataFrame({'L195F': rmsd_beg_L195F})
beg_F280Y_df = pd.DataFrame({'F280Y': rmsd_beg_F280Y})
beg_E276F_df = pd.DataFrame({'E276F': rmsd_beg_E276F})
beg_F196A_df = pd.DataFrame({'F196A': rmsd_beg_F196A})
beg_V287T_df = pd.DataFrame({'V287T': rmsd_beg_V287T})

df = pd.concat([beg_Apo_open_df, beg_Apo_close_df, beg_L192F_df, beg_E276F_df, beg_V287T_df, beg_F196A_df, beg_F280Y_df, beg_L195F_df])

sns.kdeplot(data = df, fill=True, alpha=0.5, common_grid = True)
plt.xlabel(r'RMSD($\AA$)', fontsize = 12)
plt.xticks(fontsize = 10)
plt.yticks(fontsize = 10)
plt.ylabel(r'Normalized Density', fontsize = 12)
plt.title(r'L1 RMSD Compared to WT Closed', fontsize = 14)
plt.savefig('mutate_RMSD_beg_all.png')
plt.close()

df = pd.concat([beg_L192F_df, beg_E276F_df, beg_V287T_df, beg_F196A_df, beg_F280Y_df, beg_L195F_df])

sns.kdeplot(data = df, fill=True, alpha=0.5, common_grid = True)
plt.xlabel(r'RMSD($\AA$)', fontsize = 12)
plt.xticks(fontsize = 10)
plt.yticks(fontsize = 10)
plt.ylabel(r'Normalized Density', fontsize = 12)
plt.title(r'L1 RMSD Compared to WT Closed', fontsize = 14)
plt.savefig('mutate_RMSD_beg_mut_all.png')
plt.close()

#Determine p-values for each of the sections of focus
file_p = open('p_values_mut.txt', 'w')
p = np.zeros((5, 7))
st, p[0,0] = stats.ttest_ind(rmsd_a3_top_apo, rmsd_a3_top_L192F, equal_var = False) #Welch's t-test
st, p[0,1] = stats.ttest_ind(rmsd_a3_top_apo, rmsd_a3_top_E276F, equal_var = False) #Welch's t-test
st, p[0,2] = stats.ttest_ind(rmsd_a3_top_apo, rmsd_a3_top_F280Y, equal_var = False) #Welch's t-test
st, p[0,3] = stats.ttest_ind(rmsd_a3_top_apo, rmsd_a3_top_L195F, equal_var = False) #Welch's t-test
st, p[0,4] = stats.ttest_ind(rmsd_a3_top_apo, rmsd_a3_top_F196A, equal_var = False) #Welch's t-test
st, p[0,5] = stats.ttest_ind(rmsd_a3_top_apo, rmsd_a3_top_V287T, equal_var = False) #Welch's t-test
st, p[0,6] = stats.ttest_ind(rmsd_a3_top_apo, rmsd_a3_top_1sug, equal_var = False) #Welch's t-test

st, p[1,0] = stats.ttest_ind(rmsd_a3_apo, rmsd_a3_L192F, equal_var = False) #Welch's t-test
st, p[1,1] = stats.ttest_ind(rmsd_a3_apo, rmsd_a3_E276F, equal_var = False) #Welch's t-test
st, p[1,2] = stats.ttest_ind(rmsd_a3_apo, rmsd_a3_F280Y, equal_var = False) #Welch's t-test
st, p[1,3] = stats.ttest_ind(rmsd_a3_apo, rmsd_a3_L195F, equal_var = False) #Welch's t-test
st, p[1,4] = stats.ttest_ind(rmsd_a3_apo, rmsd_a3_F196A, equal_var = False) #Welch's t-test
st, p[1,5] = stats.ttest_ind(rmsd_a3_apo, rmsd_a3_V287T, equal_var = False) #Welch's t-test
st, p[1,6] = stats.ttest_ind(rmsd_a3_apo, rmsd_a3_1sug, equal_var = False) #Welch's t-test

st, p[1,0] = stats.ttest_ind(rmsd_a4_apo, rmsd_a4_L192F, equal_var = False) #Welch's t-test
st, p[1,1] = stats.ttest_ind(rmsd_a4_apo, rmsd_a4_E276F, equal_var = False) #Welch's t-test
st, p[1,2] = stats.ttest_ind(rmsd_a4_apo, rmsd_a4_F280Y, equal_var = False) #Welch's t-test
st, p[1,3] = stats.ttest_ind(rmsd_a4_apo, rmsd_a4_L195F, equal_var = False) #Welch's t-test
st, p[1,4] = stats.ttest_ind(rmsd_a4_apo, rmsd_a4_F196A, equal_var = False) #Welch's t-test
st, p[1,5] = stats.ttest_ind(rmsd_a4_apo, rmsd_a4_V287T, equal_var = False) #Welch's t-test
st, p[1,6] = stats.ttest_ind(rmsd_a4_apo, rmsd_a4_1sug, equal_var = False) #Welch's t-test

st, p[2,0] = stats.ttest_ind(rmsd_a6_bot_apo, rmsd_a6_bot_L192F, equal_var = False) #Welch's t-test
st, p[2,1] = stats.ttest_ind(rmsd_a6_bot_apo, rmsd_a6_bot_E276F, equal_var = False) #Welch's t-test
st, p[2,2] = stats.ttest_ind(rmsd_a6_bot_apo, rmsd_a6_bot_F280Y, equal_var = False) #Welch's t-test
st, p[2,3] = stats.ttest_ind(rmsd_a6_bot_apo, rmsd_a6_bot_L195F, equal_var = False) #Welch's t-test
st, p[2,4] = stats.ttest_ind(rmsd_a6_bot_apo, rmsd_a6_bot_F196A, equal_var = False) #Welch's t-test
st, p[2,5] = stats.ttest_ind(rmsd_a6_bot_apo, rmsd_a6_bot_V287T, equal_var = False) #Welch's t-test
st, p[2,6] = stats.ttest_ind(rmsd_a6_bot_apo, rmsd_a6_bot_1sug, equal_var = False) #Welch's t-test

st, p[3,0] = stats.ttest_ind(rmsd_CYS_apo, rmsd_CYS_L192F, equal_var = False) #Welch's t-test
st, p[3,1] = stats.ttest_ind(rmsd_CYS_apo, rmsd_CYS_E276F, equal_var = False) #Welch's t-test
st, p[3,2] = stats.ttest_ind(rmsd_CYS_apo, rmsd_CYS_F280Y, equal_var = False) #Welch's t-test
st, p[3,3] = stats.ttest_ind(rmsd_CYS_apo, rmsd_CYS_L195F, equal_var = False) #Welch's t-test
st, p[3,4] = stats.ttest_ind(rmsd_CYS_apo, rmsd_CYS_F196A, equal_var = False) #Welch's t-test
st, p[3,5] = stats.ttest_ind(rmsd_CYS_apo, rmsd_CYS_V287T, equal_var = False) #Welch's t-test
st, p[3,6] = stats.ttest_ind(rmsd_CYS_apo, rmsd_CYS_1sug, equal_var = False) #Welch's t-test

sections_mini = ['a3_top', 'a3', 'a4', 'a6_bot']
Labels_mut = ['L192F', 'E276F', 'F280Y', 'L195F', 'F196A', 'V287T', 'Apo Closed']
file_p.write('P values of RMSD with Apo closed reference structure Relative to Apo Open \n')
for i in range(len(sections_mini)):
    file_p.write(str(sections_mini[i]) + '\n') 
    for j in range(len(Labels_mut)):
        file_p.write(Labels_mut[j] + ': ' + str(p[i,j]) + '\n')

p = np.zeros((5, 12))
st, p[0,0] = stats.ttest_ind(rmsd_a3_top_L192F, rmsd_a3_top_L192F_AD, equal_var = False) #Welch's t-test
st, p[0,1] = stats.ttest_ind(rmsd_a3_top_E276F, rmsd_a3_top_E276F_AD, equal_var = False) #Welch's t-test
st, p[0,2] = stats.ttest_ind(rmsd_a3_top_F280Y, rmsd_a3_top_F280Y_AD, equal_var = False) #Welch's t-test
st, p[0,3] = stats.ttest_ind(rmsd_a3_top_L195F, rmsd_a3_top_L195F_AD, equal_var = False) #Welch's t-test
st, p[0,4] = stats.ttest_ind(rmsd_a3_top_F196A, rmsd_a3_top_F196A_AD, equal_var = False) #Welch's t-test
st, p[0,5] = stats.ttest_ind(rmsd_a3_top_V287T, rmsd_a3_top_V287T_AD, equal_var = False) #Welch's t-test
st, p[0,6] = stats.ttest_ind(rmsd_a3_top_L192F, rmsd_a3_top_L192F_BBR, equal_var = False) #Welch's t-test
st, p[0,7] = stats.ttest_ind(rmsd_a3_top_E276F, rmsd_a3_top_E276F_BBR, equal_var = False) #Welch's t-test
st, p[0,8] = stats.ttest_ind(rmsd_a3_top_F280Y, rmsd_a3_top_F280Y_BBR, equal_var = False) #Welch's t-test
st, p[0,9] = stats.ttest_ind(rmsd_a3_top_L195F, rmsd_a3_top_L195F_BBR, equal_var = False) #Welch's t-test
st, p[0,10] = stats.ttest_ind(rmsd_a3_top_F196A, rmsd_a3_top_F196A_BBR, equal_var = False) #Welch's t-test
st, p[0,11] = stats.ttest_ind(rmsd_a3_top_V287T, rmsd_a3_top_V287T_BBR, equal_var = False) #Welch's t-test

st, p[1,0] = stats.ttest_ind(rmsd_a3_L192F, rmsd_a3_L192F_AD, equal_var = False) #Welch's t-test
st, p[1,1] = stats.ttest_ind(rmsd_a3_E276F, rmsd_a3_E276F_AD, equal_var = False) #Welch's t-test
st, p[1,2] = stats.ttest_ind(rmsd_a3_F280Y, rmsd_a3_F280Y_AD, equal_var = False) #Welch's t-test
st, p[1,3] = stats.ttest_ind(rmsd_a3_L195F, rmsd_a3_L195F_AD, equal_var = False) #Welch's t-test
st, p[1,4] = stats.ttest_ind(rmsd_a3_F196A, rmsd_a3_F196A_AD, equal_var = False) #Welch's t-test
st, p[1,5] = stats.ttest_ind(rmsd_a3_V287T, rmsd_a3_V287T_AD, equal_var = False) #Welch's t-test
st, p[1,6] = stats.ttest_ind(rmsd_a3_L192F, rmsd_a3_L192F_BBR, equal_var = False) #Welch's t-test
st, p[1,7] = stats.ttest_ind(rmsd_a3_E276F, rmsd_a3_E276F_BBR, equal_var = False) #Welch's t-test
st, p[1,8] = stats.ttest_ind(rmsd_a3_F280Y, rmsd_a3_F280Y_BBR, equal_var = False) #Welch's t-test
st, p[1,9] = stats.ttest_ind(rmsd_a3_L195F, rmsd_a3_L195F_BBR, equal_var = False) #Welch's t-test
st, p[1,10] = stats.ttest_ind(rmsd_a3_F196A, rmsd_a3_F196A_BBR, equal_var = False) #Welch's t-test
st, p[1,11] = stats.ttest_ind(rmsd_a3_V287T, rmsd_a3_V287T_BBR, equal_var = False) #Welch's t-test

st, p[2,0] = stats.ttest_ind(rmsd_a4_L192F, rmsd_a4_L192F_AD, equal_var = False) #Welch's t-test
st, p[2,1] = stats.ttest_ind(rmsd_a4_E276F, rmsd_a4_E276F_AD, equal_var = False) #Welch's t-test
st, p[2,2] = stats.ttest_ind(rmsd_a4_F280Y, rmsd_a4_F280Y_AD, equal_var = False) #Welch's t-test
st, p[2,3] = stats.ttest_ind(rmsd_a4_L195F, rmsd_a4_L195F_AD, equal_var = False) #Welch's t-test
st, p[2,4] = stats.ttest_ind(rmsd_a4_F196A, rmsd_a4_F196A_AD, equal_var = False) #Welch's t-test
st, p[2,5] = stats.ttest_ind(rmsd_a4_V287T, rmsd_a4_V287T_AD, equal_var = False) #Welch's t-test
st, p[2,6] = stats.ttest_ind(rmsd_a4_L192F, rmsd_a4_L192F_BBR, equal_var = False) #Welch's t-test
st, p[2,7] = stats.ttest_ind(rmsd_a4_E276F, rmsd_a4_E276F_BBR, equal_var = False) #Welch's t-test
st, p[2,8] = stats.ttest_ind(rmsd_a4_F280Y, rmsd_a4_F280Y_BBR, equal_var = False) #Welch's t-test
st, p[2,9] = stats.ttest_ind(rmsd_a4_L195F, rmsd_a4_L195F_BBR, equal_var = False) #Welch's t-test
st, p[2,10] = stats.ttest_ind(rmsd_a4_F196A, rmsd_a4_F196A_BBR, equal_var = False) #Welch's t-test
st, p[2,11] = stats.ttest_ind(rmsd_a4_V287T, rmsd_a4_V287T_BBR, equal_var = False) #Welch's t-test

sections_mini = ['a3_top', 'a3', 'a4']
Labels_mut = ['L192F', 'E276F', 'F280Y', 'L195F', 'F196A', 'V287T']
file_p.write('P values of RMSD with Apo closed reference structure Relative to Apo Mut \n')
for i in range(len(sections_mini)):
    file_p.write(str(sections_mini[i]) + '\n') 
    for j in range(len(Labels_mut)):
        n = j+6
        file_p.write(Labels_mut[j] + ' AD: ' + str(p[i,j]) + '\n')
        file_p.write(Labels_mut[j] + ' BBR: ' + str(p[i,n]) + '\n')

