#Import Necessary Pacjages
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

def adjacent_values(vals, q1, q3):
    upper_adjacent_value = q3 + (q3 - q1) * 1.5
    upper_adjacent_value = np.clip(upper_adjacent_value, q3, vals[-1])

    lower_adjacent_value = q1 - (q3 - q1) * 1.5
    lower_adjacent_value = np.clip(lower_adjacent_value, vals[0], q1)
    return lower_adjacent_value, upper_adjacent_value

def set_axis_style(ax, labels):
    ax.xaxis.set_tick_params(direction='out')
    ax.xaxis.set_ticks_position('bottom')
    ax.set_xticks(np.arange(1, len(labels) + 1))
    ax.set_xticklabels(labels)
    ax.set_xlim(0.25, len(labels) + 0.75)
    ax.set_xlabel('Sample name')

def load_txt(folder_path, file_name):
    x,y = [],[]
    #Load data
    with open('../../' + folder_path + '/rmsf_full.xvg', 'r') as f:
        for _ in range(18):
            next(f)
        for line in f:
            cols = line.split()
            x.append(float(cols[0]))
            y.append(float(cols[1])*10)

    return x, y

#Input Files
t_a7, r_a7 = load_txt('rebuild_a7/analysis', 'a7')
t_a7_AD, r_a7_AD = load_txt('AD_a7/analysis', 'AD_a7')
t_dis7, r_dis7 = load_txt('rebuild_a7_high/config7/analysis', 'Apo_dis7')
t_dis9, r_dis9 = load_txt('rebuild_a7_high/config9/analysis', 'Apo_dis9')
t_dis11, r_dis11 = load_txt('rebuild_a7_high/config11/analysis', 'Apo_dis11')
t_dis_alt, r_dis_alt = load_txt('Apo_dis/analysis', 'Apo_dis')
t_dis7_AD, r_dis7_AD = load_txt('AD_dis/analysis/config7', 'AD_dis7')
t_dis9_AD, r_dis9_AD = load_txt('AD_dis/analysis/config9', 'AD_dis9')
t_dis11_AD, r_dis11_AD = load_txt('AD_dis/analysis/config11', 'AD_dis11')
t_dis_alt_AD, r_dis_alt_AD = load_txt('AD_dis/analysis/config_alt', 'AD_dis_alt')
t_1sug, r_1sug = load_txt('Apo_1SUG/analysis/1sug', '1sug')
t_1sug_AD, r_1sug_AD = load_txt('1sug_AD/analysis', '1sug_AD')
t_1sug_dis7, r_1sug_dis7 = load_txt('1sug_dis/analysis/config7', '1sug_dis7')
t_1sug_dis9, r_1sug_dis9 = load_txt('1sug_dis/analysis/config9', '1sug_dis9')
t_1sug_dis11, r_1sug_dis11 = load_txt('1sug_dis/analysis/config11', '1sug_dis11')
t_1sug_dis7_AD, r_1sug_dis7_AD = load_txt('1sug_dis_AD/analysis/config7', '1sug_dis7_AD')
t_1sug_dis9_AD, r_1sug_dis9_AD = load_txt('1sug_dis_AD/analysis/config9', '1sug_dis9_AD')
t_1sug_dis11_AD, r_1sug_dis11_AD = load_txt('1sug_dis_AD/analysis/config11', '1sug_dis11_AD')
t_1sug_dis_alt_AD, r_1sug_dis_alt_AD = load_txt('1sug_dis_AD/analysis/config_alt', '1sug_dis_alt_AD')
t_1sug_dis_alt2_AD, r_1sug_dis_alt2_AD = load_txt('1sug_dis_AD/analysis/config_alt2', '1sug_dis_alt2')
t_BBR_a7, r_BBR_a7 = load_txt('BBR_a7/analysis', 'BBR_a7')
t_BBR_1sug, r_BBR_1sug = load_txt('BBR_1sug/analysis', 'BBR_1sug')
t_BBR_dis7, r_BBR_dis7 = load_txt('BBR_dis/analysis/config7', 'BBR_dis7')
t_BBR_dis9, r_BBR_dis9 = load_txt('BBR_dis/analysis/config9', 'BBR_dis9')
t_BBR_dis11, r_BBR_dis11 = load_txt('BBR_dis/analysis/config11', 'BBR_dis11')
t_BBR_1sug_dis7, r_BBR_1sug_dis7 = load_txt('BBR_1sug_dis/analysis/config7', 'BBR_1sug_dis7')
t_BBR_1sug_dis11, r_BBR_1sug_dis11 = load_txt('BBR_1sug_dis/analysis/config11', 'BBR_1sug_dis11')
t_AD_BBR, r_AD_BBR = load_txt('AD_BBR/analysis', 'AD_BBR')

#WPD RMSF
res_w = np.array([177, 178, 179, 180, 181, 182, 183, 184, 185]) #Residues in WPD loop
w_Apo_open, w_Apo_close, w_AD_open, w_AD_close, w_BBR_open, w_BBR_close = [],[],[],[],[],[] #Empty array for RMSF of WPD residues
for i in res_w:
    j = i-1
    w_Apo_open.append(np.mean([r_a7[j], r_dis7[j], r_1sug_dis7[j], r_dis9[j], r_1sug_dis9[j]]))
    w_Apo_close.append(np.mean([r_1sug[j], r_dis11[j], r_1sug_dis11[j]]))
    w_AD_open.append(np.mean([r_1sug_dis_alt_AD[j], r_1sug_dis_alt2_AD[j], r_dis_alt_AD[j]]))
    w_BBR_open.append(np.mean([r_BBR_a7[j], r_BBR_dis9[j]]))

#Plot WPD RMSF
fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.set_title("Comparison of WPD RMSF")    
ax1.set_xlabel('Residue ID')
ax1.set_ylabel(r'RMSF($\AA$)')
ax1.plot(res_w,w_Apo_open, label='Apo Open')
ax1.plot(res_w,w_Apo_close, label='Apo Close')
ax1.plot(res_w,w_AD_open, label='AD Open')
leg = ax1.legend(loc='upper right')
fig.savefig('WPD_RMSF_compare.png') 

#Plot WPD RMSF
fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.set_title("Comparison of WPD RMSF")    
ax1.set_xlabel('Residue ID')
ax1.set_ylabel(r'RMSF($\AA$)')
ax1.plot(res_w,w_Apo_open, label='Apo Open')
ax1.plot(res_w,w_Apo_close, label='Apo Close')
ax1.plot(res_w,w_AD_open, label='AD Open')
ax1.plot(res_w,w_BBR_open, label='BBR Open')
leg = ax1.legend(loc='upper right')
fig.savefig('WPD_RMSF_compare_BBR.png') 

#Maximum RMSF WPD
max_wpd = []
index_Apo_open = w_Apo_open.index(max(w_Apo_open))
index_Apo_close = w_Apo_close.index(max(w_Apo_close))
index_AD_open = w_AD_open.index(max(w_AD_open))
index_BBR_open = w_BBR_open.index(max(w_BBR_open))

max_wpd.append((w_Apo_open[index_Apo_open-1] + w_Apo_open[index_Apo_open] + w_Apo_open[index_Apo_open+1])/3)
max_wpd.append((w_Apo_close[index_Apo_close-1] + w_Apo_close[index_Apo_close])/2)
max_wpd.append((w_AD_open[index_AD_open-1] + w_AD_open[index_AD_open] + w_AD_open[index_AD_open+1])/3)
max_wpd.append((w_BBR_open[index_BBR_open-1] + w_BBR_open[index_BBR_open] + w_BBR_open[index_BBR_open+1])/3)

name = ['Apo Open', 'Apo Closed', 'AD Open', 'BBR Open']
num = np.array([1, 2, 3, 4])
fig4 = plt.figure()
ax1 = fig4.add_subplot(111)
ax1.set_title("Comparison of WPD RMSF Peaj")
ax1.set_ylabel(r'RMSF($\AA$)')
ax1.bar(num, max_wpd)
plt.xticks(num, name)
leg = ax1.legend()
fig4.savefig('WPD_peaj_compare.png')

#Substrate Binding loop
res_rl = np.array([112, 113, 114, 115, 116, 117])
rl_Apo_open, rl_Apo_close, rl_AD_open, rl_AD_close, rl_BBR_open, rl_BBR_close = [],[],[],[],[],[] #Empty array for RMSF of WPD residues
for i in res_rl:
    j = i-1
    rl_Apo_open.append(np.mean([r_a7[j], r_dis7[j], r_1sug_dis7[j], r_dis9[j], r_1sug_dis9[j]]))
    rl_Apo_close.append(np.mean([r_1sug[j], r_dis11[j], r_1sug_dis11[j]]))
    rl_AD_open.append(np.mean([r_1sug_dis_alt_AD[j], r_1sug_dis_alt2_AD[j], r_dis_alt_AD[j]]))
    rl_BBR_open.append(np.mean([r_BBR_a7[j], r_BBR_dis9[j]]))
res_rl = np.array([113, 114, 115, 116, 117, 118])

#Plot substrate binding loop RMSF
figr = plt.figure()
ax1 = figr.add_subplot(111)
ax1.set_title("Comparison of Substrate Binding Loop RMSF")    
ax1.set_xlabel('Residue ID')
ax1.set_ylabel(r'RMSF($\AA$)')
ax1.plot(res_rl,rl_Apo_open, label='Apo Open')
ax1.plot(res_rl,rl_Apo_close, label='Apo Closed')
ax1.plot(res_rl,rl_AD_open, label='AD')
ax1.plot(res_rl,rl_BBR_open, label='BBR')
leg = ax1.legend(loc='upper right')
figr.savefig('R_RMSF_compare.png') 


#a3 RMAF
res_a3 = np.array([186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199, 200]) #Residues in WPD loop
a3_Apo_open, a3_Apo_close, a3_AD_open, a3_AD_close, a3_BBR_open, a3_BBR_close = [],[],[],[],[],[] #Empty array for RMSF of WPD residues
for i in res_a3:
    j = i-1
    a3_Apo_open.append(np.mean([r_a7[j], r_dis7[j], r_1sug_dis7[j], r_dis9[j], r_1sug_dis9[j]]))
    a3_Apo_close.append(np.mean([r_1sug[j], r_dis11[j], r_1sug_dis11[j]]))
    a3_AD_open.append(np.mean([r_1sug_dis_alt_AD[j], r_1sug_dis_alt2_AD[j], r_dis_alt_AD[j]]))
    a3_BBR_open.append(np.mean([r_BBR_a7[j], r_BBR_dis9[j]]))

#Plot a3 RMSF
fig2 = plt.figure()
ax1 = fig2.add_subplot(111)
ax1.set_title("Comparison of a3 helix RMSF")    
ax1.set_xlabel('Residue ID')
ax1.set_ylabel(r'RMSF($\AA$)')
ax1.plot(res_a3,a3_Apo_open, label='Apo Open')
ax1.plot(res_a3,a3_Apo_close, label='Apo Closed')
ax1.plot(res_a3,a3_AD_open, label='AD Open')
ax1.plot(res_a3,a3_BBR_open, label='BBR Open')
leg = ax1.legend(loc='upper right')
fig2.savefig('a3_RMSF_compare.png') 

#a6 RMSF
res_a6 = np.array([264, 265, 266, 267, 268, 269, 270, 271, 272, 273, 274, 275, 276, 277, 278, 279, 280, 281]) #Residues in WPD loop
a6_Apo_open, a6_Apo_close, a6_AD_open, a6_AD_close, a6_BBR_open, a6_BBR_close = [],[],[],[],[],[] #Empty array for RMSF of WPD residues
for i in res_a6:
    j = i-1
    a6_Apo_open.append(np.mean([r_a7[j], r_dis7[j], r_1sug_dis7[j], r_dis9[j], r_1sug_dis9[j]]))
    a6_Apo_close.append(np.mean([r_1sug[j], r_dis11[j], r_1sug_dis11[j]]))
    a6_AD_open.append(np.mean([r_1sug_dis_alt_AD[j], r_1sug_dis_alt2_AD[j], r_dis_alt_AD[j]]))
    a6_BBR_open.append(np.mean([r_BBR_a7[j], r_BBR_dis9[j]]))

#Plot a6 RMSF
fig3 = plt.figure()
ax1 = fig3.add_subplot(111)
ax1.set_title("Comparison of a6 helix RMSF")    
ax1.set_xlabel('Residue ID')
ax1.set_ylabel(r'RMSF($\AA$)')
ax1.plot(res_a6,a6_Apo_open, label='Apo Open')
ax1.plot(res_a6,a6_Apo_close, label='Apo Closed')
ax1.plot(res_a6,a6_AD_open, label='AD Open')
ax1.plot(res_a6,a6_BBR_open, label='BBR Open')
leg = ax1.legend(loc='upper right')
fig3.savefig('a6_RMSF_compare.png') 

#a7 RMSF
res_a7 = np.array([287, 288, 289, 290, 291, 292, 293, 294, 295, 296, 297]) #Residues in WPD loop
a7_Apo_open, a7_Apo_close, a7_AD_open, a7_AD_close, a7_AD_unb, a7_BBR_open, a7_BBR_close, a7_BBR_unb = [],[],[],[],[],[],[],[] #Empty array for RMSF of WPD residues
a7_Apo_open_sem, a7_Apo_close_sem, a7_AD_open_sem, a7_AD_close_sem, a7_AD_unb_sem, a7_BBR_open_sem, a7_BBR_close_sem, a7_BBR_unb_sem = [],[],[],[],[],[],[],[] #Empty array for Error for RMSF of WPD residues
a7_unb1, a7_unb2, a7_unb3, a7_unb4, a7_unb5, a7_unb6, a7_unb7 = [],[],[],[],[],[],[]
a7_dis7, a7_dis9, a7_dis11, a7_dis_unb, a7_ord = [],[],[],[],[]
a7_lig_open, a7_lig_open_sem, a7_dis11_sem, a7_dis_unb_sem = [],[],[],[]
a7_lig_unb, a7_lig_unb_sem = [],[]
a7_lig_bound, a7_lig_bound_sem = [],[]
a7_AD_BBR = []
for i in res_a7:
    j = i-1
    a7_Apo_open.append(np.mean([r_a7[j], r_dis7[j], r_1sug_dis7[j], r_dis9[j], r_1sug_dis9[j]]))
    a7_Apo_close.append(np.mean([r_1sug[j], r_dis11[j], r_1sug_dis11[j]]))
    a7_AD_open.append(np.mean([r_1sug_dis_alt_AD[j], r_1sug_dis_alt2_AD[j], r_dis_alt_AD[j]]))
    a7_BBR_open.append(np.mean([r_BBR_a7[j], r_BBR_dis9[j]]))
    a7_lig_unb.append(np.mean([r_dis7[j], r_1sug_dis7[j], r_dis9[j], r_1sug_dis9[j], r_1sug_dis11[j], r_dis11[j]]))
    a7_lig_bound.append(r_dis_alt[j])
    a7_AD_BBR.append(r_AD_BBR[j])

    a7_Apo_open_sem.append(stats.sem([r_a7[j], r_dis7[j], r_1sug_dis7[j], r_dis9[j], r_1sug_dis9[j]]))
    a7_Apo_close_sem.append(stats.sem([r_1sug[j], r_dis11[j], r_1sug_dis11[j]]))
    a7_AD_open_sem.append(stats.sem([r_1sug_dis_alt_AD[j], r_1sug_dis_alt2_AD[j], r_dis_alt_AD[j]]))
    a7_BBR_open_sem.append(stats.sem([r_BBR_a7[j], r_BBR_dis9[j]]))
    a7_lig_unb_sem.append(stats.sem([r_dis7[j], r_1sug_dis7[j], r_dis9[j], r_1sug_dis9[j], r_1sug_dis11[j], r_dis11[j]]))
    a7_lig_bound_sem.append(0)

#Plot a7 RMSF
fig4 = plt.figure()
ax1 = fig4.add_subplot(111)
ax1.set_title(r'Comparison of $\alpha$-7 helix RMSF')    
ax1.set_xlabel('Residue ID')
ax1.set_ylabel(r'RMSF($\AA$)')
ax1.errorbar(res_a7, a7_Apo_open, yerr = a7_Apo_open_sem, label='Apo Open')
ax1.errorbar(res_a7, a7_Apo_close, yerr = a7_Apo_close_sem, label='Apo Closed')
ax1.errorbar(res_a7, a7_AD_open, yerr = a7_AD_open_sem, label='AD')
ax1.errorbar(res_a7, a7_BBR_open, yerr = a7_BBR_open_sem, label='BBR')
ax1.legend(loc='upper right')
fig4.savefig('a7_RMSF_compare_err.png') 

#Plot a7 RMSF
fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.set_title(r'Comparison of $\alpha$-7 helix RMSF')    
ax1.set_xlabel('Residue ID')
ax1.set_ylabel(r'RMSF($\AA$)')
ax1.plot(res_a7, a7_Apo_open, label='Apo Open', color = 'blue')
ax1.plot(res_a7, a7_Apo_close, label='Apo Closed', color = 'orange')
ax1.legend(loc='upper right')
fig.savefig('a7_RMSF_Apo.png') 
plt.close(fig)

fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.set_title(r"Comparison of $\alpha$-7 helix RMSF")    
ax1.set_xlabel('Residue ID')
ax1.set_ylabel(r'RMSF($\AA$)')
ax1.plot(res_a7, a7_Apo_open, label='Apo Open', color = 'blue')
ax1.plot(res_a7, a7_Apo_close, label='Apo Closed', color = 'orange')
ax1.plot(res_a7, a7_AD_open, label='AD', color = 'purple')
ax1.legend(loc='upper right')
fig.savefig('a7_RMSF_AD.png') 
plt.close(fig)

fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.set_title(r"Comparison of $\alpha$-7 helix RMSF")    
ax1.set_xlabel('Residue ID')
ax1.set_ylabel(r'RMSF($\AA$)')
ax1.plot(res_a7, a7_Apo_open, label='Apo Open', color = 'blue')
ax1.plot(res_a7, a7_AD_open, label='AD', color = 'blue', linestyle = 'dashed')
ax1.plot(res_a7, a7_BBR_open, label='BBR', color = 'blue', linestyle = 'dotted')
ax1.plot(res_a7, a7_Apo_close, label='Apo Closed', color = 'red')
ax1.legend(loc='best')
fig.savefig('a7_RMSF_compare.png') 

#convert to numpy array
a7_AD_open = np.array(a7_AD_open)
a7_BBR_open = np.array(a7_BBR_open)

fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.set_title(r'Comparison of $\alpha$-7 helix RMSF', fontsize=15)
ax1.set_xlabel('Residue ID', fontsize=13)
ax1.set_ylabel(r'RMSF($\AA$)', fontsize=13)
ax1.plot(res_a7, a7_AD_open, label='AD', color = 'blue', linestyle = 'dashed')
ax1.fill_between(res_a7, a7_AD_open-a7_AD_open_sem, a7_AD_open+a7_AD_open_sem, alpha=0.1, facecolor = 'blue', edgecolor = 'blue')
ax1.plot(res_a7, a7_BBR_open, label='BBR', color = 'purple', linestyle = 'dotted')
ax1.fill_between(res_a7, a7_BBR_open-a7_BBR_open_sem, a7_BBR_open+a7_BBR_open_sem, alpha=0.1, facecolor = 'purple', edgecolor = 'purple')
ax1.plot(res_a7, a7_AD_BBR, label='AD BBR Combo', color = 'green', linestyle = 'dashed')
ax1.axvspan(290, 295, alpha = 0.2, color = 'gray')
plt.xticks(fontsize=11)
plt.yticks(fontsize=11)
ax1.legend(loc='best', fontsize=12)
fig.savefig('a7_RMSF_combo_compare.png') 

a7_Apo_open = np.array(a7_Apo_open)
a7_Apo_open_sem = np.array(a7_Apo_open_sem)
a7_AD_open = np.array(a7_AD_open)
a7_AD_open_sem = np.array(a7_AD_open_sem)
a7_BBR_open = np.array(a7_BBR_open)
a7_BBR_open_sem = np.array(a7_BBR_open_sem)
a7_Apo_close = np.array(a7_Apo_close)
a7_Apo_close_sem = np.array(a7_Apo_close_sem)

fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.set_title(r'Comparison of $\alpha$-7 helix RMSF', fontsize=15)
ax1.set_xlabel('Residue ID', fontsize=13)
ax1.set_ylabel(r'RMSF($\AA$)', fontsize=13)
ax1.plot(res_a7, a7_Apo_open, label='Apo Open', color = 'blue')
ax1.fill_between(res_a7, a7_Apo_open-a7_Apo_open_sem, a7_Apo_open+a7_Apo_open_sem, alpha=0.2, facecolor = 'blue', edgecolor = 'blue')
ax1.plot(res_a7, a7_AD_open, label='AD', color = 'blue', linestyle = 'dashed')
ax1.fill_between(res_a7, a7_AD_open-a7_AD_open_sem, a7_AD_open+a7_AD_open_sem, alpha=0.15, facecolor = 'blue', edgecolor = 'blue')
ax1.plot(res_a7, a7_BBR_open, label='BBR', color = 'blue', linestyle = 'dotted')
ax1.fill_between(res_a7, a7_BBR_open-a7_BBR_open_sem, a7_BBR_open+a7_BBR_open_sem, alpha=0.1, facecolor = 'blue', edgecolor = 'blue')
ax1.plot(res_a7, a7_Apo_close, label='Apo Closed', color = 'red')
ax1.fill_between(res_a7, a7_Apo_close-a7_Apo_close_sem, a7_Apo_close+a7_Apo_close_sem, alpha=0.2, facecolor = 'red', edgecolor = 'red')
plt.xticks(fontsize=11)
plt.yticks(fontsize=11)
ax1.legend(loc='best', fontsize=12)
fig.savefig('a7_RMSF_compare_err.png') 

#Plot a7 Unbound RMSF
a7_lig_bound = np.array(a7_lig_bound)
a7_lig_bound_sem = np.array(a7_lig_bound_sem)
a7_lig_unb = np.array(a7_lig_unb)
a7_lig_unb_sem = np.array(a7_lig_unb_sem)

fig4 = plt.figure()
ax1 = fig4.add_subplot(111)
ax1.set_title(r"Comparison of $\alpha$-7 helix Flexibility", fontsize = 16)
ax1.set_xlabel('Residue ID', fontsize = 12)
ax1.set_ylabel(r'RMSF($\AA$)', fontsize = 12)
ax1.plot(res_a7, a7_lig_bound, label='Ligand Binds', color = 'blue')
ax1.fill_between(res_a7, a7_lig_bound-a7_lig_bound_sem, a7_lig_bound+a7_lig_bound_sem, alpha=0.5, facecolor = 'blue', edgecolor = 'blue')
ax1.plot(res_a7, a7_lig_unb, label='Ligand Does Not Bind', color = 'red')
ax1.fill_between(res_a7, a7_lig_unb-a7_lig_unb_sem, a7_lig_unb+a7_lig_unb_sem, alpha=0.5, facecolor = 'red', edgecolor = 'red')
ax1.legend(loc='upper left', fontsize = 12)
fig4.savefig('a7_RMSF_lig_unbound_err.png') 
plt.close(fig4)

