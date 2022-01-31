#Import Necessary Packages
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

#Make open arrays for residues and rmsf values
t_a7, t_a7_AD, t_Apo, t_AD, t_1sug, t_1sug_AD, t_1sug_na7, t_1sug_na7_AD = [],[],[],[],[],[],[],[] #time
t_dis7, t_dis9, t_dis11, t_dis7_AD, t_dis9_AD, t_dis11_AD = [],[],[],[],[],[]
t_1sug_dis7, t_1sug_dis9, t_1sug_dis11, t_1sug_dis7_AD, t_1sug_dis9_AD, t_1sug_dis11_AD, t_1sug_dis11_2AD, t_1sug_dis_alt_AD, t_1sug_dis_alt2_AD = [],[],[],[],[],[],[],[],[]
t_BBR_a7, t_BBR_1sug, t_BBR_dis7, t_BBR_dis9, t_BBR_dis11, t_BBR_1sug_dis7, t_BBR_1sug_dis11 = [],[],[],[],[],[],[]
r_a7, r_a7_AD, r_Apo, r_AD, r_1sug, r_1sug_AD, r_1sug_na7, r_1sug_na7_AD = [],[],[],[],[],[],[],[] #rmsf for trajectory residues
r_dis7, r_dis9, r_dis11, r_dis7_AD, r_dis9_AD, r_dis11_AD = [],[],[],[],[],[]
r_1sug_dis7, r_1sug_dis9, r_1sug_dis11, r_1sug_dis7_AD, r_1sug_dis9_AD, r_1sug_dis11_AD, r_1sug_dis11_2AD, r_1sug_dis_alt_AD, r_1sug_dis_alt2_AD = [],[],[],[],[],[],[],[],[]
r_BBR_a7, r_BBR_1sug, r_BBR_dis7, r_BBR_dis9, r_BBR_dis11, r_BBR_1sug_dis7, r_BBR_1sug_dis11 = [],[],[],[],[],[],[]

#Input Files
with open("../../rebuild_a7/rmsf_a7_equil.xvg") as f:
    for _ in range(17):
        next(f)
    for line in f:
        cols = line.split()
        if len(cols) == 2:
            t_a7.append(float(cols[0]))
            r_a7.append(float(cols[1]))
with open("../../AD_rebuild_a7/rmsf_a7_AD_equil.xvg") as f:
    for _ in range(17):
        next(f)
    for line in f:
        cols = line.split()
        if len(cols) == 2:
            t_a7_AD.append(float(cols[0]))
            r_a7_AD.append(float(cols[1]))
with open("../../Apo_dis/config7/rmsf_Apo_dis7_equil.xvg") as f:
    for _ in range(17):
        next(f)
    for line in f:
        cols = line.split()
        if len(cols) == 2:
            t_dis7.append(float(cols[0]))
            r_dis7.append(float(cols[1]))
with open("../../AD_dis/config7/rmsf_AD_dis7_equil.xvg") as f:
    for _ in range(17):
        next(f)
    for line in f:
        cols = line.split()
        if len(cols) == 2:
            t_dis7_AD.append(float(cols[0]))
            r_dis7_AD.append(float(cols[1]))
with open("../../Apo_dis/config9/rmsf_Apo_dis9_equil.xvg") as f:
    for _ in range(17):
        next(f)
    for line in f:
        cols = line.split()
        if len(cols) == 2:
            t_dis9.append(float(cols[0]))
            r_dis9.append(float(cols[1]))
with open("../../AD_dis/config9/rmsf_AD_dis9_equil.xvg") as f:
    for _ in range(17):
        next(f)
    for line in f:
        cols = line.split()
        if len(cols) == 2:
            t_dis9_AD.append(float(cols[0]))
            r_dis9_AD.append(float(cols[1]))
with open("../../Apo_dis/config11/rmsf_Apo_dis11_equil.xvg") as f:
    for _ in range(17):
        next(f)
    for line in f:
        cols = line.split()
        if len(cols) == 2:
            t_dis11.append(float(cols[0]))
            r_dis11.append(float(cols[1]))
with open("../../AD_dis/config11/rmsf_AD_dis11_equil.xvg") as f:
    for _ in range(17):
        next(f)
    for line in f:
        cols = line.split()
        if len(cols) == 2:
            t_dis11_AD.append(float(cols[0]))
            r_dis11_AD.append(float(cols[1]))
with open("../../1sug/rmsf_1sug_equil.xvg") as f:
    for _ in range(17):
        next(f)
    for line in f:
        cols = line.split()
        if len(cols) == 2:
            t_1sug.append(float(cols[0]))
            r_1sug.append(float(cols[1]))
with open("../../1sug_AD/rmsf_1sug_AD_equil.xvg") as f:
    for _ in range(17):
        next(f)
    for line in f:
        cols = line.split()
        if len(cols) == 2:
            t_1sug_AD.append(float(cols[0]))
            r_1sug_AD.append(float(cols[1]))
with open("../../1sug_dis/config7/rmsf_1sug_dis7_equil.xvg") as f:
    for _ in range(17):
        next(f)
    for line in f:
        cols = line.split()
        if len(cols) == 2:
            t_1sug_dis7.append(float(cols[0]))
            r_1sug_dis7.append(float(cols[1]))
with open("../../1sug_dis_AD/config7/rmsf_1sug_dis7_AD_equil.xvg") as f:
    for _ in range(17):
        next(f)
    for line in f:
        cols = line.split()
        if len(cols) == 2:
            t_1sug_dis7_AD.append(float(cols[0]))
            r_1sug_dis7_AD.append(float(cols[1]))
with open("../../1sug_dis/config9/rmsf_1sug_dis9_equil.xvg") as f:
    for _ in range(17):
        next(f)
    for line in f:
        cols = line.split()
        if len(cols) == 2:
            t_1sug_dis9.append(float(cols[0]))
            r_1sug_dis9.append(float(cols[1]))
with open("../../1sug_dis_AD/config9/rmsf_1sug_dis9_AD_equil.xvg") as f:
    for _ in range(17):
        next(f)
    for line in f:
        cols = line.split()
        if len(cols) == 2:
            t_1sug_dis9_AD.append(float(cols[0]))
            r_1sug_dis9_AD.append(float(cols[1]))
with open("../../1sug_dis/config11/rmsf_1sug_dis11_equil.xvg") as f:
    for _ in range(17):
        next(f)
    for line in f:
        cols = line.split()
        if len(cols) == 2:
            t_1sug_dis11.append(float(cols[0]))
            r_1sug_dis11.append(float(cols[1]))
with open("../../1sug_dis_AD/config11/rmsf_1sug_dis11_AD_equil.xvg") as f:
    for _ in range(17):
        next(f)
    for line in f:
        cols = line.split()
        if len(cols) == 2:
            t_1sug_dis11_AD.append(float(cols[0]))
            r_1sug_dis11_AD.append(float(cols[1]))
with open("../../1sug_dis_AD/config11_2/rmsf_1sug_dis11_2_AD_equil.xvg") as f:
    for _ in range(17):
        next(f)
    for line in f:
        cols = line.split()
        if len(cols) == 2:
            t_1sug_dis11_2AD.append(float(cols[0]))
            r_1sug_dis11_2AD.append(float(cols[1]))
with open("../../1sug_AD_dis_alt/run_1/rmsf_1sug_dis_alt_equil.xvg") as f:
    for _ in range(17):
        next(f)
    for line in f:
        cols = line.split()
        if len(cols) == 2:
            t_1sug_dis_alt_AD.append(float(cols[0]))
            r_1sug_dis_alt_AD.append(float(cols[1]))
with open("../../1sug_AD_dis_alt/run_2/rmsf_1sug_dis_alt2_equil.xvg") as f:
    for _ in range(17):
        next(f)
    for line in f:
        cols = line.split()
        if len(cols) == 2:
            t_1sug_dis_alt2_AD.append(float(cols[0]))
            r_1sug_dis_alt2_AD.append(float(cols[1]))
with open("../../BBR_a7/rmsf_BBR_a7_equil.xvg") as f:
    for _ in range(17):
        next(f)
    for line in f:
        cols = line.split()
        if len(cols) == 2:
            t_BBR_a7.append(float(cols[0]))
            r_BBR_a7.append(float(cols[1]))
with open("../../1sug_BBR/rmsf_BBR_1sug_equil.xvg") as f:
    for _ in range(17):
        next(f)
    for line in f:
        cols = line.split()
        if len(cols) == 2:
            t_BBR_1sug.append(float(cols[0]))
            r_BBR_1sug.append(float(cols[1]))
with open("../../BBR_dis/config7/rmsf_BBR_dis7_equil.xvg") as f:
    for _ in range(17):
        next(f)
    for line in f:
        cols = line.split()
        if len(cols) == 2:
            t_BBR_dis7.append(float(cols[0]))
            r_BBR_dis7.append(float(cols[1]))
with open("../../BBR_dis/config9/rmsf_BBR_dis9_equil.xvg") as f:
    for _ in range(17):
        next(f)
    for line in f:
        cols = line.split()
        if len(cols) == 2:
            t_BBR_dis9.append(float(cols[0]))
            r_BBR_dis9.append(float(cols[1]))
with open("../../BBR_dis/config11/rmsf_BBR_dis11_equil.xvg") as f:
    for _ in range(17):
        next(f)
    for line in f:
        cols = line.split()
        if len(cols) == 2:
            t_BBR_dis11.append(float(cols[0]))
            r_BBR_dis11.append(float(cols[1]))
with open("../../1sug_dis_BBR/config7/rmsf_BBR_1sug_dis7_equil.xvg") as f:
    for _ in range(17):
        next(f)
    for line in f:
        cols = line.split()
        if len(cols) == 2:
            t_BBR_1sug_dis7.append(float(cols[0]))
            r_BBR_1sug_dis7.append(float(cols[1]))
with open("../../1sug_dis_BBR/config11/rmsf_BBR_1sug_dis11_equil.xvg") as f:
    for _ in range(17):
        next(f)
    for line in f:
        cols = line.split()
        if len(cols) == 2:
            t_BBR_1sug_dis11.append(float(cols[0]))
            r_BBR_1sug_dis11.append(float(cols[1]))

#WPD RMSF
res_w = np.array([177, 178, 179, 180, 181, 182, 183, 184, 185]) #Residues in WPD loop
w_Apo_open, w_Apo_close, w_AD_open, w_AD_close, w_BBR_open, w_BBR_close = [],[],[],[],[],[] #Empty array for RMSF of WPD residues
for i in res_w:
    j = i-1
    k = i-2
    w_Apo_open.append(np.mean([r_a7[j], r_dis7[j], r_1sug_dis7[k], r_dis9[j], r_1sug_dis9[k]]))
    w_Apo_close.append(np.mean([r_1sug[k], r_dis11[j], r_1sug_dis11[k]]))
    w_AD_open.append(np.mean([r_1sug_dis11_AD[k], r_1sug_dis_alt_AD[k], r_1sug_dis_alt2_AD[k]]))
    w_AD_close.append(np.mean([r_dis11_AD[j]]))
    w_BBR_open.append(np.mean([r_BBR_a7[j], r_BBR_dis9[j]]))
    w_BBR_close.append(np.mean([r_BBR_1sug[k], r_BBR_1sug_dis7[k]]))

#Plot WPD RMSF
fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.set_title("Comparison of WPD RMSF")    
ax1.set_xlabel('Residues')
ax1.set_ylabel('RMSF(nm)')
ax1.plot(res_w,w_Apo_open, label='Apo Open')
ax1.plot(res_w,w_Apo_close, label='Apo Close')
ax1.plot(res_w,w_AD_open, label='AD Open')
ax1.plot(res_w,w_AD_close, label='AD Close')
leg = ax1.legend(loc='upper right')
fig.savefig('WPD_RMSF_compare.png') 

#Plot WPD RMSF
fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.set_title("Comparison of WPD RMSF")    
ax1.set_xlabel('Residues')
ax1.set_ylabel('RMSF(nm)')
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
index_AD_close = w_AD_close.index(max(w_AD_close))
index_BBR_open = w_BBR_open.index(max(w_BBR_open))
index_BBR_close = w_BBR_close.index(max(w_BBR_close))

max_wpd.append((w_Apo_open[index_Apo_open-1] + w_Apo_open[index_Apo_open] + w_Apo_open[index_Apo_open+1])/3)
max_wpd.append((w_Apo_close[index_Apo_close-1] + w_Apo_close[index_Apo_close])/2)
max_wpd.append((w_AD_open[index_AD_open-1] + w_AD_open[index_AD_open] + w_AD_open[index_AD_open+1])/3)
max_wpd.append((w_AD_close[index_AD_close-1] + w_AD_close[index_AD_close] + w_AD_close[index_AD_close+1])/3)
max_wpd.append((w_BBR_open[index_BBR_open-1] + w_BBR_open[index_BBR_open] + w_BBR_open[index_BBR_open+1])/3)
max_wpd.append((w_BBR_close[index_BBR_close-1] + w_BBR_close[index_BBR_close])/2)

name = ['Apo Open', 'Apo Closed', 'AD Open', 'AD Closed', 'BBR Open', 'BBR Closed']
num = np.array([1, 2, 3, 4, 5, 6])
fig4 = plt.figure()
ax1 = fig4.add_subplot(111)
ax1.set_title("Comparison of WPD RMSF Peak")
ax1.set_ylabel('RMSF(nm)')
ax1.bar(num, max_wpd)
plt.xticks(num, name)
leg = ax1.legend()
fig4.savefig('WPD_peak_compare.png')

#Substrate Binding loop
res_rl = np.array([112, 113, 114, 115, 116, 117])
rl_Apo_open, rl_Apo_close, rl_AD_open, rl_AD_close, rl_BBR_open, rl_BBR_close = [],[],[],[],[],[] #Empty array for RMSF of WPD residues
for i in res_rl:
    j = i-1
    k = i-2
    rl_Apo_open.append(np.mean([r_a7[j], r_dis7[j], r_1sug_dis7[k], r_dis9[j], r_1sug_dis9[k]]))
    rl_Apo_close.append(np.mean([r_1sug[k], r_dis11[j], r_1sug_dis11[k]]))
    rl_AD_open.append(np.mean([r_1sug_dis11_AD[k], r_1sug_dis_alt_AD[k], r_1sug_dis_alt2_AD[k]]))
    rl_AD_close.append(np.mean([r_dis11_AD[j]]))
    rl_BBR_open.append(np.mean([r_BBR_a7[j], r_BBR_dis9[j]]))
    rl_BBR_close.append(np.mean([r_BBR_1sug[k], r_BBR_1sug_dis7[k]]))
res_rl = np.array([113, 114, 115, 116, 117, 118])

#Plot substrate binding loop RMSF
figr = plt.figure()
ax1 = figr.add_subplot(111)
ax1.set_title("Comparison of Substrate Binding Loop RMSF")    
ax1.set_xlabel('Residues')
ax1.set_ylabel('RMSF(nm)')
ax1.plot(res_rl,rl_Apo_open, label='Apo Open')
ax1.plot(res_rl,rl_Apo_close, label='Apo Closed')
ax1.plot(res_rl,rl_AD_open, label='AD')
#ax1.plot(res_rl,rl_AD_close, label='AD Closed')
ax1.plot(res_rl,rl_BBR_open, label='BBR')
#ax1.plot(res_rl,rl_BBR_close, label='BBR Closed')
leg = ax1.legend(loc='upper right')
figr.savefig('R_RMSF_compare.png') 


#a3 RMAF
res_a3 = np.array([186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199, 200]) #Residues in WPD loop
a3_Apo_open, a3_Apo_close, a3_AD_open, a3_AD_close, a3_BBR_open, a3_BBR_close = [],[],[],[],[],[] #Empty array for RMSF of WPD residues
for i in res_a3:
    j = i-1
    k = i-2
    a3_Apo_open.append(np.mean([r_a7[j], r_dis7[j], r_1sug_dis7[k], r_dis9[j], r_1sug_dis9[k]]))
    a3_Apo_close.append(np.mean([r_1sug[k], r_dis11[j], r_1sug_dis11[k]]))
    a3_AD_open.append(np.mean([r_1sug_dis11_AD[k], r_1sug_dis_alt_AD[k], r_1sug_dis_alt2_AD[k]]))
    a3_AD_close.append(np.mean([r_dis11_AD[j]]))
    a3_BBR_open.append(np.mean([r_BBR_a7[j], r_BBR_dis9[j]]))
    a3_BBR_close.append(np.mean([r_BBR_1sug[k], r_BBR_1sug_dis7[k]]))

#Plot a3 RMSF
fig2 = plt.figure()
ax1 = fig2.add_subplot(111)
ax1.set_title("Comparison of a3 helix RMSF")    
ax1.set_xlabel('Residues')
ax1.set_ylabel('RMSF(nm)')
ax1.plot(res_a3,a3_Apo_open, label='Apo Open')
ax1.plot(res_a3,a3_Apo_close, label='Apo Closed')
ax1.plot(res_a3,a3_AD_open, label='AD Open')
ax1.plot(res_a3,a3_AD_close, label='AD Closed')
ax1.plot(res_a3,a3_BBR_open, label='BBR Open')
ax1.plot(res_a3,a3_BBR_close, label='BBR Closed')
leg = ax1.legend(loc='upper right')
fig2.savefig('a3_RMSF_compare.png') 

#a6 RMSF
res_a6 = np.array([264, 265, 266, 267, 268, 269, 270, 271, 272, 273, 274, 275, 276, 277, 278, 279, 280, 281]) #Residues in WPD loop
a6_Apo_open, a6_Apo_close, a6_AD_open, a6_AD_close, a6_BBR_open, a6_BBR_close = [],[],[],[],[],[] #Empty array for RMSF of WPD residues
for i in res_a6:
    j = i-1
    k = i-2
    a6_Apo_open.append(np.mean([r_a7[j], r_dis7[j], r_1sug_dis7[k], r_dis9[j], r_1sug_dis9[k]]))
    a6_Apo_close.append(np.mean([r_1sug[k], r_dis11[j], r_1sug_dis11[k]]))
    a6_AD_open.append(np.mean([r_1sug_dis11_AD[k], r_1sug_dis_alt_AD[k], r_1sug_dis_alt2_AD[k]]))
    a6_AD_close.append(np.mean([r_dis11_AD[j]]))
    a6_BBR_open.append(np.mean([r_BBR_a7[j], r_BBR_dis9[j]]))
    a6_BBR_close.append(np.mean([r_BBR_1sug[k], r_BBR_1sug_dis7[k]]))

#Plot a6 RMSF
fig3 = plt.figure()
ax1 = fig3.add_subplot(111)
ax1.set_title("Comparison of a6 helix RMSF")    
ax1.set_xlabel('Residues')
ax1.set_ylabel('RMSF(nm)')
ax1.plot(res_a6,a6_Apo_open, label='Apo Open')
ax1.plot(res_a6,a6_Apo_close, label='Apo Closed')
ax1.plot(res_a6,a6_AD_open, label='AD Open')
ax1.plot(res_a6,a6_AD_close, label='AD Closed')
ax1.plot(res_a6,a6_BBR_open, label='BBR Open')
ax1.plot(res_a6,a6_BBR_close, label='BBR Closed')
leg = ax1.legend(loc='upper right')
fig3.savefig('a6_RMSF_compare.png') 

#a7 RMSF
res_a7 = np.array([285, 286, 287, 288, 289, 290, 291, 292, 293, 294, 295, 296, 297]) #Residues in WPD loop
a7_Apo_open, a7_Apo_close, a7_AD_open, a7_AD_close, a7_AD_unb, a7_BBR_open, a7_BBR_close, a7_BBR_unb = [],[],[],[],[],[],[],[] #Empty array for RMSF of WPD residues
a7_Apo_open_sem, a7_Apo_close_sem, a7_AD_open_sem, a7_AD_close_sem, a7_AD_unb_sem, a7_BBR_open_sem, a7_BBR_close_sem, a7_BBR_unb_sem = [],[],[],[],[],[],[],[] #Empty array for Error for RMSF of WPD residues
a7_unb1, a7_unb2, a7_unb3, a7_unb4, a7_unb5, a7_unb6, a7_unb7 = [],[],[],[],[],[],[]
a7_dis7, a7_dis9, a7_dis11, a7_dis_unb, a7_ord = [],[],[],[],[]
a7_lig_open, a7_lig_open_sem = [],[]
a7_lig_unb, a7_lig_unb_sem = [],[]
for i in res_a7:
    j = i-1
    k = i-2
    a7_Apo_open.append(np.mean([r_a7[j], r_dis7[j], r_1sug_dis7[k], r_dis9[j], r_1sug_dis9[k]]))
    a7_Apo_close.append(np.mean([r_1sug[k], r_dis11[j], r_1sug_dis11[k]]))
    a7_AD_open.append(np.mean([r_1sug_dis11_AD[k], r_1sug_dis_alt_AD[k], r_1sug_dis_alt2_AD[k]]))
    a7_AD_close.append(np.mean([r_dis11_AD[j]]))
    a7_AD_unb.append(np.mean([r_dis7_AD[j], r_1sug_dis7_AD[k], r_1sug_dis9_AD[k], r_1sug_dis11_2AD[k]]))
    a7_BBR_open.append(np.mean([r_BBR_a7[j], r_BBR_dis9[j]]))
    a7_BBR_close.append(np.mean([r_BBR_1sug[k], r_BBR_1sug_dis7[k]]))
    a7_BBR_unb.append(np.mean([r_BBR_dis7[j], r_BBR_dis11[j], r_BBR_1sug_dis11[k]]))
    a7_lig_open.append(np.mean([r_1sug_dis11_AD[k], r_1sug_dis_alt_AD[k], r_1sug_dis_alt2_AD[k], r_BBR_a7[j], r_BBR_dis9[j]]))
    a7_lig_unb.append(np.mean([r_dis7_AD[j], r_1sug_dis7_AD[k], r_1sug_dis9_AD[k], r_1sug_dis11_2AD[k], r_BBR_dis7[j], r_BBR_dis11[j], r_BBR_1sug_dis11[k]]))

    a7_unb1.append(r_dis7_AD[j])
    a7_unb2.append(r_1sug_dis7_AD[k])
    a7_unb3.append(r_1sug_dis9_AD[k])
    a7_unb4.append(r_1sug_dis11_2AD[k])
    a7_unb5.append(r_BBR_dis7[j])
    a7_unb6.append(r_BBR_dis11[j])
    a7_unb7.append(r_BBR_1sug_dis11[k])
    a7_dis7.append(np.mean([r_dis7[j], r_1sug_dis7[k]]))
    a7_dis9.append(np.mean([r_dis9[j], r_1sug_dis9[k]]))
    a7_dis_unb.append(np.mean([r_dis7[j], r_1sug_dis7[k], r_dis9[j], r_1sug_dis9[k]]))
    a7_dis11.append(np.mean([r_dis11[j], r_1sug_dis11[k]]))
    a7_ord.append(np.mean([r_a7[j], r_1sug[k]]))

    a7_Apo_open_sem.append(stats.sem([r_a7[j], r_dis7[j], r_1sug_dis7[k], r_dis9[j], r_1sug_dis9[k]]))
    a7_Apo_close_sem.append(stats.sem([r_1sug[k], r_dis11[j], r_1sug_dis11[k]]))
    a7_AD_open_sem.append(stats.sem([r_1sug_dis11_AD[k], r_1sug_dis_alt_AD[k], r_1sug_dis_alt2_AD[k]]))
    a7_AD_close_sem.append(0)
    a7_AD_unb_sem.append(stats.sem([r_dis7_AD[j], r_1sug_dis7_AD[k], r_1sug_dis9_AD[k], r_1sug_dis11_2AD[k]]))
    a7_BBR_open_sem.append(stats.sem([r_BBR_a7[j], r_BBR_dis9[j]]))
    a7_BBR_close_sem.append(stats.sem([r_BBR_1sug[k], r_BBR_1sug_dis7[k]]))
    a7_BBR_unb_sem.append(stats.sem([r_BBR_dis7[j], r_BBR_dis11[j], r_BBR_1sug_dis11[k]]))
    a7_lig_open_sem.append(stats.sem([r_1sug_dis11_AD[k], r_1sug_dis_alt_AD[k], r_1sug_dis_alt2_AD[k], r_BBR_a7[j], r_BBR_dis9[j]]))
    a7_lig_unb_sem.append(stats.sem([r_dis7_AD[j], r_1sug_dis7_AD[k], r_1sug_dis9_AD[k], r_1sug_dis11_2AD[k], r_BBR_dis7[j], r_BBR_dis11[j], r_BBR_1sug_dis11[k]]))

#Plot a7 RMSF
fig4 = plt.figure()
ax1 = fig4.add_subplot(111)
ax1.set_title("Comparison of a7 helix RMSF")    
ax1.set_xlabel('Residues')
ax1.set_ylabel('RMSF(nm)')
ax1.errorbar(res_a7, a7_Apo_open, yerr = a7_Apo_open_sem, label='Apo Open')
ax1.errorbar(res_a7, a7_Apo_close, yerr = a7_Apo_close_sem, label='Apo Closed')
ax1.errorbar(res_a7, a7_AD_open, yerr = a7_AD_open_sem, label='AD')
#ax1.errorbar(res_a7, a7_AD_close, yerr = a7_AD_close_sem, label='AD Closed')
ax1.errorbar(res_a7, a7_BBR_open, yerr = a7_BBR_open_sem, label='BBR')
#ax1.errorbar(res_a7, a7_BBR_close, yerr = a7_BBR_close_sem, label='BBR Closed')
ax1.legend(loc='upper right')
fig4.savefig('a7_RMSF_compare_err.png') 

#Plot a7 RMSF
fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.set_title("Comparison of a7 helix RMSF")    
ax1.set_xlabel('Residues')
ax1.set_ylabel('RMSF(nm)')
ax1.plot(res_a7, a7_Apo_open, label='Apo Open', color = 'blue')
ax1.plot(res_a7, a7_Apo_close, label='Apo Closed', color = 'orange')
ax1.legend(loc='upper right')
fig.savefig('a7_RMSF_Apo.png') 
plt.close(fig)

fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.set_title("Comparison of a7 helix RMSF")    
ax1.set_xlabel('Residues')
ax1.set_ylabel('RMSF(nm)')
ax1.plot(res_a7, a7_Apo_open, label='Apo Open', color = 'blue')
ax1.plot(res_a7, a7_Apo_close, label='Apo Closed', color = 'orange')
ax1.plot(res_a7, a7_AD_open, label='AD', color = 'purple')
ax1.legend(loc='upper right')
fig.savefig('a7_RMSF_AD.png') 
plt.close(fig)

fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.set_title("Comparison of a7 helix RMSF")    
ax1.set_xlabel('Residues')
ax1.set_ylabel('RMSF(nm)')
ax1.plot(res_a7, a7_Apo_open, label='Apo Open', color = 'blue')
ax1.plot(res_a7, a7_AD_open, label='AD', color = 'blue', linestyle = 'dashed')
ax1.plot(res_a7, a7_BBR_open, label='BBR', color = 'blue', linestyle = 'dotted')
ax1.plot(res_a7, a7_Apo_close, label='Apo Closed', color = 'red')
ax1.legend(loc='best')
fig.savefig('a7_RMSF_compare.png') 

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
ax1.set_title("Comparison of a7 helix RMSF")    
ax1.set_xlabel('Residues')
ax1.set_ylabel('RMSF(nm)')
ax1.plot(res_a7, a7_Apo_open, label='Apo Open', color = 'blue')
ax1.fill_between(res_a7, a7_Apo_open-a7_Apo_open_sem, a7_Apo_open+a7_Apo_open_sem, alpha=0.3, facecolor = 'blue', edgecolor = 'blue')
ax1.plot(res_a7, a7_AD_open, label='AD', color = 'blue', linestyle = 'dashed')
ax1.fill_between(res_a7, a7_AD_open-a7_AD_open_sem, a7_AD_open+a7_AD_open_sem, alpha=0.25, facecolor = 'blue', edgecolor = 'blue')
ax1.plot(res_a7, a7_BBR_open, label='BBR', color = 'blue', linestyle = 'dotted')
ax1.fill_between(res_a7, a7_BBR_open-a7_BBR_open_sem, a7_BBR_open+a7_BBR_open_sem, alpha=0.2, facecolor = 'blue', edgecolor = 'blue')
ax1.plot(res_a7, a7_Apo_close, label='Apo Closed', color = 'red')
ax1.fill_between(res_a7, a7_Apo_close-a7_Apo_close_sem, a7_Apo_close+a7_Apo_close_sem, alpha=0.3, facecolor = 'red', edgecolor = 'red')
ax1.legend(loc='best')
fig.savefig('a7_RMSF_compare_err.png') 

#Plot a7 Unbound RMSF
a7_lig_open = np.array(a7_lig_open)
a7_lig_open_sem = np.array(a7_lig_open_sem)
a7_lig_unb = np.array(a7_lig_unb)
a7_lig_unb_sem = np.array(a7_lig_unb_sem)

fig4 = plt.figure()
ax1 = fig4.add_subplot(111)
ax1.set_title("Comparison of alpha-7 helix Flexibility", fontsize = 16)
ax1.set_xlabel('Residue ID', fontsize = 12)
ax1.set_ylabel('RMSF(nm)', fontsize = 12)
ax1.plot(res_a7, a7_lig_open, label='Ligand Binds', color = 'blue')
ax1.fill_between(res_a7, a7_lig_open-a7_lig_open_sem, a7_lig_open+a7_lig_open_sem, alpha=0.5, facecolor = 'blue', edgecolor = 'blue')
ax1.plot(res_a7, a7_lig_unb, label='Ligand Does Not Bind', color = 'red')
ax1.fill_between(res_a7, a7_lig_unb-a7_lig_unb_sem, a7_lig_unb+a7_lig_unb_sem, alpha=0.5, facecolor = 'red', edgecolor = 'red')
ax1.legend(loc='upper left', fontsize = 12)
fig4.savefig('a7_RMSF_lig_unbound_err.png') 
plt.close(fig4)

fig4 = plt.figure()
ax1 = fig4.add_subplot(111)
ax1.set_title("Comparison of alpha-7 helix Flexibility", fontsize = 16)
ax1.set_xlabel('Residue ID', fontsize = 12)
ax1.set_ylabel('RMSF(nm)', fontsize = 12)
ax1.plot(res_a7, a7_AD_open, label='AD Binds', color = 'blue')
ax1.plot(res_a7, a7_AD_unb, label='AD Does Not Bind', color = 'blue', linestyle='dashed')
ax1.plot(res_a7, a7_BBR_open, label='BBR Binds', color = 'purple')
ax1.plot(res_a7, a7_BBR_unb, label='BBR Does Not Bind', color = 'purple', linestyle='dashed')
ax1.legend(loc='upper left', fontsize = 12)
fig4.savefig('a7_RMSF_unbound_compare.png') 

#Plot a7 Unbound vs Bound
a7_dis = [a7_AD_open, a7_AD_unb, a7_BBR_open, a7_BBR_unb]
fig, ax2 = plt.subplots(figsize=(18, 8))
ax2.set_title('Customized violin plot')
parts = ax2.violinplot(a7_dis, showmeans=False, showmedians=False, showextrema=False)

for pc in parts['bodies']:
    pc.set_facecolor('#D43F3A')
    pc.set_edgecolor('black')
    pc.set_alpha(1)

quartile1, medians, quartile3 = np.percentile(a7_dis, [25, 50, 75], axis=1)
whiskers = np.array([
    adjacent_values(sorted_array, q1, q3)
    for sorted_array, q1, q3 in zip(a7_dis, quartile1, quartile3)])
whiskers_min, whiskers_max = whiskers[:, 0], whiskers[:, 1]

inds = np.arange(1, len(medians) + 1)
ax2.scatter(inds, medians, marker='o', color='white', s=30, zorder=3)
ax2.vlines(inds, quartile1, quartile3, color='k', linestyle='-', lw=5)
ax2.vlines(inds, whiskers_min, whiskers_max, color='k', linestyle='-', lw=1)

# set style for the axes
labels = ['AD Binds', 'AD Does Not Bind', 'BBR Binds', 'BBR Does Not Bind']
set_axis_style(ax2, labels)
fig.savefig('a7_avg_RMSF_binding_compare.png') 

#All unbound structures
fig3 = plt.figure()
ax1 = fig3.add_subplot(111)
ax1.set_title("Comparison of a7 helix RMSF")    
ax1.set_xlabel('Residues')
ax1.set_ylabel('RMSF(nm)')
ax1.plot(res_a7,a7_unb1, label='config 7')
ax1.plot(res_a7,a7_unb2, label='config 7')
ax1.plot(res_a7,a7_unb3, label='config 9')
#ax1.plot(res_a7,a7_unb4, label='config 11')
ax1.plot(res_a7,a7_unb5, label='config 7')
#ax1.plot(res_a7,a7_unb6, label='config 11')
#ax1.plot(res_a7,a7_unb7, label='config 11')
leg = ax1.legend(loc='upper right')
fig3.savefig('a7_RMSF_unbound_all_compare.png') 

#All unbound structures
fig3 = plt.figure()
ax1 = fig3.add_subplot(111)
ax1.set_title("Comparison of a7 helix RMSF for Apo Structures")    
ax1.set_xlabel('Residues')
ax1.set_ylabel('RMSF(nm)')
ax1.plot(res_a7,a7_dis_unb, label='Ligand Does not Bind')
ax1.plot(res_a7,a7_dis11, label='Ligand Binds')
ax1.plot(res_a7,a7_ord, label='Ligand Binds Alternative Location')
leg = ax1.legend(loc='upper right')
fig3.savefig('a7_RMSF_Apo_dis_compare.png') 

