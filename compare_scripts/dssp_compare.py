#!/ usr / bin / env python

#Import packages
from matplotlib import pyplot as plt
import numpy as np
from scipy import stats

def per_helx(dssp):
    char_num = np.arange(2,15,2)

    num = 0
    time_tot = int(len(dssp))
    
    alpha = np.zeros(len(char_num))
    struct = np.zeros(len(char_num))

    per_ind_tot = round(time_tot/10)
    alpha_per = np.zeros([len(char_num), per_ind_tot])
    struct_per = np.zeros([len(char_num), per_ind_tot])

    for i in dssp:
        char = i
        c = 0
        #Determine DSSP Values
        for n in char_num:
            if char[n]=='H':
                alpha[c] += 1
            if char[n] == 'T' or char[0] == 'G' or char[0] == 'I':
                struct[c] += 1
            
            #Every 20 time steps take a running percentage
            if num % 10 == 0 and num != 0 and num < per_ind_tot*10:
                t = int(num/10)
                alpha_per[c][t] = 100 * alpha[c] / 10
                struct_per[c][t] = 100 * struct[c] / 10
                alpha[c] = 0
                struct[c] = 0
            c += 1
        #Iterate time step
        num += 1
    #Determine overall percent for each residue
    alpha_per_mean = np.zeros(len(char_num))
    alpha_per_sem = np.zeros(len(char_num))
    struct_per_mean = np.zeros(len(char_num))
    struct_per_sem = np.zeros(len(char_num))

    for i in range(len(char_num)):
        alpha_per_mean[i] = np.mean(alpha_per[i][:])
        struct_per_mean[i] = np.mean(struct_per[i][:])
        alpha_per_sem[i] = stats.sem(alpha_per[i][:])
        struct_per_sem[i] = stats.sem(struct_per[i][:])    
    return alpha_per, struct_per, alpha_per_mean, struct_per_mean, alpha_per_sem, struct_per_sem

def time_avg(alpha_per, struct_per, eq_time, tot_time):
    #Determine number of characters and time indices
    res, ind = np.shape(alpha_per)

    #Set time vector
    time = np.linspace(eq_time, tot_time, num=ind)

    #Avg over residues for each time point
    alpha_t = np.zeros(ind)
    alpha_t_err = np.zeros(ind)

    struct_t = np.zeros(ind)
    struct_t_err = np.zeros(ind)

    for i in range(ind):
        alpha_t[i] = np.mean(alpha_per[:,i])
        alpha_t_err[i] = stats.sem(alpha_per[:,i])
        struct_t[i] = np.mean(struct_per[:,i])
        struct_t_err[i] = stats.sem(struct_per[:,i])

    return time, alpha_t, alpha_t_err, struct_t, struct_t_err

#Open files listing H_bonds in read-only format
list_a7 = open('../../rebuild_a7/analysis/DSSP_a7.txt','r').readlines() #WPD open, a7 ordered, no AD
list_a7_AD = open('../../AD_rebuild_a7/analysis/DSSP_a7_AD.txt', 'r') .readlines()#WPD open, a7 ordered, AD present
list_1sug = open('../../Apo_1SUG/analysis/1sug/DSSP_1sug.txt','r').readlines() #WPD closed, a7 ordered, no AD
list_1sug_AD = open ('../../1sug_AD/analysis/DSSP_1sug_AD.txt','r').readlines() #WPD closed/open, a7 ordered, AD present
list_1sug_dis_AD_7 = open ('../../1sug_dis_AD/analysis/config7/DSSP_1sug_dis_AD7.txt','r').readlines() #WPD closed/open, a7 disordered, AD present
list_1sug_dis_AD_9 = open ('../../1sug_dis_AD/analysis/config9/DSSP_1sug_dis_AD9.txt','r').readlines() #WPD closed, a7 disordered, AD present
list_1sug_dis_AD_11 = open ('../../1sug_dis_AD/analysis/config11/DSSP_1sug_dis_AD11.txt','r').readlines() #WPD closed, a7 disordered, AD present
list_1sug_dis_AD_11_2 = open ('../../1sug_dis_AD/analysis/config11_2/DSSP_1sug_dis_AD11_2.txt','r').readlines() #WPD closed, a7 disordered, AD present
list_1sug_dis_AD_alt = open ('../../1sug_dis_AD/analysis/config_alt/DSSP_1sug_AD_alt.txt','r').readlines() #WPD closed, a7 disordered, AD present
list_1sug_dis_AD_alt2 = open ('../../1sug_dis_AD/analysis/config_alt2/DSSP_alt2.txt','r').readlines() #WPD closed, a7 disordered, AD present
list_AD_dis_7 = open ('../../AD_dis/analysis/config7/DSSP_AD_dis7.txt','r').readlines() #WPD open/closed, a7 disordered, AD present
list_AD_dis_9 = open ('../../AD_dis/analysis/config9/DSSP_AD_dis9.txt','r').readlines() #WPD open/closed, a7 disordered, AD present
list_AD_dis_11 = open ('../../AD_dis/analysis/config11/DSSP_AD_dis11.txt','r').readlines() #WPD open/closed, a7 disordered, AD present
list_1sug_dis_7 = open ('../../1sug_dis/analysis/config7/DSSP_1sug_dis7.txt','r').readlines() #WPD closed/open, a7 disordered, no AD
list_1sug_dis_9 = open ('../../1sug_dis/analysis/config9/DSSP_1sug_dis9.txt','r').readlines() #WPD closed, a7 disordered, no AD
list_1sug_dis_11 = open ('../../1sug_dis/analysis/config11/DSSP_1sug_dis11.txt','r').readlines() #WPD closed, a7 disordered, no AD
list_dis_7 = open ('../../rebuild_a7_high/config7/analysis/DSSP_Apo_dis7.txt','r').readlines()#WPD open/closed, a7 disordered, no AD
list_dis_9 = open ('../../rebuild_a7_high/config9/analysis/DSSP_Apo_dis9.txt','r').readlines() #WPD open/closed, a7 disordered, no AD
list_dis_11 = open ('../../rebuild_a7_high/config11/analysis/DSSP_Apo_dis11.txt','r').readlines() #WPD open/closed, a7 disordered, no AD
list_BBR_a7 = open('../../BBR_a7/analysis/DSSP_BBR_a7.txt', 'r').readlines() #WPD open, a7 ordered, BBR present
list_BBR_1sug = open('../../BBR_1sug/analysis/DSSP_BBR_1sug.txt', 'r').readlines() #WPD close, a7 ordered, BBR present


#Seperate Characters in the string and record number that are in an alpha helix
a7_alpha_per, a7_struct_per, a7_alpha_mean, a7_struct_mean, a7_alpha_sem, a7_struct_sem = per_helx(list_a7)
a7_AD_alpha_per, a7_AD_struct_per, a7_AD_alpha_mean, a7_AD_struct_mean, a7_AD_alpha_sem, a7_AD_struct_sem = per_helx(list_a7_AD)
sug_alpha_per, sug_struct_per, sug_alpha_mean, sug_struct_mean, sug_alpha_sem, sug_struct_sem = per_helx(list_1sug)
sug_AD_alpha_per, sug_AD_struct_per, sug_AD_alpha_mean, sug_AD_struct_mean, sug_AD_alpha_sem, sug_AD_struct_sem = per_helx(list_1sug_AD)
sug_AD7_alpha_per, sug_AD7_struct_per, sug_AD7_alpha_mean, sug_AD7_struct_mean, sug_AD7_alpha_sem, sug_AD7_struct_sem = per_helx(list_1sug_dis_AD_7)
sug_AD9_alpha_per, sug_AD9_struct_per, sug_AD9_alpha_mean, sug_AD9_struct_mean, sug_AD9_alpha_sem, sug_AD9_struct_sem = per_helx(list_1sug_dis_AD_9)
sug_AD11_alpha_per, sug_AD11_struct_per, sug_AD11_alpha_mean, sug_AD11_struct_mean, sug_AD11_alpha_sem, sug_AD11_struct_sem = per_helx(list_1sug_dis_AD_11)
sug_AD11_2_alpha_per, sug_AD11_2_struct_per, sug_AD11_2_alpha_mean, sug_AD11_2_struct_mean, sug_AD11_2_alpha_sem, sug_AD11_2_struct_sem = per_helx(list_1sug_dis_AD_11_2)
sug_AD_alt_alpha_per, sug_AD_alt_struct_per, sug_AD_alt_alpha_mean, sug_AD_alt_struct_mean, sug_AD_alt_alpha_sem, sug_AD_alt_struct_sem = per_helx(list_1sug_dis_AD_alt)
sug_AD_alt2_alpha_per, sug_AD_alt2_struct_per, sug_AD_alt2_alpha_mean, sug_AD_alt2_struct_mean, sug_AD_alt2_alpha_sem, sug_AD_alt2_struct_sem = per_helx(list_1sug_dis_AD_alt2)
AD7_alpha_per, AD7_struct_per, AD7_alpha_mean, AD7_struct_mean, AD7_alpha_sem, AD7_struct_sem = per_helx(list_AD_dis_7)
AD9_alpha_per, AD9_struct_per, AD9_alpha_mean, AD9_struct_mean, AD9_alpha_sem, AD9_struct_sem = per_helx(list_AD_dis_9)
AD11_alpha_per, AD11_struct_per, AD11_alpha_mean, AD11_struct_mean, AD11_alpha_sem, AD11_struct_sem = per_helx(list_AD_dis_11)
sug7_alpha_per, sug7_struct_per, sug7_alpha_mean, sug7_struct_mean, sug7_alpha_sem, sug7_struct_sem = per_helx(list_1sug_dis_7)
sug9_alpha_per, sug9_struct_per, sug9_alpha_mean, sug9_struct_mean, sug9_alpha_sem, sug9_struct_sem = per_helx(list_1sug_dis_9)
sug11_alpha_per, sug11_struct_per, sug11_alpha_mean, sug11_struct_mean, sug11_alpha_sem, sug11_struct_sem = per_helx(list_1sug_dis_11)
BBR_a7_alpha_per, BBR_a7_struct_per, BBR_a7_alpha_mean, BBR_a7_struct_mean, BBR_a7_alpha_sem, BBR_a7_struct_sem = per_helx(list_BBR_a7)
BBR_1sug_alpha_per, BBR_1sug_struct_per, BBR_1sug_alpha_mean, BBR_1sug_struct_mean, BBR_1sug_alpha_sem, BBR_1sug_struct_sem = per_helx(list_BBR_1sug)

#Plot comparing degree of helicity between corresponding Apo and AD bound structures
num = np.linspace(287, 294, num=7)
fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.set_title('Degree of Helicity')    
ax1.set_ylabel('% Residue was in Alpha Helix')
ax1.set_ylim(0,100)
ax1.plot(num, a7_alpha_mean, label='Open-Apo', color = 'blue')
ax1.plot(num, a7_AD_alpha_mean, label='Open-AD', color = 'orange')
ax1.plot(num, BBR_a7_alpha_mean, label='Open-BBR', color = 'red')
ax1.plot(num, sug_alpha_mean, label='Close-Apo', color = 'green')
ax1.plot(num, sug_AD_alpha_mean, label='Close/Open-AD', color = 'blue')
leg = ax1.legend()
fig.savefig('Helicity_disordering.png')
plt.close(fig)

#Average AD open and AD closed/open
fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.set_title('Degree of Helicity')    
ax1.set_ylabel('% Residue was in Alpha Helix')
ax1.set_ylim(0,100)
ax1.plot(num, a7_alpha_mean, label='Apo Open', color = 'blue')
ax1.plot(num, a7_AD_alpha_mean, label='AD', linestyle = 'dashed', color = 'blue')
ax1.plot(num, BBR_a7_alpha_mean, label='BBR', linestyle = 'dotted', color = 'blue')
ax1.plot(num, sug_alpha_mean, label='Apo Closed', color = 'red')
ax1.set_xlabel('Residue ID')
leg = ax1.legend()
fig.savefig('Helicity_disordering_simp.png')
plt.close(fig)

fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.set_title('Degree of Helicity')    
ax1.set_ylabel('% Residue was in Alpha Helix')
ax1.set_ylim(0,100)
ax1.plot(num, a7_alpha_mean, label='Apo Open', color = 'blue')
ax1.fill_between(num, a7_alpha_mean - a7_alpha_sem, a7_alpha_mean + a7_alpha_sem, alpha=0.3, facecolor = 'blue', edgecolor = 'blue')
ax1.plot(num, a7_AD_alpha_mean, label='AD', linestyle = 'dashed', color = 'blue')
ax1.fill_between(num, a7_AD_alpha_mean - a7_AD_alpha_sem, a7_AD_alpha_mean + a7_AD_alpha_sem, alpha=0.3, facecolor = 'blue', edgecolor = 'blue')
ax1.plot(num, BBR_a7_alpha_mean, label='BBR', linestyle = 'dotted', color = 'blue')
ax1.fill_between(num, BBR_a7_alpha_mean - BBR_a7_alpha_sem, BBR_a7_alpha_mean + BBR_a7_alpha_sem, alpha=0.3, facecolor = 'blue', edgecolor = 'blue')
ax1.plot(num, sug_alpha_mean, label='Apo Closed', color = 'red')
ax1.fill_between(num, sug_alpha_mean - sug_alpha_sem, sug_alpha_mean + sug_alpha_sem, alpha=0.3, facecolor = 'red', edgecolor = 'red')
ax1.set_xlabel('Residue ID')
leg = ax1.legend()
fig.savefig('Helicity_disordering_simp_err.png')
plt.close(fig)

#Compute time averages
time_1sug, alpha_1sug, alpha_err_1sug, struct_1sug, struct_err_1sug = time_avg(sug_alpha_per, sug_struct_per, 5, 300)
time_a7, alpha_a7, alpha_err_a7, struct_a7, struct_err_a7 = time_avg(a7_alpha_per, a7_struct_per, 50, 300)
time_1sug_AD, alpha_1sug_AD, alpha_err_1sug_AD, struct_1sug_AD, struct_err_1sug_AD = time_avg(sug_AD_alpha_per, sug_AD_struct_per, 80, 300)
time_a7_AD, alpha_a7_AD, alpha_err_a7_AD, struct_a7_AD, struct_err_a7_AD = time_avg(a7_AD_alpha_per, a7_AD_struct_per, 80, 300)
time_BBR_a7, alpha_BBR_a7, alpha_err_BBR_a7, struct_BBR_a7, struct_err_BBR_a7 = time_avg(BBR_a7_alpha_per, BBR_a7_struct_per, 70, 300)
time_BBR_1sug, alpha_BBR_1sug, alpha_err_BBR_1sug, struct_BBR_1sug, struct_err_BBR_1sug = time_avg(BBR_1sug_alpha_per, BBR_1sug_struct_per, 5, 300)

fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.set_title('Degree of Helicity')    
ax1.set_ylabel('% Residue was in Alpha Helix')
ax1.set_ylim(0,100)
ax1.set_xlim(80,300)
ax1.plot(time_a7, alpha_a7, label='Apo Open', color = 'blue')
ax1.fill_between(time_a7, alpha_a7 - alpha_err_a7, alpha_a7 + alpha_err_a7, alpha=0.3, facecolor = 'blue', edgecolor = 'blue')
ax1.plot(time_a7_AD, alpha_a7_AD, label='AD', linestyle = 'dashed', color = 'blue')
ax1.fill_between(time_a7_AD, alpha_a7_AD - alpha_err_a7_AD, alpha_a7_AD + alpha_err_a7_AD, alpha=0.3, facecolor = 'blue', edgecolor = 'blue')
ax1.plot(time_BBR_a7, alpha_BBR_a7, label='BBR', linestyle = 'dotted', color = 'blue')
ax1.fill_between(time_BBR_a7, alpha_BBR_a7 - alpha_err_BBR_a7, alpha_BBR_a7 + alpha_err_BBR_a7, alpha=0.3, facecolor = 'blue', edgecolor = 'blue')
ax1.plot(time_1sug, alpha_1sug, label='Apo Closed', color = 'red')
ax1.fill_between(time_1sug, alpha_1sug - alpha_err_1sug, alpha_1sug + alpha_err_1sug, alpha=0.3, facecolor = 'red', edgecolor = 'red')
ax1.set_xlabel('Time(ns)')
leg = ax1.legend()
fig.savefig('Helicity_disordering_time_err.png')
plt.close(fig)

fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.set_title('Degree of Helicity')    
ax1.set_ylabel('% Residue was in Alpha Helix')
ax1.set_ylim(0,100)
ax1.set_xlim(80,300)
ax1.plot(time_a7, alpha_a7, label='Apo Open', color = 'blue')
ax1.plot(time_a7_AD, alpha_a7_AD, label='AD', linestyle = 'dashed', color = 'blue')
ax1.plot(time_BBR_a7, alpha_BBR_a7, label='BBR', linestyle = 'dotted', color = 'blue')
ax1.plot(time_1sug, alpha_1sug, label='Apo Closed', color = 'red')
ax1.set_xlabel('Time(ns)')
leg = ax1.legend()
fig.savefig('Helicity_disordering_time.png')
plt.close(fig)

