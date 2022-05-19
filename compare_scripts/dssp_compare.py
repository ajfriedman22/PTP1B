#!/ usr / bin / env python

#Import packages
from matplotlib import pyplot as plt
import numpy as np
from scipy import stats
import sys

#Import custom modules
sys.path.insert(1, '/ocean/projects/cts160011p/afriedma/code/PTP1B/util/')
import mdfunc
import plot

#Load file path
file_path = open('dssp_file_path.txt', 'r').readlines()

#Open files listing H_bonds in read-only format
list_a7 = open(file_path[0].strip(), 'r').readlines() #WPD open, a7 ordered, no AD
list_a7_AD = open(file_path[1].strip(), 'r') .readlines()#WPD open, a7 ordered, AD present
list_1sug = open(file_path[2].strip(),'r').readlines() #WPD closed, a7 ordered, no AD
list_1sug_AD = open (file_path[3].strip(),'r').readlines() #WPD closed/open, a7 ordered, AD present
list_1sug_dis_AD_7 = open (file_path[4].strip(),'r').readlines() #WPD closed/open, a7 disordered, AD present
list_1sug_dis_AD_9 = open (file_path[5].strip(),'r').readlines() #WPD closed, a7 disordered, AD present
list_1sug_dis_AD_11 = open (file_path[6].strip(),'r').readlines() #WPD closed, a7 disordered, AD present
list_1sug_dis_AD_11_2 = open (file_path[7].strip(),'r').readlines() #WPD closed, a7 disordered, AD present
list_1sug_dis_AD_alt = open (file_path[8].strip(),'r').readlines() #WPD closed, a7 disordered, AD present
list_1sug_dis_AD_alt2 = open (file_path[9].strip(),'r').readlines() #WPD closed, a7 disordered, AD present
list_AD_dis_7 = open (file_path[10].strip(),'r').readlines() #WPD open/closed, a7 disordered, AD present
list_AD_dis_9 = open (file_path[11].strip(),'r').readlines() #WPD open/closed, a7 disordered, AD present
list_AD_dis_11 = open (file_path[12].strip(),'r').readlines() #WPD open/closed, a7 disordered, AD present
list_1sug_dis_7 = open (file_path[13].strip(),'r').readlines() #WPD closed/open, a7 disordered, no AD
list_1sug_dis_9 = open (file_path[14].strip(),'r').readlines() #WPD closed, a7 disordered, no AD
list_1sug_dis_11 = open (file_path[15].strip(),'r').readlines() #WPD closed, a7 disordered, no AD
list_dis_7 = open (file_path[16].strip(),'r').readlines()#WPD open/closed, a7 disordered, no AD
list_dis_9 = open (file_path[17].strip(),'r').readlines() #WPD open/closed, a7 disordered, no AD
list_dis_11 = open (file_path[18].strip(),'r').readlines() #WPD open/closed, a7 disordered, no AD
list_BBR_a7 = open(file_path[19].strip(), 'r').readlines() #WPD open, a7 ordered, BBR present
list_BBR_1sug = open(file_path[20].strip(), 'r').readlines() #WPD close, a7 ordered, BBR present
list_AD_BBR = open(file_path[21].strip(), 'r').readlines() #AD and BBR combo simulation

#Seperate Characters in the string and record number that are in an alpha helix
a7_alpha_per, a7_struct_per, a7_alpha_mean, a7_struct_mean, a7_alpha_sem, a7_struct_sem = mdfunc.per_helx(list_a7, False)
a7_AD_alpha_per, a7_AD_struct_per, a7_AD_alpha_mean, a7_AD_struct_mean, a7_AD_alpha_sem, a7_AD_struct_sem = mdfunc.per_helx(list_a7_AD, False)
sug_alpha_per, sug_struct_per, sug_alpha_mean, sug_struct_mean, sug_alpha_sem, sug_struct_sem = mdfunc.per_helx(list_1sug, False)
sug_AD_alpha_per, sug_AD_struct_per, sug_AD_alpha_mean, sug_AD_struct_mean, sug_AD_alpha_sem, sug_AD_struct_sem = mdfunc.per_helx(list_1sug_AD, False)
sug_AD7_alpha_per, sug_AD7_struct_per, sug_AD7_alpha_mean, sug_AD7_struct_mean, sug_AD7_alpha_sem, sug_AD7_struct_sem = mdfunc.per_helx(list_1sug_dis_AD_7, False)
sug_AD9_alpha_per, sug_AD9_struct_per, sug_AD9_alpha_mean, sug_AD9_struct_mean, sug_AD9_alpha_sem, sug_AD9_struct_sem = mdfunc.per_helx(list_1sug_dis_AD_9, False)
sug_AD11_alpha_per, sug_AD11_struct_per, sug_AD11_alpha_mean, sug_AD11_struct_mean, sug_AD11_alpha_sem, sug_AD11_struct_sem = mdfunc.per_helx(list_1sug_dis_AD_11, False)
sug_AD11_2_alpha_per, sug_AD11_2_struct_per, sug_AD11_2_alpha_mean, sug_AD11_2_struct_mean, sug_AD11_2_alpha_sem, sug_AD11_2_struct_sem = mdfunc.per_helx(list_1sug_dis_AD_11_2, False)
sug_AD_alt_alpha_per, sug_AD_alt_struct_per, sug_AD_alt_alpha_mean, sug_AD_alt_struct_mean, sug_AD_alt_alpha_sem, sug_AD_alt_struct_sem = mdfunc.per_helx(list_1sug_dis_AD_alt, False)
sug_AD_alt2_alpha_per, sug_AD_alt2_struct_per, sug_AD_alt2_alpha_mean, sug_AD_alt2_struct_mean, sug_AD_alt2_alpha_sem, sug_AD_alt2_struct_sem = mdfunc.per_helx(list_1sug_dis_AD_alt2, False)
AD7_alpha_per, AD7_struct_per, AD7_alpha_mean, AD7_struct_mean, AD7_alpha_sem, AD7_struct_sem = mdfunc.per_helx(list_AD_dis_7, False)
AD9_alpha_per, AD9_struct_per, AD9_alpha_mean, AD9_struct_mean, AD9_alpha_sem, AD9_struct_sem = mdfunc.per_helx(list_AD_dis_9, False)
AD11_alpha_per, AD11_struct_per, AD11_alpha_mean, AD11_struct_mean, AD11_alpha_sem, AD11_struct_sem = mdfunc.per_helx(list_AD_dis_11, False)
sug7_alpha_per, sug7_struct_per, sug7_alpha_mean, sug7_struct_mean, sug7_alpha_sem, sug7_struct_sem = mdfunc.per_helx(list_1sug_dis_7, False)
sug9_alpha_per, sug9_struct_per, sug9_alpha_mean, sug9_struct_mean, sug9_alpha_sem, sug9_struct_sem = mdfunc.per_helx(list_1sug_dis_9, False)
sug11_alpha_per, sug11_struct_per, sug11_alpha_mean, sug11_struct_mean, sug11_alpha_sem, sug11_struct_sem = mdfunc.per_helx(list_1sug_dis_11, False)
BBR_a7_alpha_per, BBR_a7_struct_per, BBR_a7_alpha_mean, BBR_a7_struct_mean, BBR_a7_alpha_sem, BBR_a7_struct_sem = mdfunc.per_helx(list_BBR_a7, False)
BBR_1sug_alpha_per, BBR_1sug_struct_per, BBR_1sug_alpha_mean, BBR_1sug_struct_mean, BBR_1sug_alpha_sem, BBR_1sug_struct_sem = mdfunc.per_helx(list_BBR_1sug, False)
AD_BBR_alpha_per, AD_BBR_struct_per, AD_BBR_alpha_mean, AD_BBR_struct_mean, AD_BBR_alpha_sem, AD_BBR_struct_sem = mdfunc.per_helx(list_AD_BBR, False)

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
ax1.set_title('Degree of Helicity', fontsize=15)    
ax1.set_ylabel('% Residue was in Alpha Helix', fontsize=13)
ax1.set_ylim(0,100)
ax1.plot(num, a7_alpha_mean, label='Apo Open', color = 'blue')
ax1.fill_between(num, a7_alpha_mean - a7_alpha_sem, a7_alpha_mean + a7_alpha_sem, alpha=0.2, facecolor = 'blue', edgecolor = 'blue')
ax1.plot(num, a7_AD_alpha_mean, label='AD', linestyle = 'dashed', color = 'blue')
ax1.fill_between(num, a7_AD_alpha_mean - a7_AD_alpha_sem, a7_AD_alpha_mean + a7_AD_alpha_sem, alpha=0.2, facecolor = 'blue', edgecolor = 'blue')
ax1.plot(num, BBR_a7_alpha_mean, label='BBR', linestyle = 'dotted', color = 'blue')
ax1.fill_between(num, BBR_a7_alpha_mean - BBR_a7_alpha_sem, BBR_a7_alpha_mean + BBR_a7_alpha_sem, alpha=0.2, facecolor = 'blue', edgecolor = 'blue')
ax1.plot(num, sug_alpha_mean, label='Apo Closed', color = 'red')
ax1.fill_between(num, sug_alpha_mean - sug_alpha_sem, sug_alpha_mean + sug_alpha_sem, alpha=0.3, facecolor = 'red', edgecolor = 'red')
ax1.set_xlabel('Residue ID', fontsize=13)
plt.xticks(fontsize=11)
plt.yticks(fontsize=11)
leg = ax1.legend(fontsize=11)
fig.savefig('Helicity_disordering_simp_err.png')
plt.close(fig)

fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.set_title('Degree of Helicity')    
ax1.set_ylabel('% Residue was in Alpha Helix')
ax1.set_ylim(0,100)
ax1.plot(num, a7_AD_alpha_mean, label='AD', linestyle = 'dashed', color = 'blue')
ax1.fill_between(num, a7_AD_alpha_mean - a7_AD_alpha_sem, a7_AD_alpha_mean + a7_AD_alpha_sem, alpha=0.3, facecolor = 'blue', edgecolor = 'blue')
ax1.plot(num, BBR_a7_alpha_mean, label='BBR', linestyle = 'dotted', color = 'purple')
ax1.fill_between(num, BBR_a7_alpha_mean - BBR_a7_alpha_sem, BBR_a7_alpha_mean + BBR_a7_alpha_sem, alpha=0.3, facecolor = 'blue', edgecolor = 'blue')
ax1.plot(num, AD_BBR_alpha_mean, label='AD BBR Combo', linestyle = 'dotted', color = 'green')
ax1.fill_between(num, AD_BBR_alpha_mean - AD_BBR_alpha_sem, AD_BBR_alpha_mean + AD_BBR_alpha_sem, alpha=0.3, facecolor = 'blue', edgecolor = 'blue')
ax1.set_xlabel('Residue ID')
leg = ax1.legend()
fig.savefig('Helicity_disordering_combo_simp_err.png')
plt.close(fig)

#Load DSSP for full trajectories
dssp_a7_all = open(file_path[22].strip(), 'r').readlines()
dssp_a7_AD_all = open(file_path[23].strip(), 'r').readlines()
dssp_1sug_all = open(file_path[24].strip(), 'r').readlines()
dssp_1sug_AD_all = open(file_path[25].strip(), 'r').readlines()
dssp_BBR_a7_all = open(file_path[26].strip(), 'r').readlines()
dssp_BBR_1sug_all = open(file_path[27].strip(), 'r').readlines()

#Detemine % alpha helicity at each time point
alpha_a7 = mdfunc.per_helx(dssp_a7_all, True) 
alpha_a7_AD = mdfunc.per_helx(dssp_a7_AD_all, True) 
alpha_1sug = mdfunc.per_helx(dssp_1sug_all, True)
alpha_1sug_AD = mdfunc.per_helx(dssp_1sug_AD_all, True) 
alpha_BBR_a7 = mdfunc.per_helx(dssp_BBR_a7_all, True)
alpha_BBR_1sug = mdfunc.per_helx(dssp_BBR_1sug_all, True)

#Compute running time averages
alpha_a7_avg = mdfunc.moving_average(alpha_a7, 25)
alpha_a7_AD_avg = mdfunc.moving_average(alpha_a7_AD, 25)
alpha_1sug_avg = mdfunc.moving_average(alpha_1sug, 25)
alpha_1sug_AD_avg = mdfunc.moving_average(alpha_1sug_AD, 25)
alpha_BBR_a7_avg = mdfunc.moving_average(alpha_BBR_a7, 25)
alpha_BBR_1sug_avg = mdfunc.moving_average(alpha_BBR_1sug, 25)

time_1sug = np.linspace(0, 300, num=len(alpha_1sug_avg))
time_a7 = np.linspace(0, 300, num=len(alpha_a7_avg))
time_1sug_AD = np.linspace(0, 300, num=len(alpha_1sug_AD_avg))
time_a7_AD = np.linspace(0, 300, num=len(alpha_a7_AD_avg))
time_BBR_a7 = np.linspace(0, 300, num=len(alpha_BBR_a7_avg))
time_BBR_1sug = np.linspace(0, 300, num=len(alpha_BBR_1sug_avg))

fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.set_title(r'Disordering of the $\alpha$7 Helix with Ligand Binding', fontsize = 15)
ax1.set_ylabel(r'% $\alpha$ Helicity', fontsize = 13)
ax1.set_ylim(0,100)
ax1.set_xlim(0,300)
ax1.plot(time_a7, alpha_a7_avg, label='Apo Open', color = 'gray')
ax1.plot(time_a7_AD, alpha_a7_AD_avg, label='AD', color = 'blue')
ax1.plot(time_BBR_a7, alpha_BBR_a7_avg, label='BBR', color = 'purple')
ax1.plot(time_1sug, alpha_1sug_avg, label='Apo Closed', color = 'red')
ax1.set_xlabel('Time(ns)', fontsize = 13)
plt.xticks(fontsize = 11)
plt.yticks(fontsize=11)
leg = ax1.legend(fontsize=12)
fig.savefig('Helicity_disordering_time.png')
plt.close(fig)

