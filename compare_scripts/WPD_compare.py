#Import Necessary Packages
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator)

def load_data(file_dir, file_name, eq_time):
    t, w = [],[]
    #Equilibrium time 
    eq_time = eq_time * 1000
    #Load data
    with open('../../' + file_dir + '/' + file_name + '_WPD.xvg') as f:
        for _ in range(17):
            next(f)
        for line in f:
            cols = line.split()
            if len(cols) == 2:
                t.append(float(cols[0]))
                w.append(float(cols[1]))
    #Count variable for number of time steps the WPD loop is open
    count = 0
    count_eq = 0 #equilibrium section only
    count_eq_tot = 0 #total for equilibrium section

    for i in range(len(t)):
        if w[i]>1:
            count +=1
            if t[i] >= eq_time:
                count_eq += 1
        if t[i] >= eq_time:
            count_eq_tot += 1
    per_eq = 100*count_eq/count_eq_tot
    per = 100*count/len(t)
    
    return w, per, per_eq

#Input Data and Determine the Percent of Time the WPD loop is open
all_dir = ['rebuild_a7/analysis', 'AD_rebuild_a7/analysis', 'Apo/analysis', 'AD/analysis', 'Apo_1SUG/analysis/1sug', 'Apo_1SUG/analysis/1sug2', 'Apo_1SUG/analysis/1sug3', '1sug_AD/analysis', '1sug_no_a7/analysis', '1sug_no_a7_AD/analysis', '1sug_dis/analysis/config7', '1sug_dis/analysis/config9', 
        '1sug_dis/analysis/config11', '1sug_dis_AD/analysis/config7', '1sug_dis_AD/analysis/config9', '1sug_dis_AD/analysis/config11', '1sug_dis_AD/analysis/config11_2', '1sug_dis_AD/analysis/config_alt', '1sug_dis_AD/analysis/config_alt2', 'rebuild_a7_high/config7/analysis', 
        'rebuild_a7_high/config9/analysis', 'rebuild_a7_high/config11/analysis', 'Apo_dis/analysis', 'AD_dis/analysis/config7', 'AD_dis/analysis/config9', 'AD_dis/analysis/config11', 'BBR_a7/analysis', 'BBR_1sug/analysis', 'BBR_dis/analysis/config7', 'BBR_dis/analysis/config9', 'BBR_dis/analysis/config11', 
        'BBR_1sug_dis/analysis/config7', 'BBR_1sug_dis/analysis/config11']
all_name = ['a7', 'a7_AD', 'Apo', 'AD', '1sug', '1sug2', '1sug3', '1sug_AD', '1sug_no_a7', '1sug_no_a7_AD', 'config7', 'config9', 'config11', 'complex7', 'complex9', 'complex11', 'complex11_2', 'alt', 'alt2', 'config7', 'config9', 'config11', 'Apo_dis', 'complex7', 'complex9', 'complex11', 'BBR_a7', 
        'BBR_1sug', 'BBR_dis7', 'BBR_dis9', 'BBR_dis11', 'BBR_dis7', 'BBR_dis11']
eq_time_list = [50, 80, 50, 25, 5, 5, 5, 80, 75, 5, 10, 5, 5, 5, 5, 5, 5, 60, 30, 150, 5, 75, 5, 15, 25, 25, 70, 5, 70, 5, 75, 60, 60, 10, 10]

#empty vector for the percent of time the WPD loop is open
per = np.zeros(len(all_dir))
per_eq = np.zeros(len(all_dir))
per_apo_no_a7_o, per_apo_order_o, per_apo_dis_o = [],[],[]
per_apo_no_a7_c, per_apo_order_c, per_apo_dis_c = [],[],[]
per_AD_no_a7, per_AD_order, per_AD_dis = [],[],[]
per_BBR_order, per_BBR_dis = [],[]

file_per = open('Per_all.txt', 'w') #File for all percentages

for i in range(len(all_dir)):
    data, per[i], per_eq[i] = load_data(all_dir[i], all_name[i], eq_time_list[i])
    
    #Save all percentages to files
    file_per.write(all_name[i] + ':\n' + 'Full: ' + str(per[i]) + '\n' + 'Equilibrated: ' + str(per_eq[i]) + '\n')
    
    if i == 0:
        w_a7 = data
        WPD_all = data
        per_apo_order_o.append(per_eq[i])
    if i == 1:
        w_a7_AD = data
        per_AD_order.append(per_eq[i])
    if i == 2:
        w_Apo = data
        per_apo_no_a7_o.append(per_eq[i])
    if i == 3:
        w_AD = data
        per_AD_no_a7.append(per_eq[i])
    if i == 4:
        w_1sug = data
        per_apo_order_c.append(per_eq[i])
    if i == 5:
        w_1sug2 = data
        per_apo_order_c.append(per_eq[i])
    if i == 6:
        w_1sug3 = data
        per_apo_order_c.append(per_eq[i])
    if i == 7:
        w_1sug_AD = data
        per_AD_order.append(per_eq[i])
    if i == 8:
        w_1sug_na7 = data
        per_apo_no_a7_c.append(per_eq[i])
    if i == 9:
        w_1sug_na7_AD = data
        per_AD_no_a7.append(per_eq[i])
    if i == 10:
        w_1sug_dis7 = data
        per_apo_dis_c.append(per_eq[i])
    if i == 11:
        w_1sug_dis9 = data
        per_apo_dis_c.append(per_eq[i])
    if i == 12:
        w_1sug_dis11 = data
        per_apo_dis_c.append(per_eq[i])
    if i == 13:
        w_1sug_dis7_AD = data
    if i == 14:
        w_1sug_dis9_AD = data
    if i == 15:
        w_1sug_dis11_AD = data
        per_AD_dis.append(per_eq[i])
    if i == 16:
        w_1sug_dis11_2_AD = data
        per_AD_dis.append(per_eq[i])
    if i == 17:
        w_1sug_alt_AD = data
        per_AD_dis.append(per_eq[i])
    if i == 18:
        w_1sug_alt2_AD = data
        per_AD_dis.append(per_eq[i])
    if i == 19:
        w_dis7 = data
        per_apo_dis_o.append(per_eq[i])
    if i == 20:
        w_dis9 = data
        per_apo_dis_o.append(per_eq[i])
    if i == 21:
        w_dis11 = data
        per_apo_dis_o.append(per_eq[i])
    if i == 22:
        w_apo_dis = data
        per_apo_dis_o.append(per_eq[i])
    if i == 23:
        w_AD_dis7 = data
    if i == 24:
        w_AD_dis9 = data
    if i == 25:
        w_AD_dis11 = data
        per_AD_dis.append(per_eq[i])
    if i == 26:
        w_BBR_a7 = data
        per_BBR_order.append(per_eq[i])
    if i == 27:
        w_BBR_1sug = data
    if i == 28:
        w_BBR_dis7 = data
        per_BBR_dis.append(per_eq[i])
    if i == 29:
        w_BBR_dis9 = data
        per_BBR_dis.append(per_eq[i])
    if i == 30:
        w_BBR_dis11 = data
    if i == 31:
        w_BBR_1sug_dis7 = data
    if i == 32:
        w_BBR_1sug_dis11 = data
        per_BBR_dis.append(per_eq[i])

    if i != 0:
        WPD_all.extend(data) 
file_per.close()

#Plot Percentages
num = [1, 9, 17]
num2 = [3, 11, 19]
num3 = [5, 13, 21]
num4 = [7, 15, 23]
Method = ['Absent', 'Ordered', 'Disordered']
#Apo Open WPD Loop dependance on a7 helix conformation
per_ao = [np.mean(per_apo_no_a7_o), np.mean(per_apo_order_o), np.mean(per_apo_dis_o)]
per_ac = [np.mean(per_apo_no_a7_c), np.mean(per_apo_order_c), np.mean(per_apo_dis_c)]
per_AD = [np.mean(per_AD_no_a7), np.mean(per_AD_order), np.mean(per_AD_dis)]
per_BBR = [0 , np.mean(per_BBR_order), np.mean(per_BBR_dis)]

per_err_ao = [stats.sem(per_apo_no_a7_o), stats.sem(per_apo_order_o), stats.sem(per_apo_dis_o)]
per_err_ac = [stats.sem(per_apo_no_a7_c), stats.sem(per_apo_order_c), stats.sem(per_apo_dis_c)]
per_err_AD = [stats.sem(per_AD_no_a7), stats.sem(per_AD_order), stats.sem(per_AD_dis)]
per_err_BBR = [0 , stats.sem(per_BBR_order), stats.sem(per_BBR_dis)]

fig = plt.figure(figsize = (7,7))
ax1 = fig.add_subplot(111)
ax1.set_title('Comparison of WPD Loop Conformation', fontsize = 18)
ax1.set_ylabel('% Time WPD loop Open ($D_{PtoWPD}$ $\geq$ 1nm)', fontsize = 14)
ax1.set_xlabel(r'Conformation of the $\alpha$-7 Helix', fontsize = 14)
ax1.bar(num, per_ao, yerr = per_err_ao, color = 'gray', width=1.9, label = 'Apo Open')
ax1.bar(num2, per_ac, yerr = per_err_ac, color = 'red', width=1.9, label = 'Apo Closed')
ax1.bar(num3, per_AD, yerr = per_err_AD, color = 'blue', width=1.9, label = 'AD')
ax1.bar(num4, per_BBR, yerr = per_err_BBR, color = 'purple', width=1.9, label = 'BBR')
plt.xticks(num2, Method, fontsize = 14)
plt.yticks(fontsize = 12)
plt.legend(loc='lower right', fontsize = 12)
fig.savefig('WPD_per_all.png')

#Apo Closed WPD Loop dependance on a7 helix conformation
per = [np.mean(per_apo_no_a7_c), np.mean(per_apo_order_c), np.mean(per_apo_dis_c)]
per_err = [stats.sem(per_apo_no_a7_c), stats.sem(per_apo_order_c), stats.sem(per_apo_dis_c)]
fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.set_title('WPD Loop with Closed Initial Conformation')
ax1.set_ylabel('% Time WPD loop Open ($D_{PtoWPD}$ $\geq$ 1nm)')
ax1.set_xlabel(r'Conformation of the $\alpha$-7 Helix')
ax1.bar(num, per, color = ['gray', 'blue', 'red'], width=4.5)
plt.xticks(num, Method, fontsize=12)
fig.savefig('WPD_closed_per.png')

#WPD Loop Distance for 1sug_AD opening
time = np.linspace(0, 300, num = len(w_1sug_AD))
time2 = np.linspace(0, 300, num = len(w_BBR_a7))
fig2 = plt.figure()
ax1 = fig2.add_subplot(111)
ax1.set_title("WPD Loop Opening with Ordered a7", fontsize = 18)
ax1.set_xlabel('Time (ns)', fontsize = 14)
ax1.set_ylabel('Residue Distances (nm)', fontsize = 14)
ax1.plot(time,w_1sug_AD, label='AD', color = 'blue')
ax1.plot(time2, w_BBR_a7, label='BBR', color = 'purple')
plt.xticks(fontsize = 12)
plt.yticks(fontsize = 12)
ax1.tick_params(axis='x', which = 'major', length = 7, color='black', width = 2)
ax1.tick_params(axis='x', which = 'minor', length = 4, color='black', width = 1)
ax1.tick_params(axis='y', which = 'major', length = 7, color='black', width = 2)
ax1.tick_params(axis='y', which = 'minor', length = 4, color='black', width = 1)
plt.legend(loc = 'best')
plt.axhline(y=1.0, color='r', linestyle='-')
fig2.savefig('order_opening.png')

#WPD Loop Distance for 1sug_dis opening
time = np.linspace(0, 300, num = len(w_1sug_dis11_AD))
time2 = np.linspace(0, 300, num = len(w_BBR_1sug_dis11))
fig4 = plt.figure()
ax1 = fig4.add_subplot(111)
ax1.set_title("WPD Loop Opening with Disordered a7", fontsize = 18)    
ax1.set_xlabel('Time (ns)', fontsize = 14)
ax1.set_ylabel('Residue Distances (nm)', fontsize = 14)
ax1.plot(time,w_1sug_dis11_AD, label='AD', color = 'blue')
ax1.plot(time2,w_BBR_1sug_dis11, label='BBR', color = 'purple')
plt.xticks(fontsize = 12)
plt.yticks(fontsize = 12)
ax1.tick_params(axis='x', which = 'major', length = 7, color='black', width = 2)
ax1.tick_params(axis='x', which = 'minor', length = 4, color='black', width = 1)
ax1.tick_params(axis='y', which = 'major', length = 7, color='black', width = 2)
ax1.tick_params(axis='y', which = 'minor', length = 4, color='black', width = 1)
plt.axhline(y=1.0, color='r', linestyle='-')
plt.legend(loc='best')
fig4.savefig('dis_opening.png')

#Histograph of all WPD lengths
fig,ax1 = plt.subplots()

ax1.xaxis.set_minor_locator(AutoMinorLocator())
ax1.yaxis.set_minor_locator(AutoMinorLocator())

ax1.set_title('Histogram of WPD Loop Distances', fontsize = 18)
ax1.set_xlabel('Distance (nm)', fontsize = 14)
ax1.set_ylabel('Occurance', fontsize = 14)
n, bin1, patches = ax1.hist(WPD_all, bins = 35, color = 'gray', density=True)
for i in range(9):
    patches[i].set_fc('red')
plt.xticks(fontsize = 12)
plt.yticks(fontsize = 12)
ax1.tick_params(axis='x', which = 'major', length = 7, color='black', width = 2)
ax1.tick_params(axis='x', which = 'minor', length = 4, color='black', width = 1)
ax1.tick_params(axis='y', which = 'major', length = 7, color='black', width = 2)
ax1.tick_params(axis='y', which = 'minor', length = 4, color='black', width = 1)
fig.savefig('WPD_hist.png')

#Histograph of Apo vs AD open lengths
wpd_apo = w_a7
wpd_apo.extend(w_1sug_dis7)
wpd_apo.extend(w_dis9)

wpd_AD = w_1sug_dis11_AD
wpd_AD.extend(w_1sug_alt_AD)
wpd_AD.extend(w_1sug_alt2_AD)

wpd_BBR = w_BBR_a7
wpd_BBR.extend(w_BBR_dis7)
wpd_BBR.extend(w_BBR_dis9)
wpd_BBR.extend(w_BBR_1sug_dis11)

fig = plt.figure(figsize=(8,6))
ax1 = fig.add_subplot(111)
ax1.set_title('WPD Loop Distances: Apo vs Ligand Bound', fontsize = 18)
ax1.set_xlabel('Distance (nm)', fontsize = 14)
ax1.set_ylabel('Occurance', fontsize = 14)
ax1.hist(wpd_BBR, bins = 25, alpha = 0.3, color = 'purple', density = True, label = 'BBR')
ax1.hist(wpd_AD, bins = 25, alpha = 0.3, color = 'blue', density = True, label = 'AD')
ax1.hist(wpd_apo, bins = 25, alpha = 0.5, color = 'gray', density = True, label = 'Apo')
ax1.legend(loc = 'upper right', fontsize = 18)
plt.xticks(fontsize = 12)
plt.yticks(fontsize = 12)
ax1.tick_params(axis='x', which = 'major', length = 7, color='black', width = 2)
ax1.tick_params(axis='x', which = 'minor', length = 4, color='black', width = 1)
ax1.tick_params(axis='y', which = 'major', length = 7, color='black', width = 2)
ax1.tick_params(axis='y', which = 'minor', length = 4, color='black', width = 1)
fig.savefig('WPD_AD_cmpr_hist.png')

