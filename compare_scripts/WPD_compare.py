#Import Necessary Packages
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

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
        'rebuild_a7_high/config9/analysis', 'rebuild_a7_high/config11/analysis', 'Apo_dis/analysis', 'AD_dis/analysis/config7', 'AD_dis/analysis/config9', 'AD_dis/analysis/config11', 'BBR_a7/analysis', 'BBR_1sug/analysis']
all_name = ['a7', 'a7_AD', 'Apo', 'AD', '1sug', '1sug2', '1sug3', '1sug_AD', '1sug_no_a7', '1sug_no_a7_AD', 'config7', 'config9', 'config11', 'complex7', 'complex9', 'complex11', 'complex11_2', 'alt', 'alt2', 'config7', 'config9', 'config11', 'Apo_dis', 'complex7', 'complex9', 'complex11', 'BBR_a7', 'BBR_1sug']
eq_time_list = [50, 80, 50, 25, 5, 5, 5, 80, 75, 5, 10, 5, 5, 5, 5, 5, 5, 60, 30, 150, 5, 75, 5, 15, 25, 25, 70, 5]

#empty vector for the percent of time the WPD loop is open
per = np.zeros(len(all_dir))
per_eq = np.zeros(len(all_dir))
per_apo_no_a7_o, per_apo_order_o, per_apo_dis_o = [],[],[]
per_apo_no_a7_c, per_apo_order_c, per_apo_dis_c = [],[],[]
per_AD_no_a7, per_AD_order, per_AD_dis = [],[],[]
per_BBR_order, per_BBR_dis = [],[]

for i in range(len(all_dir)):
    data, per[i], per_eq[i] = load_data(all_dir[i], all_name[i], eq_time_list[i])

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
        per_AD_dis.append(per_eq[i])
    if i == 14:
        w_1sug_dis9_AD = data
        per_AD_dis.append(per_eq[i])
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
        w_AD_dis7 = data
        per_AD_dis.append(per_eq[i])
    if i == 23:
        w_AD_dis9 = data
        per_AD_dis.append(per_eq[i])
    if i == 24:
        w_AD_dis11 = data
        per_AD_dis.append(per_eq[i])
    if i == 25:
        w_BBR_a7 = data
        per_BBR_order.append(per_eq[i])
    if i == 26:
        w_BBR_1sug = data
        per_BBR_order.append(per_eq[i])
    if i != 0:
        WPD_all.extend(data) 

#Plot Percentages
num = [5, 10, 15]
Method = ['Absent', 'Ordered', 'Disordered']
#Apo Open WPD Loop dependance on a7 helix conformation
per = [np.mean(per_apo_no_a7_o), np.mean(per_apo_order_o), np.mean(per_apo_dis_o)]
per_err = [stats.sem(per_apo_no_a7_o), stats.sem(per_apo_order_o), stats.sem(per_apo_dis_o)]
fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.set_title('WPD Loop with Open Initial Conformation')
ax1.set_ylabel('% Time WPD loop Open ($D_{PtoWPD}$ $\geq$ 1nm)')
ax1.set_xlabel(r'Conformation of the $\alpha$-7 Helix')
ax1.bar(num, per, color = ['gray', 'blue', 'red'], width=4.5)
plt.xticks(num, Method, fontsize=12)
fig.savefig('WPD_open_per.png')

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
fig2 = plt.figure()
ax1 = fig2.add_subplot(111)
ax1.set_title("WPD Loop Opening with AD and Ordered a7")    
ax1.set_xlabel('Time (ns)')
ax1.set_ylabel('Residue Distances (nm)')
ax1.plot(time,w_1sug_AD, label='1sug_AD')
plt.axhline(y=1.0, color='r', linestyle='-')
fig2.savefig('1sug_opening.png')

#WPD Loop Distance for 1sug_dis_AD opening
time = np.linspace(0, 300, num = len(w_1sug_dis9_AD))
fig3 = plt.figure()
ax1 = fig3.add_subplot(111)
ax1.set_title("WPD Loop Opening with AD and Disordered a7")    
ax1.set_xlabel('Time (ns)')
ax1.set_ylabel('Residue Distances (nm)')
ax1.plot(time,w_1sug_dis9_AD, label='1sug_dis9_AD')
plt.axhline(y=1.0, color='r', linestyle='-')
fig3.savefig('1sug_dis_AD_opening.png')

#WPD Loop Distance for 1sug_dis opening
time = np.linspace(0, 300, num = len(w_1sug_dis11_AD))
fig4 = plt.figure()
ax1 = fig4.add_subplot(111)
ax1.set_title("WPD Loop Opening with AD and Disordered a7")    
ax1.set_xlabel('Time (ns)')
ax1.set_ylabel('Residue Distances (nm)')
ax1.plot(time,w_1sug_dis11_AD, label='1sug_dis11')
plt.axhline(y=1.0, color='r', linestyle='-')
fig4.savefig('1sug_dis_opening.png')

#Histograph of all WPD lengths
fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.set_title('Histogram of WPD Loop Distances')
ax1.set_xlabel('Distance (nm)')
n, bin1, patches = ax1.hist(WPD_all, bins = 30, color = 'blue', density=True)
for i in range(8):
    patches[i].set_fc('red')
fig.savefig('WPD_hist.png')

#Histograph of Apo vs AD open lengths
wpd_apo = w_a7
wpd_apo.extend(w_1sug_dis7)
wpd_apo.extend(w_dis9)

wpd_AD = w_1sug_dis11_AD
wpd_AD.extend(w_1sug_alt_AD)
wpd_AD.extend(w_1sug_alt2_AD)

fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.set_title('WPD Loop Distances Compared Between Apo and Ligand Bound')
ax1.set_xlabel('Distance (nm)', fontsize = 14)
ax1.hist(wpd_AD, bins = 25, alpha = 0.5, color = 'black', density = True, label = 'AD')
ax1.hist(wpd_apo, bins = 25, alpha = 0.5, color = 'blue', density = True, label = 'Apo')
ax1.legend(loc = 'upper right')
fig.savefig('WPD_AD_cmpr_hist.png')

