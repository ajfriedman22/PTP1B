#Import Necessary Packages
import numpy as np
import matplotlib.pyplot as plt

#Make open arrays for time and atomic distances
t_a7, t_a7_AD, t_Apo, t_AD, t_1sug, t_1sug_AD, t_1sug_na7, t_1sug_na7_AD, t_1sug_dis7, t_1sug_dis7_AD, t_1sug_dis9, t_1sug_dis9_AD, t_1sug_dis11, t_1sug_dis11_AD, t_Apo_dis7, t_AD_dis7, t_Apo_dis9, t_AD_dis9, t_Apo_dis11, t_AD_dis11, t_1sug_AD_dis11_2 = [],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[] #time
w_a7, w_a7_AD, w_Apo, w_AD, w_1sug, w_1sug_AD, w_1sug_na7, w_1sug_na7_AD, w_1sug_dis7, w_1sug_dis7_AD, w_1sug_dis9, w_1sug_dis9_AD, w_1sug_dis11, w_1sug_dis11_AD, w_Apo_dis7, w_AD_dis7, w_Apo_dis9, w_AD_dis9, w_Apo_dis11, w_AD_dis11, w_1sug_AD_dis11_2 = [],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[] #atom distance
t_1sug_alt, t_1sug_alt2 = [],[]
w_1sug_alt, w_1sug_alt2 = [],[]

w_a7_eq, w_a7_AD_eq, w_Apo_eq, w_AD_eq, w_1sug_eq, w_1sug_AD_eq, w_1sug_na7_eq, w_1sug_na7_AD_eq, w_1sug_dis7_eq, w_1sug_dis7_AD_eq, w_1sug_dis9_eq, w_1sug_dis9_AD_eq, w_1sug_dis11_eq, w_1sug_dis11_AD_eq, w_Apo_dis7_eq, w_AD_dis7_eq, w_Apo_dis9_eq, w_AD_dis9_eq, w_Apo_dis11_eq, w_AD_dis11_eq, w_1sug_AD_dis11_2_eq = [],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[] #atom distance
w_1sug_alt_eq, w_1sug_alt2_eq = [],[]

#Count variable for number of time steps the WPD loop is open
count = 0
count_eq = 0 #equilibrium section only
count_eq_tot = 0 #total for equilibrium section

#empty vector for the percent of time the WPD loop is open
per, per_eq = [],[]

#Input Data
with open("../../rebuild_a7/a7_WPD.xvg") as f:
    for _ in range(17):
        next(f)
    for line in f:
        cols = line.split()
        if len(cols) == 2:
            t_a7.append(float(cols[0]))
            w_a7.append(float(cols[1]))
with open("../../AD_rebuild_a7/a7_AD_WPD.xvg") as f:
    for _ in range(17):
        next(f)
    for line in f:
        cols = line.split()
        if len(cols) == 2:
            t_a7_AD.append(float(cols[0]))
            w_a7_AD.append(float(cols[1]))
with open("../../Apo/Apo_WPD.xvg") as f:
    for _ in range(17):
        next(f)
    for line in f:
        cols = line.split()
        if len(cols) == 2:
            t_Apo.append(float(cols[0]))
            w_Apo.append(float(cols[1]))
with open("../../AD/AD_WPD.xvg") as f:
    for _ in range(17):
        next(f)
    for line in f:
        cols = line.split()
        if len(cols) == 2:
            t_AD.append(float(cols[0]))
            w_AD.append(float(cols[1]))
with open("../../1sug/1sug_WPD.xvg") as f:
    for _ in range(17):
        next(f)
    for line in f:
        cols = line.split()
        if len(cols) == 2:
            t_1sug.append(float(cols[0]))
            w_1sug.append(float(cols[1]))
with open("../../1sug_AD/1sug_AD_WPD.xvg") as f:
    for _ in range(17):
        next(f)
    for line in f:
        cols = line.split()
        if len(cols) == 2:
            t_1sug_AD.append(float(cols[0]))
            w_1sug_AD.append(float(cols[1]))
with open("../../1sug_no_a7/1sug_no_a7_WPD.xvg") as f:
    for _ in range(17):
        next(f)
    for line in f:
        cols = line.split()
        if len(cols) == 2:
            t_1sug_na7.append(float(cols[0]))
            w_1sug_na7.append(float(cols[1]))
with open("../../1sug_no_a7_AD/1sug_no_a7_AD_WPD.xvg") as f:
    for _ in range(17):
        next(f)
    for line in f:
        cols = line.split()
        if len(cols) == 2:
            t_1sug_na7_AD.append(float(cols[0]))
            w_1sug_na7_AD.append(float(cols[1]))
with open("../../1sug_dis/config7/config7_WPD.xvg") as f:
    for _ in range(17):
        next(f)
    for line in f:
        cols = line.split()
        if len(cols) == 2:
            t_1sug_dis7.append(float(cols[0]))
            w_1sug_dis7.append(float(cols[1]))
with open("../../1sug_dis/config9/config9_WPD.xvg") as f:
    for _ in range(17):
        next(f)
    for line in f:
        cols = line.split()
        if len(cols) == 2:
            t_1sug_dis9.append(float(cols[0]))
            w_1sug_dis9.append(float(cols[1]))
with open("../../1sug_dis/config11/config11_WPD.xvg") as f:
    for _ in range(17):
        next(f)
    for line in f:
        cols = line.split()
        if len(cols) == 2:
            t_1sug_dis11.append(float(cols[0]))
            w_1sug_dis11.append(float(cols[1]))
with open("../../1sug_dis_AD/config7/complex7_WPD.xvg") as f:
    for _ in range(17):
        next(f)
    for line in f:
        cols = line.split()
        if len(cols) == 2:
            t_1sug_dis7_AD.append(float(cols[0]))
            w_1sug_dis7_AD.append(float(cols[1]))
with open("../../1sug_dis_AD/config9/complex9_WPD.xvg") as f:
    for _ in range(17):
        next(f)
    for line in f:
        cols = line.split()
        if len(cols) == 2:
            t_1sug_dis9_AD.append(float(cols[0]))
            w_1sug_dis9_AD.append(float(cols[1]))
with open("../../1sug_dis_AD/config11/complex11_WPD.xvg") as f:
    for _ in range(17):
        next(f)
    for line in f:
        cols = line.split()
        if len(cols) == 2:
            t_1sug_dis11_AD.append(float(cols[0]))
            w_1sug_dis11_AD.append(float(cols[1]))
with open("../../Apo_dis/config7/config7_WPD.xvg") as f:
    for _ in range(17):
        next(f)
    for line in f:
        cols = line.split()
        if len(cols) == 2:
            t_Apo_dis7.append(float(cols[0]))
            w_Apo_dis7.append(float(cols[1]))
with open("../../AD_dis/config7/complex7_WPD.xvg") as f:
    for _ in range(17):
        next(f)
    for line in f:
        cols = line.split()
        if len(cols) == 2:
            t_AD_dis7.append(float(cols[0]))
            w_AD_dis7.append(float(cols[1]))
with open("../../Apo_dis/config9/config9_WPD.xvg") as f:
    for _ in range(17):
        next(f)
    for line in f:
        cols = line.split()
        if len(cols) == 2:
            t_Apo_dis9.append(float(cols[0]))
            w_Apo_dis9.append(float(cols[1]))
with open("../../AD_dis/config9/complex9_WPD.xvg") as f:
    for _ in range(17):
        next(f)
    for line in f:
        cols = line.split()
        if len(cols) == 2:
            t_AD_dis9.append(float(cols[0]))
            w_AD_dis9.append(float(cols[1]))
with open("../../Apo_dis/config11/config11_WPD.xvg") as f:
    for _ in range(17):
        next(f)
    for line in f:
        cols = line.split()
        if len(cols) == 2:
            t_Apo_dis11.append(float(cols[0]))
            w_Apo_dis11.append(float(cols[1]))
with open("../../AD_dis/config11/complex11_WPD.xvg") as f:
    for _ in range(17):
        next(f)
    for line in f:
        cols = line.split()
        if len(cols) == 2:
            t_AD_dis11.append(float(cols[0]))
            w_AD_dis11.append(float(cols[1]))
with open("../../1sug_dis_AD/config11_2/complex11_2_WPD.xvg") as f:
    for _ in range(17):
        next(f)
    for line in f:
        cols = line.split()
        if len(cols) == 2:
            t_1sug_AD_dis11_2.append(float(cols[0]))
            w_1sug_AD_dis11_2.append(float(cols[1]))
with open("../../1sug_AD_dis_alt/run_1/alt_WPD.xvg") as f:
    for _ in range(17):
        next(f)
    for line in f:
        cols = line.split()
        if len(cols) == 2:
            t_1sug_alt.append(float(cols[0]))
            w_1sug_alt.append(float(cols[1]))
with open("../../1sug_AD_dis_alt/run_2/alt2_WPD.xvg") as f:
    for _ in range(17):
        next(f)
    for line in f:
        cols = line.split()
        if len(cols) == 2:
            t_1sug_alt2.append(float(cols[0]))
            w_1sug_alt2.append(float(cols[1]))

#Determine the Percent of Time the WPD loop is open
for i in range(len(t_a7)):
    if w_a7[i]>1:
        count +=1
        if t_a7[i] >= 50000:
            count_eq += 1
    if t_a7[i] >= 50000:
        count_eq_tot += 1
        w_a7_eq.append(w_a7[i])
per_eq.append(100*count_eq/count_eq_tot)
per.append(100*count/len(t_a7))
count=0
count_eq=0
count_eq_tot=0
for i in range(len(t_a7_AD)):
    if w_a7_AD[i]>1:
        count +=1
        if t_a7_AD[i] >= 80000:
            count_eq += 1
    if t_a7_AD[i] >= 80000:
        count_eq_tot += 1
        w_a7_AD_eq.append(w_a7_AD[i])
per_eq.append(100*count_eq/count_eq_tot)
per.append(100*count/len(t_a7_AD))
count=0
count_eq=0
count_eq_tot=0
for i in range(len(t_Apo)):
    if w_Apo[i]>1:
        count +=1
        if t_Apo[i] >= 50000:
            count_eq += 1
    if t_Apo[i] >= 50000:
        count_eq_tot += 1
        w_Apo_eq.append(w_Apo[i])
per_eq.append(100*count_eq/count_eq_tot)
per.append(100*count/len(t_Apo))
count=0
count_eq=0
count_eq_tot=0
for i in range(len(t_AD)):
    if w_AD[i]>1:
        count +=1
        if t_AD[i] >= 25000:
            count_eq += 1
    if t_AD[i] >= 25000:
        count_eq_tot += 1
        w_AD_eq.append(w_AD[i])
per_eq.append(100*count_eq/count_eq_tot)
per.append(100*count/len(t_AD))
count=0
count_eq=0
count_eq_tot=0
for i in range(len(t_1sug)):
    if w_1sug[i]>1:
        count +=1
        if t_1sug[i] >= 5000:
            count_eq += 1
    if t_1sug[i] >= 5000:
        count_eq_tot += 1
        w_1sug_eq.append(w_1sug[i])
per_eq.append(100*count_eq/count_eq_tot)
per.append(100*count/len(t_1sug))
count=0
count_eq=0
count_eq_tot=0
for i in range(len(t_1sug_AD)):
    if w_1sug_AD[i]>1:
        count +=1
        if t_1sug_AD[i] >= 80000:
            count_eq += 1
    if t_1sug_AD[i] >= 80000:
        count_eq_tot += 1
        w_1sug_AD_eq.append(w_1sug_AD[i])

per_eq.append(100*count_eq/count_eq_tot)
per.append(100*count/len(t_1sug_AD))
count=0
count_eq=0
count_eq_tot=0
for i in range(len(t_1sug_na7)):
    if w_1sug_na7[i]>1:
        count +=1
        if t_1sug_na7[i] >= 75000:
            count_eq += 1
    if t_1sug_na7[i] >= 75000:
        count_eq_tot += 1
per_eq.append(100*count_eq/count_eq_tot)
per.append(100*count/len(t_1sug_na7))
count=0
count_eq=0
count_eq_tot=0
for i in range(len(t_1sug_na7_AD)):
    if w_1sug_na7_AD[i]>1:
        count +=1
        if t_1sug_na7_AD[i] >= 5000:
            count_eq += 1
    if t_1sug_na7_AD[i] >= 5000:
        count_eq_tot += 1
per_eq.append(100*count_eq/count_eq_tot)
per.append(100*count/len(t_1sug_na7_AD))
count=0
count_eq=0
count_eq_tot=0
for i in range(len(t_1sug_dis7)):
    if w_1sug_dis7[i]>1:
        count +=1
        if t_1sug_dis7[i] >= 5000:
            count_eq += 1
    if t_1sug_dis7[i] >= 5000:
        count_eq_tot += 1
        w_1sug_dis7_eq.append(w_1sug_dis7[i])
per_eq.append(100*count_eq/count_eq_tot)
per.append(100*count/len(t_1sug_dis7))
count=0
count_eq=0
count_eq_tot=0
for i in range(len(t_1sug_dis7_AD)):
    if w_1sug_dis7_AD[i]>1:
        count +=1
        if t_1sug_dis7_AD[i] >= 5000:
            count_eq += 1
    if t_1sug_dis7_AD[i] >= 5000:
        count_eq_tot += 1
per_eq.append(100*count_eq/count_eq_tot)
per.append(100*count/len(t_1sug_dis7_AD))
count=0
count_eq=0
count_eq_tot=0
for i in range(len(t_1sug_dis9)):
    if w_1sug_dis9[i]>1:
        count +=1
        if t_1sug_dis9[i] >= 5000:
            count_eq += 1
    if t_1sug_dis9[i] >= 5000:
        count_eq_tot += 1
per_eq.append(100*count_eq/count_eq_tot)
per.append(100*count/len(t_1sug_dis9))
count=0
count_eq=0
count_eq_tot=0
for i in range(len(t_1sug_dis9_AD)):
    if w_1sug_dis9_AD[i]>1:
        count +=1
        if t_1sug_dis9_AD[i] >= 5000:
            count_eq += 1
    if t_1sug_dis9_AD[i] >= 5000:
        count_eq_tot += 1
per_eq.append(100*count_eq/count_eq_tot)
per.append(100*count/len(t_1sug_dis9_AD))
count=0
count_eq=0
count_eq_tot=0
for i in range(len(t_1sug_dis11)):
    if w_1sug_dis11[i]>1:
        count +=1
        if t_1sug_dis11[i] >= 5000:
            count_eq += 1
    if t_1sug_dis11[i] >= 5000:
        count_eq_tot += 1
per_eq.append(100*count_eq/count_eq_tot)
per.append(100*count/len(t_1sug_dis11))
count=0
count_eq=0
count_eq_tot=0
for i in range(len(t_1sug_dis11_AD)):
    if w_1sug_dis11_AD[i]>1:
        count +=1
        if t_1sug_dis11_AD[i] >= 5000:
            count_eq += 1
    if t_1sug_dis11_AD[i] >= 5000:
        count_eq_tot += 1
        w_1sug_dis11_AD_eq.append(w_1sug_dis11_AD[i])
per_eq.append(100*count_eq/count_eq_tot)
per.append(100*count/len(t_1sug_dis11_AD))
count =0
count_eq=0
count_eq_tot=0
for i in range(len(t_Apo_dis7)):
    if w_Apo_dis7[i]>1:
        count +=1
        if t_Apo_dis7[i] >= 5000:
            count_eq += 1
    if t_Apo_dis7[i] >= 5000:
        count_eq_tot += 1
per_eq.append(100*count_eq/count_eq_tot)
per.append(100*count/len(t_Apo_dis7))
count =0
count_eq=0
count_eq_tot=0
for i in range(len(t_AD_dis7)):
    if w_AD_dis7[i]>1:
        count +=1
        if t_AD_dis7[i] >= 15000:
            count_eq += 1
    if t_AD_dis7[i] >= 15000:
        count_eq_tot += 1
per_eq.append(100*count_eq/count_eq_tot)
per.append(100*count/len(t_AD_dis7))
count =0
count_eq=0
count_eq_tot=0
for i in range(len(t_Apo_dis9)):
    if w_Apo_dis9[i]>1:
        count +=1
        if t_Apo_dis9[i] >= 5000:
            count_eq += 1
    if t_Apo_dis9[i] >= 5000:
        count_eq_tot += 1
        w_Apo_dis9_eq.append(w_Apo_dis9[i])
per_eq.append(100*count_eq/count_eq_tot)
per.append(100*count/len(t_Apo_dis9))
count =0
count_eq=0
count_eq_tot=0
for i in range(len(t_AD_dis9)):
    if w_AD_dis9[i]>1:
        count +=1
        if t_AD_dis9[i] >= 25000:
            count_eq += 1
    if t_AD_dis9[i] >= 25000:
        count_eq_tot += 1
per_eq.append(100*count_eq/count_eq_tot)
per.append(100*count/len(t_AD_dis9))
count =0
count_eq=0
count_eq_tot=0
for i in range(len(t_Apo_dis11)):
    if w_Apo_dis11[i]>1:
        count +=1
        if t_Apo_dis11[i] >= 75000:
            count_eq += 1
    if t_Apo_dis11[i] >= 75000:
        count_eq_tot += 1
per_eq.append(100*count_eq/count_eq_tot)
per.append(100*count/len(t_Apo_dis11))
count =0
count_eq=0
count_eq_tot=0
for i in range(len(t_AD_dis11)):
    if w_AD_dis11[i]>1:
        count +=1
        if t_AD_dis11[i] >= 25000:
            count_eq += 1
    if t_AD_dis11[i] >= 25000:
        count_eq_tot += 1
per_eq.append(100*count_eq/count_eq_tot)
per.append(100*count/len(t_AD_dis11))
count =0
count_eq=0
count_eq_tot=0
for i in range(len(t_1sug_AD_dis11_2)):
    if w_1sug_AD_dis11_2[i]>1:
        count +=1
        if t_1sug_AD_dis11_2[i] >= 5000:
            count_eq += 1
    if t_1sug_AD_dis11_2[i] >= 5000:
        count_eq_tot += 1
per_eq.append(100*count_eq/count_eq_tot)
per.append(100*count/len(t_1sug_AD_dis11_2))
count =0
count_eq=0
count_eq_tot=0
for i in range(len(t_1sug_alt)):
    if w_1sug_alt[i]>1:
        count +=1
        if t_1sug_alt[i] >= 60000:
            count_eq += 1
    if t_1sug_alt[i] >= 60000:
        count_eq_tot += 1
        w_1sug_alt_eq.append(w_1sug_alt[i])
per_eq.append(100*count_eq/count_eq_tot)
per.append(100*count/len(t_1sug_alt))
count =0
count_eq=0
count_eq_tot=0
for i in range(len(t_1sug_alt2)):
    if w_1sug_alt2[i]>1:
        count +=1
        if t_1sug_alt2[i] >= 30000:
            count_eq += 1
    if t_1sug_alt2[i] >= 30000:
        count_eq_tot += 1
        w_1sug_alt2_eq.append(w_1sug_alt2[i])
per_eq.append(100*count_eq/count_eq_tot)
per.append(100*count/len(t_1sug_alt2))

#Ordered a7
num = [5, 10, 15, 20]
Method = ['Apo Open\n Start', 'AD Open\n Start', 'Apo Closed\n Start', 'AD Closed\n Start']
fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.set_title("Time Spent Open with Ordered Helix")    
ax1.set_ylabel('% Greater than 1nm')
ax1.bar(num, [per[0], per[1], per[4], per[5]], color = ['darkblue', 'blue', 'darkred', 'red'], width=4.5)
plt.xticks(num, Method, fontsize=12)
fig.savefig('WPD_a7_order_per.png')

#a7 Absent
fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.set_title("Time Spent Open with No a7 Helix")    
ax1.set_ylabel('% Greater than 1nm')
ax1.bar(num, [per[2], per[3], per[6], per[7]], color = ['darkblue', 'blue', 'darkred', 'red'], width=4.5)
plt.xticks(num, Method, fontsize=12)
fig.savefig('WPD_na7_order_per.png')

#Disordered a7 config7
fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.set_title("Time Spent Open with Disordered Config7")    
ax1.set_ylabel('% Greater than 1nm')
ax1.bar(num, [per[14], per[15], per[8], per[9]], color = ['darkblue', 'blue', 'darkred', 'red'], width=4.5)
plt.xticks(num, Method, fontsize=12)
fig.savefig('WPD_a7_dis7_per.png')

#disordered a7 config9
fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.set_title("Time Spent Open with Disordered Config9")
ax1.set_ylabel('% Greater than 1nm')
ax1.bar(num, [per[16], per[17], per[10], per[11]], color = ['darkblue', 'blue', 'darkred', 'red'], width=4.5)
plt.xticks(num, Method, fontsize=12)
fig.savefig('WPD_a7_dis9_per.png')

#Disordered a7 config11
num = [5, 10, 15, 20, 25]
Method = ['Apo Open\n Start', 'AD Open\n Start', 'Apo Closed\n Start', 'AD Closed\n Start', 'AD Closed\n Start']
fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.set_title("Time Spent Open with Disordered Config11")    
ax1.set_ylabel('% Greater than 1nm')
ax1.bar(num, [per[18], per[19], per[12], per[13], per[20]], color = ['darkblue', 'blue', 'darkred', 'red'], width=4.5)
plt.xticks(num, Method, fontsize=12)
fig.savefig('WPD_a7_dis11_per.png')

#Disordered a7 config11
num = [5, 10]
Method = ['AD Closed\n Start', 'AD Closed\n Start']
fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.set_title("Time Spent Open with Disordered Config from BBR")    
ax1.set_ylabel('% Greater than 1nm')
ax1.bar(num, [per[21], per[22]], color = ['blue', 'blue'], width=4.5)
plt.xticks(num, Method, fontsize=12)
fig.savefig('WPD_a7_dis_alt_per.png')

#Plot Bar Graph for the percent of time WPD loop is open
#num = [5,10,15,20,25,30,35,40,45,50,55,60,65,70]
#Method = ['Open Start\n Ordered a7\n Apo', 'Open Start\n Ordered a7\n AD', 'Open Start\n Absent a7\n Apo', 'Open Start\n Absent a7\n AD', 'Close Start\n Ordered a7\n Apo', 'Close Start\n Ordered a7\n AD', 'Close Start\n Absent a7\n Apo', 'Close Start\n Absent a7\n AD', 'Close Start\n Disordered a7\n Apo', 'Close Start\n Disordered a7\n AD', 'Open Start\n Disordered a7\n Apo', 'Open Start\n Disordered a7\n AD', 'Open Start\n Disordered a7\n Apo', 'Open Start\n Disordered a7\n AD']
#fig = plt.figure(figsize=(25,15))
#ax1 = fig.add_subplot(111)
#ax1.set_title("Time Spent Open")    
#ax1.set_ylabel('% Greater than 1nm')
#ax1.bar(num,per, color = ['darkblue', 'blue', 'darkgreen', 'green', 'darkred', 'red', 'darkorange', 'orange', 'purple', 'magenta', 'brown', 'firebrick', 'gold', 'yellow'], width=4.5)
#plt.xticks(num, Method, fontsize=8)
#fig.savefig('WPD_percent.png')

#Plot Bar Graph for the percent of time WPD loop is open during equilibrium
#fig = plt.figure(figsize=(25,15))
#ax1 = fig.add_subplot(111)
#ax1.set_title("Time Spent Open", fontsize=24)
#ax1.set_ylabel('% Greater than 1nm', fontsize = 18)
#ax1.bar(num, per_eq, color = ['darkblue', 'blue', 'darkgreen', 'green', 'darkred', 'red', 'darkorange', 'orange', 'purple', 'magenta', 'brown', 'firebrick', 'gold', 'yellow'], width=4.5)
#plt.xticks(num, Method, fontsize=12)
#fig.savefig('WPD_eq_per.png')

#WPD Loop Distance for 1sug_AD opening
time = np.array(t_1sug_AD)/1000
fig2 = plt.figure()
ax1 = fig2.add_subplot(111)
ax1.set_title("WPD Loop Opening with AD and Ordered a7")    
ax1.set_xlabel('Time (ns)')
ax1.set_ylabel('Residue Distances (nm)')
ax1.plot(time,w_1sug_AD, label='1sug_AD')
fig2.savefig('1sug_opening.png')

#WPD Loop Distance for 1sug_dis_AD opening
time = np.array(t_1sug_dis9_AD)/1000
fig3 = plt.figure()
ax1 = fig3.add_subplot(111)
ax1.set_title("WPD Loop Opening with AD and Disordered a7")    
ax1.set_xlabel('Time (ns)')
ax1.set_ylabel('Residue Distances (nm)')
ax1.plot(time,w_1sug_dis9_AD, label='1sug_dis9_AD')
fig3.savefig('1sug_dis_AD_opening.png')

#WPD Loop Distance for 1sug_dis opening
time = np.array(t_1sug_dis11)/1000
fig4 = plt.figure()
ax1 = fig4.add_subplot(111)
ax1.set_title("WPD Loop Opening with AD and Disordered a7")    
ax1.set_xlabel('Time (ns)')
ax1.set_ylabel('Residue Distances (nm)')
ax1.plot(time,w_1sug_dis11, label='1sug_dis11')
fig4.savefig('1sug_dis_opening.png')

#Histograph of all WPD lengths
wpd_all = w_a7
wpd_all.extend(w_a7_AD)
wpd_all.extend(w_Apo)
wpd_all.extend(w_AD)
wpd_all.extend(w_1sug)
wpd_all.extend(w_1sug_AD)
wpd_all.extend(w_1sug_na7), 
wpd_all.extend(w_1sug_na7_AD)
wpd_all.extend(w_1sug_dis7)
wpd_all.extend(w_1sug_dis7_AD)
wpd_all.extend(w_1sug_dis9)
wpd_all.extend(w_1sug_dis9_AD)
wpd_all.extend(w_1sug_dis11), 
wpd_all.extend(w_1sug_dis11_AD)
wpd_all.extend(w_Apo_dis7)
wpd_all.extend(w_AD_dis7)
wpd_all.extend(w_Apo_dis9)
wpd_all.extend(w_AD_dis9)
wpd_all.extend(w_Apo_dis11)
wpd_all.extend(w_AD_dis11)
wpd_all.extend(w_1sug_alt)
wpd_all.extend(w_1sug_alt2)

fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.set_title('Histogram of WPD Loop Distances')
ax1.set_xlabel('Distance (nm)')
n, bin1, patches = ax1.hist(wpd_all, bins = 20, color = 'blue', density=True)
for i in range(6):
    patches[i].set_fc('red')
fig.savefig('WPD_hist.png')

#Histograph of Apo vs AD open lengths
wpd_apo = w_a7_eq
wpd_apo.extend(w_1sug_dis7_eq)
wpd_apo.extend(w_Apo_dis9_eq)

wpd_AD = w_1sug_dis11_AD_eq
wpd_AD.extend(w_1sug_alt_eq)
wpd_AD.extend(w_1sug_alt2_eq)

fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.set_title('WPD Loop Distances Compared Between Apo and Ligand Bound')
ax1.set_xlabel('Distance (nm)', fontsize = 14)
ax1.hist(wpd_AD, bins = 25, alpha = 0.5, color = 'black', density = True, label = 'AD')
ax1.hist(wpd_apo, bins = 25, alpha = 0.5, color = 'blue', density = True, label = 'Apo')
ax1.legend(loc = 'upper right')
fig.savefig('WPD_AD_cmpr_hist.png')

