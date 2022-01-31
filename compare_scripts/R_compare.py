#Import Necessary Packages
import numpy as np
import matplotlib.pyplot as plt

#Make open arrays for time and atomic distances
t_a7, t_a7_AD, t_Apo, t_AD, t_1sug, t_1sug_AD, t_1sug_na7, t_1sug_na7_AD = [],[],[],[],[],[],[],[] #time
r_a7, r_a7_AD, r_Apo, r_AD, r_1sug, r_1sug_AD, r_1sug_na7, r_1sug_na7_AD = [],[],[],[],[],[],[],[] #atom distance

#Count variable for number of time steps the WPD loop is open
count = 0

#empty vector for the percent of time the WPD loop is open
per = []

#Input Data
with open("../../rebuild_a7/a7_R.xvg") as f:
    for _ in range(17):
        next(f)
    for line in f:
        cols = line.split()
        if len(cols) == 2:
            t_a7.append(float(cols[0]))
            r_a7.append(float(cols[1]))
with open("../../AD_rebuild_a7/a7_AD_R.xvg") as f:
    for _ in range(17):
        next(f)
    for line in f:
        cols = line.split()
        if len(cols) == 2:
            t_a7_AD.append(float(cols[0]))
            r_a7_AD.append(float(cols[1]))
with open("../../Apo/Apo_R.xvg") as f:
    for _ in range(17):
        next(f)
    for line in f:
        cols = line.split()
        if len(cols) == 2:
            t_Apo.append(float(cols[0]))
            r_Apo.append(float(cols[1]))
with open("../../AD/AD_R.xvg") as f:
    for _ in range(17):
        next(f)
    for line in f:
        cols = line.split()
        if len(cols) == 2:
            t_AD.append(float(cols[0]))
            r_AD.append(float(cols[1]))
with open("../../1sug/1sug_R.xvg") as f:
    for _ in range(17):
        next(f)
    for line in f:
        cols = line.split()
        if len(cols) == 2:
            t_1sug.append(float(cols[0]))
            r_1sug.append(float(cols[1]))
with open("../../1sug_AD/1sug_AD_R.xvg") as f:
    for _ in range(17):
        next(f)
    for line in f:
        cols = line.split()
        if len(cols) == 2:
            t_1sug_AD.append(float(cols[0]))
            r_1sug_AD.append(float(cols[1]))
with open("../../1sug_no_a7/1sug_na7_R.xvg") as f:
    for _ in range(17):
        next(f)
    for line in f:
        cols = line.split()
        if len(cols) == 2:
            t_1sug_na7.append(float(cols[0]))
            r_1sug_na7.append(float(cols[1]))
with open("../../1sug_no_a7_AD/1sug_na7_AD_R.xvg") as f:
    for _ in range(17):
        next(f)
    for line in f:
        cols = line.split()
        if len(cols) == 2:
            t_1sug_na7_AD.append(float(cols[0]))
            r_1sug_na7_AD.append(float(cols[1]))

#Determine the Percent of Time the WPD loop is open
for i in range(len(t_a7)):
    if r_a7[i]>1.2:
        count +=1
per.append(100*count/len(t_a7))
count=0
for i in range(len(t_a7_AD)):
    if r_a7_AD[i]>1.2:
        count +=1
per.append(100*count/len(t_a7_AD))
count=0
for i in range(len(t_Apo)):
    if r_Apo[i]>1.2:
        count +=1
per.append(100*count/len(t_Apo))
count=0
for i in range(len(t_AD)):
    if r_AD[i]>1.2:
        count +=1
per.append(100*count/len(t_AD))
count=0
for i in range(len(t_1sug)):
    if r_1sug[i]>1.2:
        count +=1
per.append(100*count/len(t_1sug))
count=0
for i in range(len(t_1sug_AD)):
    if r_1sug_AD[i]>1.2:
        count +=1
per.append(100*count/len(t_1sug_AD))
count=0
for i in range(len(t_1sug_na7)):
    if r_1sug_na7[i]>1.2:
        count +=1
per.append(100*count/len(t_1sug_na7))
count=0
for i in range(len(t_1sug_na7_AD)):
    if r_1sug_na7_AD[i]>1.2:
        count +=1
per.append(100*count/len(t_1sug_na7_AD))

#Plot Bar Graph for the percent of time WPD loop is open
num = [5,10,15,20,25,30,35,40]
Method = ['a7', 'a7_AD', 'Apo', 'AD', '1sug', '1sug_AD', '1sug_no_a7', '1sug_no_a7_AD']
fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.set_title("Time Spent Open")    
ax1.set_ylabel('% Greater than 1.2nm')
ax1.bar(num,per, color = ['darkblue', 'blue', 'darkgreen', 'green', 'darkred', 'red', 'darkorange', 'orange'], width=4.5)
plt.xticks(num, Method, fontsize=8)
fig.savefig('R_percent.png')

#WPD Loop Distance for a7_AD opening
time = np.array(t_1sug_AD)/1000
time_a7 = np.array(t_a7)/1000
time_a7_AD = np.array(t_a7_AD)/1000
time_Apo = np.array(t_Apo)/1000
time_AD = np.array(t_AD)/1000
time_1sug = np.array(t_1sug)/1000
time_1sug_na7 = np.array(t_1sug_na7)/1000
time_1sug_na7_AD = np.array(t_1sug_na7_AD)/1000

fig2 = plt.figure()
ax1 = fig2.add_subplot(111)
ax1.set_title("Comparison of R Loop Distances")    
ax1.set_xlabel('Time (ps)')
ax1.set_ylabel('Residue Distances (nm)')
ax1.plot(time_a7, r_a7, label='a7')
ax1.plot(time_a7_AD, r_a7_AD, label='a7_AD')
ax1.plot(time_Apo, r_Apo, label='Apo')
ax1.plot(time_AD, r_AD, label='AD')
ax1.plot(time_1sug, r_1sug, label='1sug')
ax1.plot(time,r_1sug_AD, label='1sug_AD')
ax1.plot(time_1sug_na7, r_1sug_na7, label='1sug_no_a7')
ax1.plot(time_1sug_na7_AD,r_1sug_na7_AD, label='1sug_no_a7_AD')
ax1.legend()
fig2.savefig('full_R_compare.png')

#Print max values
print('a7 max = ' + str(max(r_a7)))
print('a7+AD max = ' + str(max(r_a7_AD)))
print('Apo max = ' + str(max(r_Apo)))
print('AD max = ' + str(max(r_AD)))
print('1sug max = ' + str(max(r_1sug)))
print('1sug+AD max = ' + str(max(r_1sug_AD)))
print('1sug+na7 max = ' + str(max(r_1sug_na7)))
print('1sug+na7+AD max = ' + str(max(r_1sug_na7_AD)))


