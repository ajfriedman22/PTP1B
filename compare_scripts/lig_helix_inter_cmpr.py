from matplotlib import pyplot as plt
import numpy as np
from scipy import stats

#List of all directory paths for each group
dir_path_AD_crys = ['mutate/WT/AD/analysis', '1sug_dis_AD/analysis/config11']
dir_path_AD_alt = ['1sug_dis_AD/analysis/config_alt', '1sug_dis_AD/analysis/config_alt2']
dir_path_BBR = ['mutate/WT/BBR/analysis', 'BBR_a7/analysis']

#List interactions of interest
inters = ['a3', 'a4', 'a5', 'a6', 'a7']

#Set all empty arrays
a3_AD_crys, a4_AD_crys, a5_AD_crys, a6_AD_crys, a7_AD_crys = [],[],[],[],[]
a3_AD_alt, a4_AD_alt, a5_AD_alt, a6_AD_alt, a7_AD_alt = [],[],[],[],[]
a3_BBR, a4_BBR, a5_BBR, a6_BBR, a7_BBR = [],[],[],[],[]

#Load all data
for i in range(len(dir_path_AD_crys)):
    for j in range(len(inters)):
        data = open('../../../' + dir_path_AD_crys[i] + '/' + inters[j] + '_inter.txt').readlines()
        if j == 0:
            for k in data:
                a3_AD_crys.append(float(k))
        if j == 1:
            for k in data:
                a4_AD_crys.append(float(k))
        if j == 2:
            for k in data:
                a5_AD_crys.append(float(k))
        if j == 3:
            for k in data:
                a6_AD_crys.append(float(k))
        if j == 4:
            for k in data:
                a7_AD_crys.append(float(k))

for i in range(len(dir_path_AD_alt)):
    for j in range(len(inters)):
        data = open('../../../' + dir_path_AD_alt[i] + '/' + inters[j] + '_inter.txt').readlines()
        if j == 0:
            for k in data:
                a3_AD_alt.append(float(k))
        if j == 1:
            for k in data:
                a4_AD_alt.append(float(k))
        if j == 2:
            for k in data:
                a5_AD_alt.append(float(k))
        if j == 3:
            for k in data:
                a6_AD_alt.append(float(k))
        if j == 4:
            for k in data:
                a7_AD_alt.append(float(k))

for i in range(len(dir_path_BBR)):
    for j in range(len(inters)):
        data = open('../../../' + dir_path_BBR[i] + '/' + inters[j] + '_inter.txt').readlines()
        if j == 0:
            for k in data:
                a3_BBR.append(float(k))
        if j == 1:
            for k in data:
                a4_BBR.append(float(k))
        if j == 2:
            for k in data:
                a5_BBR.append(float(k))
        if j == 3:
            for k in data:
                a6_BBR.append(float(k))
        if j == 4:
            for k in data:
                a7_BBR.append(float(k))

#Determine the mean and sem for each helix interactions
a3_AD = np.append(a3_AD_alt, a3_AD_crys)
a4_AD = np.append(a4_AD_alt, a4_AD_crys)
a5_AD = np.append(a5_AD_alt, a5_AD_crys)
a6_AD = np.append(a6_AD_alt, a6_AD_crys)
a7_AD = np.append(a7_AD_alt, a7_AD_crys)

helix_a3 = np.array([np.mean(a3_AD), np.mean(a3_BBR)])
helix_a4 = np.array([np.mean(a4_AD), np.mean(a4_BBR)])
helix_a5 = np.array([np.mean(a5_AD), np.mean(a5_BBR)])
helix_a6 = np.array([np.mean(a6_AD), np.mean(a6_BBR)])
helix_a7 = np.array([np.mean(a7_AD), np.mean(a7_BBR)])

err_a3 = np.array([stats.sem(a3_AD), stats.sem(a3_BBR)])
err_a4 = np.array([stats.sem(a4_AD), stats.sem(a4_BBR)])
err_a5 = np.array([stats.sem(a5_AD), stats.sem(a5_BBR)])
err_a6 = np.array([stats.sem(a6_AD), stats.sem(a6_BBR)])
err_a7 = np.array([stats.sem(a7_AD), stats.sem(a7_BBR)])

#Run Welch's test on BBR vs AD helix interactions for a3-a7
st_a3, p_a3 = stats.ttest_ind(a3_AD, a3_BBR, equal_var = False) #Welch's t-test b/w a3 interactions with AD binding and BBR binding
st_a4, p_a4 = stats.ttest_ind(a4_AD, a4_BBR, equal_var = False) #Welch's t-test b/w a4 interactions with AD binding and BBR binding
st_a5, p_a5 = stats.ttest_ind(a5_AD, a5_BBR, equal_var = False) #Welch's t-test b/w a4 interactions with AD binding and BBR binding
st_a6, p_a6 = stats.ttest_ind(a6_AD, a6_BBR, equal_var = False) #Welch's t-test b/w a6 interactions with AD binding and BBR binding
st_a7, p_a7 = stats.ttest_ind(a7_AD, a7_BBR, equal_var = False) #Welch's t-test b/w a7 interactions with AD binding and BBR binding

#Print output to file
output = open('Helix_inter.txt', 'w')
output.write('a3 Helix Interactions\n')
output.write('Full for AD: ' + str(helix_a3[0]) + ' +/- ' + str(err_a3[0]) + '\n')
output.write('Full for BBR: ' + str(helix_a3[1]) + ' +/- ' + str(err_a3[1]) + '\n')
output.write('p-value b/w AD and BBR: ' + str(p_a3) + '\n')

output.write('a4 Helix Interactions\n')
output.write('Full for AD: ' + str(helix_a4[0]) + ' +/- ' + str(err_a4[0]) + '\n')
output.write('Full for BBR: ' + str(helix_a4[1]) + ' +/- ' + str(err_a4[1]) + '\n')
output.write('p-value b/w AD and BBR: ' + str(p_a4) + '\n')

output.write('a5 Helix Interactions\n')
output.write('AD crys: ' + str(helix_a5[0]) + ' +/- ' + str(err_a5[0]) + '\n')
output.write('BBR: ' + str(helix_a5[1]) + ' +/- ' + str(err_a5[1]) + '\n')
output.write('p-value b/w AD and BBR: ' + str(p_a5) + '\n')

output.write('a6 Helix Interactions\n')
output.write('AD crys: ' + str(helix_a6[0]) + ' +/- ' + str(err_a6[0]) + '\n')
output.write('BBR: ' + str(helix_a6[1]) + ' +/- ' + str(err_a6[1]) + '\n')
output.write('p-value b/w AD and BBR: ' + str(p_a6) + '\n')

output.write('a7 Helix Interactions\n')
output.write('AD crys: ' + str(helix_a7[0]) + ' +/- ' + str(err_a7[0]) + '\n')
output.write('BBR: ' + str(helix_a7[1]) + ' +/- ' + str(err_a7[1]) + '\n')
output.write('p-value b/w AD and BBR: ' + str(p_a7) + '\n')

#Plot comparison of helix interactions for a3-a7
num = [1, 2]
label = ['AD', 'BBR']
#Plot Bar graph comparing averages for a3 helix
fig = plt.figure()
plt.title('Interactions with a3 Helix')    
plt.ylabel('Mean Number of Interactions')
plt.bar(num, helix_a3[:], color = ['green', 'blue'], width=0.8)
plt.errorbar(num, helix_a3[:], yerr=err_a3[:], fmt='o', color='black')
plt.xticks(num, label, fontsize=8)
if p_a3 < 0.05 and p_a3 > 0.01:
    x1, x2 = 1, 2 #Columns for AD crystal and BBR
    y, h, col = (1.1*helix_a3.max()) + 2, 2, 'k'
    plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
    plt.text((x1+x2)*0.5, y+h, "*" , ha='center', va='bottom', color=col)
if p_a3 < 0.01 and p_a3 > 0.001:
    x1, x2 = 1, 2 #Columns for AD crystal and BBR
    y, h, col = (1.1*helix_a3.max()) + 2, 2, 'k'
    plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
    plt.text((x1+x2)*0.5, y+h, "**" , ha='center', va='bottom', color=col)
if p_a3 < 0.001:
    x1, x2 = 1, 2 #Columns for AD crystal and BBR
    y, h, col = (1.1*helix_a3.max()) + 2, 2, 'k'
    plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
    plt.text((x1+x2)*0.5, y+h, "***" , ha='center', va='bottom', color=col)
fig.savefig('a3_helix_inter_cmpr.png')
plt.close(fig)

#Plot Bar graph comparing averages for a4 helix
fig = plt.figure()
plt.title('Interactions with a4 Helix')    
plt.ylabel('Mean Number of Interactions')
plt.bar(num, helix_a4[:], color = ['green', 'blue'], width=0.8)
plt.errorbar(num, helix_a4[:], yerr=err_a4[:], fmt='o', color='black')
plt.xticks(num, label, fontsize=8)
if p_a4 < 0.05 and p_a4 > 0.01:
    x1, x2 = 1, 2 #Columns for AD crystal and BBR
    y, h, col = (1.1*helix_a4.max()) + 2, 2, 'k'
    plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
    plt.text((x1+x2)*0.5, y+h, "*" , ha='center', va='bottom', color=col)
if p_a4 < 0.01 and p_a4 > 0.001:
    x1, x2 = 1, 2 #Columns for AD crystal and BBR
    y, h, col = (1.1*helix_a4.max()) + 2, 2, 'k'
    plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
    plt.text((x1+x2)*0.5, y+h, "**" , ha='center', va='bottom', color=col)
if p_a4 < 0.001:
    x1, x2 = 1, 2 #Columns for AD crystal and BBR
    y, h, col = (1.1*helix_a4.max()) + 2, 2, 'k'
    plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
    plt.text((x1+x2)*0.5, y+h, "***" , ha='center', va='bottom', color=col)
fig.savefig('a4_helix_inter_cmpr.png')
plt.close(fig)

#Plot Bar graph comparing averages for a5 helix
fig = plt.figure()
plt.title('Interactions with a5 Helix')    
plt.ylabel('Mean Number of Interactions')
plt.bar(num, helix_a5[:], color = ['green', 'blue'], width=0.8)
plt.errorbar(num, helix_a5[:], yerr=err_a5[:], fmt='o', color='black')
plt.xticks(num, label, fontsize=8)
if p_a5 < 0.05 and p_a5 > 0.01:
    x1, x2 = 1, 2 #Columns for AD crystal and BBR
    y, h, col = (1.1*helix_a5.max()) + 2, 2, 'k'
    plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
    plt.text((x1+x2)*0.5, y+h, "*" , ha='center', va='bottom', color=col)
if p_a5 < 0.01 and p_a5 > 0.001:
    x1, x2 = 1, 2 #Columns for AD crystal and BBR
    y, h, col = (1.1*helix_a5.max()) + 2, 2, 'k'
    plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
    plt.text((x1+x2)*0.5, y+h, "**" , ha='center', va='bottom', color=col)
if p_a5 < 0.001:
    x1, x2 = 1, 2 #Columns for AD crystal and BBR
    y, h, col = (1.1*helix_a5.max()) + 2, 2, 'k'
    plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
    plt.text((x1+x2)*0.5, y+h, "***" , ha='center', va='bottom', color=col)
fig.savefig('a5_helix_inter_cmpr.png')
plt.close(fig)

#Plot Bar graph comparing averages for a6 helix
fig = plt.figure()
plt.title('Interactions with a6 Helix')    
plt.ylabel('Mean Number of Interactions')
plt.bar(num, helix_a6[:], color = ['green', 'blue'], width=0.8)
plt.errorbar(num, helix_a6[:], yerr=err_a6[:], fmt='o', color='black')
plt.xticks(num, label, fontsize=8)
if p_a6 < 0.05 and p_a6 > 0.01:
    x1, x2 = 1, 2 #Columns for AD crystal and BBR
    y, h, col = (1.1*helix_a6.max()) + 2, 2, 'k'
    plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
    plt.text((x1+x2)*0.5, y+h, "*" , ha='center', va='bottom', color=col)
if p_a6 < 0.01 and p_a6 > 0.001:
    x1, x2 = 1, 2 #Columns for AD crystal and BBR
    y, h, col = (1.1*helix_a6.max()) + 2, 2, 'k'
    plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
    plt.text((x1+x2)*0.5, y+h, "**" , ha='center', va='bottom', color=col)
if p_a6 < 0.001:
    x1, x2 = 1, 2 #Columns for AD crystal and BBR
    y, h, col = (1.1*helix_a6.max()) + 2, 2, 'k'
    plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
    plt.text((x1+x2)*0.5, y+h, "***" , ha='center', va='bottom', color=col)
fig.savefig('a6_helix_inter_cmpr.png')
plt.close(fig)

#Plot Bar graph comparing averages for a7 helix
fig = plt.figure()
plt.title('Interactions with a7 Helix')    
plt.ylabel('Mean Number of Interactions')
plt.bar(num, helix_a7[:], color = ['green', 'blue'], width=0.8)
plt.errorbar(num, helix_a7[:], yerr=err_a7[:], fmt='o', color='black')
plt.xticks(num, label, fontsize=8)
if p_a7 < 0.05 and p_a7 > 0.01:
    x1, x2 = 1, 2 #Columns for AD crystal and BBR
    y, h, col = (1.1*helix_a7.max()) + 2, 2, 'k'
    plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
    plt.text((x1+x2)*0.5, y+h, "*" , ha='center', va='bottom', color=col)
if p_a7 < 0.01 and p_a7 > 0.001:
    x1, x2 = 1, 2 #Columns for AD crystal and BBR
    y, h, col = (1.1*helix_a7.max()) + 2, 2, 'k'
    plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
    plt.text((x1+x2)*0.5, y+h, "**" , ha='center', va='bottom', color=col)
if p_a7 < 0.001:
    x1, x2 = 1, 2 #Columns for AD crystal and BBR
    y, h, col = (1.1*helix_a7.max()) + 2, 2, 'k'
    plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
    plt.text((x1+x2)*0.5, y+h, "***" , ha='center', va='bottom', color=col)
fig.savefig('a7_helix_inter_cmpr.png')
plt.close(fig)

