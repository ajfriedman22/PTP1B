from matplotlib import pyplot as plt
import numpy as np
from scipy import stats

#Black arrays for contacts in each file
a3_pt1_AD_crys, a3_pt2_AD_crys, a4_pt1_AD_crys, a4_pt2_AD_crys, a5_AD_crys, a6_pt1_AD_crys, a6_pt2_AD_crys, a7_pt1_AD_crys, a7_pt2_AD_crys = [],[],[],[],[],[],[],[],[]
a3_pt1_AD_alt1, a3_pt2_AD_alt1, a4_pt1_AD_alt1, a4_pt2_AD_alt1, a5_AD_alt1, a6_pt1_AD_alt1, a6_pt2_AD_alt1, a7_pt1_AD_alt1, a7_pt2_AD_alt1 = [],[],[],[],[],[],[],[],[]
a3_pt1_BBR, a3_pt2_BBR, a4_pt1_BBR, a4_pt2_BBR, a5_BBR, a6_pt1_BBR, a6_pt2_BBR, a7_pt1_BBR, a7_pt2_BBR = [],[],[],[],[],[],[],[],[]

#Open input files or contacts over time
for i in open('../../../1sug_dis_AD/config11/a3_pt1_inter.txt', 'r').readlines():
    a3_pt1_AD_crys.append(float(i))
for i in open('../../../1sug_dis_AD/config11/a3_pt2_inter.txt', 'r').readlines():
    a3_pt2_AD_crys.append(float(i))
for i in open('../../../1sug_dis_AD/config11/a4_pt1_inter.txt', 'r').readlines():
    a4_pt1_AD_crys.append(float(i))
for i in open('../../../1sug_dis_AD/config11/a4_pt2_inter.txt', 'r').readlines():
    a4_pt2_AD_crys.append(float(i))
for i in open('../../../1sug_dis_AD/config11/a5_inter.txt', 'r').readlines():
    a5_AD_crys.append(float(i))
for i in open('../../../1sug_dis_AD/config11/a6_pt1_inter.txt', 'r').readlines():
    a6_pt1_AD_crys.append(float(i))
for i in open('../../../1sug_dis_AD/config11/a6_pt2_inter.txt', 'r').readlines():
    a6_pt2_AD_crys.append(float(i))
for i in open('../../../1sug_dis_AD/config11/a7_pt1_inter.txt', 'r').readlines():
    a7_pt1_AD_crys.append(float(i))
for i in open('../../../1sug_dis_AD/config11/a7_pt2_inter.txt', 'r').readlines():
    a7_pt2_AD_crys.append(float(i))

for i in open('../../../1sug_AD_dis_alt/run_1/a3_pt1_inter.txt', 'r').readlines():
    a3_pt1_AD_alt1.append(float(i))
for i in open('../../../1sug_AD_dis_alt/run_1/a3_pt2_inter.txt', 'r').readlines():
    a3_pt2_AD_alt1.append(float(i))
for i in open('../../../1sug_AD_dis_alt/run_1/a4_pt1_inter.txt', 'r').readlines():
    a4_pt1_AD_alt1.append(float(i))
for i in open('../../../1sug_AD_dis_alt/run_1/a4_pt2_inter.txt', 'r').readlines():
    a4_pt2_AD_alt1.append(float(i))
for i in open('../../../1sug_AD_dis_alt/run_1/a5_inter.txt', 'r').readlines():
    a5_AD_alt1.append(float(i))
for i in open('../../../1sug_AD_dis_alt/run_1/a6_pt1_inter.txt', 'r').readlines():
    a6_pt1_AD_alt1.append(float(i))
for i in open('../../../1sug_AD_dis_alt/run_1/a6_pt2_inter.txt', 'r').readlines():
    a6_pt2_AD_alt1.append(float(i))
for i in open('../../../1sug_AD_dis_alt/run_1/a7_pt1_inter.txt', 'r').readlines():
    a7_pt1_AD_alt1.append(float(i))
for i in open('../../../1sug_AD_dis_alt/run_1/a7_pt2_inter.txt', 'r').readlines():
    a7_pt2_AD_alt1.append(float(i))

for i in open('../../../1sug_AD_dis_alt/run_2/a3_pt1_inter.txt', 'r').readlines():
    a3_pt1_AD_alt1.append(float(i))
for i in open('../../../1sug_AD_dis_alt/run_2/a3_pt2_inter.txt', 'r').readlines():
    a3_pt2_AD_alt1.append(float(i))
for i in open('../../../1sug_AD_dis_alt/run_2/a4_pt1_inter.txt', 'r').readlines():
    a4_pt1_AD_alt1.append(float(i))
for i in open('../../../1sug_AD_dis_alt/run_2/a4_pt2_inter.txt', 'r').readlines():
    a4_pt2_AD_alt1.append(float(i))
for i in open('../../../1sug_AD_dis_alt/run_2/a5_inter.txt', 'r').readlines():
    a5_AD_alt1.append(float(i))
for i in open('../../../1sug_AD_dis_alt/run_2/a6_pt1_inter.txt', 'r').readlines():
    a6_pt1_AD_alt1.append(float(i))
for i in open('../../../1sug_AD_dis_alt/run_2/a6_pt2_inter.txt', 'r').readlines():
    a6_pt2_AD_alt1.append(float(i))
for i in open('../../../1sug_AD_dis_alt/run_2/a7_pt1_inter.txt', 'r').readlines():
    a7_pt1_AD_alt1.append(float(i))
for i in open('../../../1sug_AD_dis_alt/run_2/a7_pt2_inter.txt', 'r').readlines():
    a7_pt2_AD_alt1.append(float(i))


for i in open('../../../BBR_a7/a3_pt1_inter.txt', 'r').readlines():
    a3_pt1_BBR.append(float(i))
for i in open('../../../BBR_a7/a3_pt2_inter.txt', 'r').readlines():
    a3_pt2_BBR.append(float(i))
for i in open('../../../BBR_a7/a4_pt1_inter.txt', 'r').readlines():
    a4_pt1_BBR.append(float(i))
for i in open('../../../BBR_a7/a4_pt2_inter.txt', 'r').readlines():
    a4_pt2_BBR.append(float(i))
for i in open('../../../BBR_a7/a5_inter.txt', 'r').readlines():
    a5_BBR.append(float(i))
for i in open('../../../BBR_a7/a6_pt1_inter.txt', 'r').readlines():
    a6_pt1_BBR.append(float(i))
for i in open('../../../BBR_a7/a6_pt2_inter.txt', 'r').readlines():
    a6_pt2_BBR.append(float(i))
for i in open('../../../BBR_a7/a7_pt1_inter.txt', 'r').readlines():
    a7_pt1_BBR.append(float(i))
for i in open('../../../BBR_a7/a7_pt2_inter.txt', 'r').readlines():
    a7_pt2_BBR.append(float(i))

#Determine the mean and sem for each helix interactions
a3_AD = np.add(a3_pt1_AD_crys, a3_pt2_AD_crys)
a4_AD = np.add(a4_pt1_AD_crys, a4_pt2_AD_crys)
a6_AD = np.add(a6_pt1_AD_crys, a6_pt2_AD_crys)
a7_AD = np.add(a7_pt1_AD_crys, a7_pt2_AD_crys)

a3_AD = np.append(a3_AD, np.add(a3_pt1_AD_alt1, a3_pt2_AD_alt1))
a4_AD = np.append(a4_AD, np.add(a4_pt1_AD_alt1, a4_pt2_AD_alt1))
a5_AD = np.append(a5_AD_crys, a5_AD_alt1)
a6_AD = np.append(a6_AD, np.add(a6_pt1_AD_alt1, a6_pt2_AD_alt1))
a7_AD = np.append(a7_AD, np.add(a7_pt1_AD_alt1, a7_pt2_AD_alt1))

a3_BBR = np.add(a3_pt1_BBR, a3_pt2_BBR)
a4_BBR = np.add(a4_pt1_BBR, a4_pt2_BBR)
a6_BBR = np.add(a6_pt1_BBR, a6_pt2_BBR)
a7_BBR = np.add(a7_pt1_BBR, a7_pt2_BBR)

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

helix_a3_pt1 = np.array([np.mean(a3_pt1_AD_crys), np.mean(a3_pt1_AD_alt1), np.mean(a3_pt1_BBR)])
helix_a4_pt1 = np.array([np.mean(a4_pt1_AD_crys), np.mean(a4_pt1_AD_alt1), np.mean(a4_pt1_BBR)])
helix_a6_pt1 = np.array([np.mean(a6_pt1_AD_crys), np.mean(a6_pt1_AD_alt1), np.mean(a6_pt1_BBR)])
helix_a7_pt1 = np.array([np.mean(a7_pt1_AD_crys), np.mean(a7_pt1_AD_alt1), np.mean(a7_pt1_BBR)])

err_a3_pt1 = np.array([stats.sem(a3_pt1_AD_crys), stats.sem(a3_pt1_AD_alt1), stats.sem(a3_pt1_BBR)])
err_a4_pt1 = np.array([stats.sem(a4_pt1_AD_crys), stats.sem(a4_pt1_AD_alt1), stats.sem(a4_pt1_BBR)])
err_a6_pt1 = np.array([stats.sem(a6_pt1_AD_crys), stats.sem(a6_pt1_AD_alt1), stats.sem(a6_pt1_BBR)])
err_a7_pt1 = np.array([stats.sem(a7_pt1_AD_crys), stats.sem(a7_pt1_AD_alt1), stats.sem(a7_pt1_BBR)])

helix_a3_pt2 = np.array([np.mean(a3_pt2_AD_crys), np.mean(a3_pt2_AD_alt1), np.mean(a3_pt2_BBR)])
helix_a4_pt2 = np.array([np.mean(a4_pt2_AD_crys), np.mean(a4_pt2_AD_alt1), np.mean(a4_pt2_BBR)])
helix_a6_pt2 = np.array([np.mean(a6_pt2_AD_crys), np.mean(a6_pt2_AD_alt1), np.mean(a6_pt2_BBR)])
helix_a7_pt2 = np.array([np.mean(a7_pt2_AD_crys), np.mean(a7_pt2_AD_alt1), np.mean(a7_pt2_BBR)])

err_a3_pt2 = np.array([stats.sem(a3_pt2_AD_crys), stats.sem(a3_pt2_AD_alt1), stats.sem(a3_pt2_BBR)])
err_a4_pt2 = np.array([stats.sem(a4_pt2_AD_crys), stats.sem(a4_pt2_AD_alt1), stats.sem(a4_pt2_BBR)])
err_a6_pt2 = np.array([stats.sem(a6_pt2_AD_crys), stats.sem(a6_pt2_AD_alt1), stats.sem(a6_pt2_BBR)])
err_a7_pt2 = np.array([stats.sem(a7_pt2_AD_crys), stats.sem(a7_pt2_AD_alt1), stats.sem(a7_pt2_BBR)])

#Run Welch's test on BBR vs AD helix interactions for a3-a7
st_a3, p_a3 = stats.ttest_ind(a3_AD, a3_BBR, equal_var = False) #Welch's t-test b/w a3 interactions with AD binding and BBR binding
st_a4, p_a4 = stats.ttest_ind(a4_AD, a4_BBR, equal_var = False) #Welch's t-test b/w a4 interactions with AD binding and BBR binding
st_a5, p_a5 = stats.ttest_ind(a5_AD, a5_BBR, equal_var = False) #Welch's t-test b/w a4 interactions with AD binding and BBR binding
st_a6, p_a6 = stats.ttest_ind(a6_AD, a6_BBR, equal_var = False) #Welch's t-test b/w a6 interactions with AD binding and BBR binding
st_a7, p_a7 = stats.ttest_ind(a7_AD, a7_BBR, equal_var = False) #Welch's t-test b/w a7 interactions with AD binding and BBR binding

st_a3_pt1, p_a3_pt1 = stats.ttest_ind(a3_pt1_AD_crys, a3_pt1_BBR, equal_var = False) #Welch's t-test b/w a3 pt1 interacts with AD crystal binding and BBR binding
st_a4_pt1, p_a4_pt1 = stats.ttest_ind(a4_pt1_AD_crys, a4_pt1_BBR, equal_var = False) #Welch's t-test b/w a4 pt1 interacts with AD crystal binding and BBR binding
st_a6_pt1, p_a6_pt1 = stats.ttest_ind(a6_pt1_AD_crys, a6_pt1_BBR, equal_var = False) #Welch's t-test b/w a6 pt1 interacts with AD crystal binding and BBR binding
st_a7_pt1, p_a7_pt1 = stats.ttest_ind(a7_pt1_AD_crys, a7_pt1_BBR, equal_var = False) #Welch's t-test b/w a7 pt1 interacts with AD crystal binding and BBR binding

st_a3_alt, p_a3_alt = stats.ttest_ind(a3_pt1_AD_alt1, a3_pt1_BBR, equal_var = False) #Welch's t-test b/w a3 pt1 inter with AD alt1 + BBR binding
st_a4_alt, p_a4_alt = stats.ttest_ind(a4_pt1_AD_alt1, a4_pt1_BBR, equal_var = False) #Welch's t-test b/w a4 pt1 inter with AD alt1 and BBR binding
st_a6_alt, p_a6_alt = stats.ttest_ind(a6_pt1_AD_alt1, a6_pt1_BBR, equal_var = False) #Welch's t-test b/w a6 pt1 inter with AD alt1 and BBR binding
st_a7_alt, p_a7_alt = stats.ttest_ind(a7_pt1_AD_alt1, a7_pt1_BBR, equal_var = False) #Welch's t-test b/w a7 pt1 inter with AD alt1 and BBR binding

st_a3_pt2, p_a3_pt2 = stats.ttest_ind(a3_pt2_AD_crys, a3_pt2_BBR, equal_var = False) #Welch's t-test b/w a3 inter with AD crystal binding and BBR binding
st_a4_pt2, p_a4_pt2 = stats.ttest_ind(a4_pt2_AD_crys, a4_pt2_BBR, equal_var = False) #Welch's t-test b/w a4 inter with AD crystal binding and BBR binding
st_a6_pt2, p_a6_pt2 = stats.ttest_ind(a6_pt2_AD_crys, a6_pt2_BBR, equal_var = False) #Welch's t-test b/w a6 inter with AD crystal binding and BBR binding
st_a7_pt2, p_a7_pt2 = stats.ttest_ind(a7_pt2_AD_crys, a7_pt2_BBR, equal_var = False) #Welch's t-test b/w a7 inter with AD crystal binding and BBR binding

st_a3_pt2_alt, p_a3_pt2_alt = stats.ttest_ind(a3_pt2_AD_alt1, a3_pt2_BBR, equal_var = False) #Welch's t-test b/w a3 inter with AD alt1 bind + BBR binding
st_a4_pt2_alt, p_a4_pt2_alt = stats.ttest_ind(a4_pt2_AD_alt1, a4_pt2_BBR, equal_var = False) #Welch's t-test b/w a4 inter with AD alt1 bind + BBR binding
st_a6_pt2_alt, p_a6_pt2_alt = stats.ttest_ind(a6_pt2_AD_alt1, a6_pt2_BBR, equal_var = False) #Welch's t-test b/w a6 inter with AD alt1 bind + BBR binding
st_a7_pt2_alt, p_a7_pt2_alt = stats.ttest_ind(a7_pt2_AD_alt1, a7_pt2_BBR, equal_var = False) #Welch's t-test b/w a7 inter with AD alt1 bind + BBR binding

#Print output to file
output = open('Helix_inter.txt', 'w')
output.write('a3 Helix Interactions\n')
output.write('Full for AD: ' + str(helix_a3[0]) + ' +/- ' + str(err_a3[0]) + '\n')
output.write('Full for BBR: ' + str(helix_a3[1]) + ' +/- ' + str(err_a3[1]) + '\n')
output.write('p-value b/w AD and BBR: ' + str(p_a3) + '\n')
output.write('First Half for AD crys: ' + str(helix_a3_pt1[0]) + ' +/- ' + str(err_a3_pt1[0]) + '\n')
output.write('First Half for AD alt1: ' + str(helix_a3_pt1[1]) + ' +/- ' + str(err_a3_pt1[1]) + '\n')
output.write('First Half for BBR: ' + str(helix_a3_pt1[2]) + ' +/- ' + str(err_a3_pt1[2]) + '\n')
output.write('p-value b/w AD crys and BBR: ' + str(p_a3_pt1) + '\n')
output.write('p-value b/w AD alt1 and BBR: ' + str(p_a3_alt) + '\n')
output.write('Second Half for AD crys: ' + str(helix_a3_pt2[0]) + ' +/- ' + str(err_a3_pt2[0]) + '\n')
output.write('Second Half for AD alt1: ' + str(helix_a3_pt2[1]) + ' +/- ' + str(err_a3_pt2[1]) + '\n')
output.write('Second Half for BBR: ' + str(helix_a3_pt2[2]) + ' +/- ' + str(err_a3_pt2[2]) + '\n')
output.write('p-value b/w AD crys and BBR: ' + str(p_a3_pt2) + '\n')
output.write('p-value b/w AD alt1 and BBR: ' + str(p_a3_pt2_alt) + '\n')

output.write('a4 Helix Interactions\n')
output.write('Full for AD: ' + str(helix_a4[0]) + ' +/- ' + str(err_a4[0]) + '\n')
output.write('Full for BBR: ' + str(helix_a4[1]) + ' +/- ' + str(err_a4[1]) + '\n')
output.write('p-value b/w AD and BBR: ' + str(p_a4) + '\n')
output.write('First Half for AD crys: ' + str(helix_a4_pt1[0]) + ' +/- ' + str(err_a4_pt1[0]) + '\n')
output.write('First Half for AD alt1: ' + str(helix_a4_pt1[1]) + ' +/- ' + str(err_a4_pt1[1]) + '\n')
output.write('First Half for BBR: ' + str(helix_a4_pt1[2]) + ' +/- ' + str(err_a4_pt1[2]) + '\n')
output.write('p-value b/w AD crys and BBR: ' + str(p_a4_pt1) + '\n')
output.write('p-value b/w AD alt1 and BBR: ' + str(p_a4_alt) + '\n')
output.write('Second Half for AD crys: ' + str(helix_a4_pt2[0]) + ' +/- ' + str(err_a4_pt2[0]) + '\n')
output.write('Second Half for AD alt1: ' + str(helix_a4_pt2[1]) + ' +/- ' + str(err_a4_pt2[1]) + '\n')
output.write('Second Half for BBR: ' + str(helix_a4_pt2[2]) + ' +/- ' + str(err_a4_pt2[2]) + '\n')
output.write('p-value b/w AD crys and BBR: ' + str(p_a4_pt2) + '\n')
output.write('p-value b/w AD alt1 and BBR: ' + str(p_a4_pt2_alt) + '\n')

output.write('a5 Helix Interactions\n')
output.write('AD crys: ' + str(helix_a5[0]) + ' +/- ' + str(err_a5[0]) + '\n')
output.write('BBR: ' + str(helix_a5[1]) + ' +/- ' + str(err_a5[1]) + '\n')
output.write('p-value b/w AD and BBR: ' + str(p_a5) + '\n')

output.write('a6 Helix Interactions\n')
output.write('AD crys: ' + str(helix_a6[0]) + ' +/- ' + str(err_a6[0]) + '\n')
output.write('BBR: ' + str(helix_a6[1]) + ' +/- ' + str(err_a6[1]) + '\n')
output.write('p-value b/w AD and BBR: ' + str(p_a6) + '\n')
output.write('First Half for AD crys: ' + str(helix_a6_pt1[0]) + ' +/- ' + str(err_a6_pt1[0]) + '\n')
output.write('First Half for AD alt1: ' + str(helix_a6_pt1[1]) + ' +/- ' + str(err_a6_pt1[1]) + '\n')
output.write('First Half for BBR: ' + str(helix_a6_pt1[2]) + ' +/- ' + str(err_a6_pt1[2]) + '\n')
output.write('p-value b/w AD crys and BBR: ' + str(p_a6_pt1) + '\n')
output.write('p-value b/w AD alt1 and BBR: ' + str(p_a6_alt) + '\n')
output.write('Second Half for AD crys: ' + str(helix_a6_pt2[0]) + ' +/- ' + str(err_a6_pt2[0]) + '\n')
output.write('Second Half for AD alt1: ' + str(helix_a6_pt2[1]) + ' +/- ' + str(err_a6_pt2[1]) + '\n')
output.write('Second Half for BBR: ' + str(helix_a6_pt2[2]) + ' +/- ' + str(err_a6_pt2[2]) + '\n')
output.write('p-value b/w AD crys and BBR: ' + str(p_a6_pt2) + '\n')
output.write('p-value b/w AD alt1 and BBR: ' + str(p_a6_pt2_alt) + '\n')

output.write('a7 Helix Interactions\n')
output.write('AD crys: ' + str(helix_a7[0]) + ' +/- ' + str(err_a7[0]) + '\n')
output.write('BBR: ' + str(helix_a7[1]) + ' +/- ' + str(err_a7[1]) + '\n')
output.write('p-value b/w AD and BBR: ' + str(p_a7) + '\n')
output.write('First Half for AD crys: ' + str(helix_a7_pt1[0]) + ' +/- ' + str(err_a7_pt1[0]) + '\n')
output.write('First Half for AD alt1: ' + str(helix_a7_pt1[1]) + ' +/- ' + str(err_a7_pt1[1]) + '\n')
output.write('First Half for BBR: ' + str(helix_a7_pt1[2]) + ' +/- ' + str(err_a7_pt1[2]) + '\n')
output.write('p-value b/w AD crys and BBR: ' + str(p_a7_pt1) + '\n')
output.write('p-value b/w AD alt1 and BBR: ' + str(p_a7_alt) + '\n')
output.write('Second Half for AD crys: ' + str(helix_a7_pt2[0]) + ' +/- ' + str(err_a7_pt2[0]) + '\n')
output.write('Second Half for AD alt1: ' + str(helix_a7_pt2[1]) + ' +/- ' + str(err_a7_pt2[1]) + '\n')
output.write('Second Half for BBR: ' + str(helix_a7_pt2[2]) + ' +/- ' + str(err_a7_pt2[2]) + '\n')
output.write('p-value b/w AD crys and BBR: ' + str(p_a7_pt2) + '\n')
output.write('p-value b/w AD alt1 and BBR: ' + str(p_a7_pt2_alt) + '\n')

#Plot comparison of helix interactions for a3-a7
num = [1, 2, 3]
label = ['AD\n Crystal', 'AD\n alt1', 'BBR']

#Plot Bar graph comparing averages for i1st half of a3 helix
fig = plt.figure()
plt.title('Interactions with First Half of a3 Helix')    
plt.ylabel('Mean Number of Interactions')
plt.bar(num, helix_a3_pt1[:], color = ['green', 'green', 'blue'], width=0.8)
plt.errorbar(num, helix_a3_pt1[:], yerr=err_a3_pt1[:], fmt='o', color='black')
plt.xticks(num, label, fontsize=8)
if p_a3_pt1 < 0.05 and p_a3_pt1 > 0.01:
    x1, x2 = 1, 3 #Columns for AD crystal and BBR
    y, h, col = (1.1*helix_a3_pt1[[0,2]].max()) + 2, 2, 'k'
    plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
    plt.text((x1+x2)*0.5, y+h, "*" , ha='center', va='bottom', color=col)
if p_a3_pt1 < 0.01 and p_a3_pt1 > 0.001:
    x1, x2 = 1, 3 #Columns for AD crystal and BBR
    y, h, col = (1.1*helix_a3_pt1[[0,2]].max()) + 2, 2, 'k'
    plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
    plt.text((x1+x2)*0.5, y+h, "**" , ha='center', va='bottom', color=col)
if p_a3_pt1 < 0.001:
    x1, x2 = 1, 3 #Columns for AD crystal and BBR
    y, h, col = (1.1*helix_a3_pt1[[0,2]].max()) + 2, 2, 'k'
    plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
    plt.text((x1+x2)*0.5, y+h, "***" , ha='center', va='bottom', color=col)
if p_a3_alt < 0.05 and p_a3_alt > 0.01:
    x1, x2 = 2, 3 #Columns for AD crystal and BBR
    y, h, col = (1.1*helix_a3_pt1[[1,2]].max()) + 2, 2, 'k'
    plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
    plt.text((x1+x2)*0.5, y+h, "*" , ha='center', va='bottom', color=col)
if p_a3_alt < 0.01 and p_a3_alt > 0.001:
    x1, x2 = 2, 3 #Columns for AD crystal and BBR
    y, h, col = (1.1*helix_a3_pt1[[1,2]].max()) + 2, 2, 'k'
    plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
    plt.text((x1+x2)*0.5, y+h, "**" , ha='center', va='bottom', color=col)
if p_a3_alt < 0.001:
    x1, x2 = 2, 3 #Columns for AD crystal and BBR
    y, h, col = (1.1*helix_a3_pt1[[1,2]].max()) + 2, 2, 'k'
    plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
    plt.text((x1+x2)*0.5, y+h, "***" , ha='center', va='bottom', color=col)
fig.savefig('a3_pt1_helix_inter_cmpr.png')
plt.close(fig)

#Plot Bar graph comparing averages for 2nd half of a3 helix
fig = plt.figure()
plt.title('Interactions with Second Half of a3 Helix')    
plt.ylabel('Mean Number of Interactions')
plt.bar(num, helix_a3_pt2[:], color = ['green', 'green', 'blue'], width=0.8)
plt.errorbar(num, helix_a3_pt2[:], yerr=err_a3_pt2[:], fmt='o', color='black')
plt.xticks(num, label, fontsize=8)
if p_a3_pt2 < 0.05 and p_a3_pt2 > 0.01:
    x1, x2 = 1, 3 #Columns for AD crystal and BBR
    y, h, col = (1.1*helix_a3_pt2[[0,2]].max()) + 2, 2, 'k'
    plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
    plt.text((x1+x2)*0.5, y+h, "*" , ha='center', va='bottom', color=col)
if p_a3_pt2 < 0.01 and p_a3_pt2 > 0.001:
    x1, x2 = 1, 3 #Columns for AD crystal and BBR
    y, h, col = (1.1*helix_a3_pt2[[0,2]].max()) + 2, 2, 'k'
    plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
    plt.text((x1+x2)*0.5, y+h, "**" , ha='center', va='bottom', color=col)
if p_a3_pt2 < 0.001:
    x1, x2 = 1, 3 #Columns for AD crystal and BBR
    y, h, col = (1.1*helix_a3_pt2[[0,2]].max()) + 2, 2, 'k'
    plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
    plt.text((x1+x2)*0.5, y+h, "***" , ha='center', va='bottom', color=col)
if p_a3_pt2_alt < 0.05 and p_a3_pt2_alt > 0.01:
    x1, x2 = 2, 3 #Columns for AD crystal and BBR
    y, h, col = (1.1*helix_a3_pt2[[1,2]].max()) + 2, 2, 'k'
    plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
    plt.text((x1+x2)*0.5, y+h, "*" , ha='center', va='bottom', color=col)
if p_a3_pt2_alt < 0.01 and p_a3_pt2_alt > 0.001:
    x1, x2 = 2, 3 #Columns for AD crystal and BBR
    y, h, col = (1.1*helix_a3_pt2[[1,2]].max()) + 2, 2, 'k'
    plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
    plt.text((x1+x2)*0.5, y+h, "**" , ha='center', va='bottom', color=col)
if p_a3_pt2_alt < 0.001:
    x1, x2 = 2, 3 #Columns for AD crystal and BBR
    y, h, col = (1.1*helix_a3_pt2[[1,2]].max()) + 2, 2, 'k'
    plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
    plt.text((x1+x2)*0.5, y+h, "***" , ha='center', va='bottom', color=col)
fig.savefig('a3_pt2_helix_inter_cmpr.png')
plt.close(fig)

#Plot Bar graph comparing averages for a4 helix
fig = plt.figure()
plt.title('Interactions with First Half of a4 Helix')    
plt.ylabel('Mean Number of Interactions')
plt.bar(num, helix_a4_pt1[:], color = ['green', 'green', 'blue'], width=0.8)
plt.errorbar(num, helix_a4_pt1[:], yerr=err_a4_pt1[:], fmt='o', color='black')
plt.xticks(num, label, fontsize=8)
if p_a4_pt1 < 0.05 and p_a4_pt1 > 0.01:
    x1, x2 = 1, 3 #Columns for AD crystal and BBR
    y, h, col = (1.1*helix_a4_pt1[[0,2]].max()) + 2, 2, 'k'
    plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
    plt.text((x1+x2)*0.5, y+h, "*" , ha='center', va='bottom', color=col)
if p_a4_pt1 < 0.01 and p_a4_pt1 > 0.001:
    x1, x2 = 1, 3 #Columns for AD crystal and BBR
    y, h, col = (1.1*helix_a4_pt1[[0,2]].max()) + 2, 2, 'k'
    plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
    plt.text((x1+x2)*0.5, y+h, "**" , ha='center', va='bottom', color=col)
if p_a4_pt1 < 0.001:
    x1, x2 = 1, 3 #Columns for AD crystal and BBR
    y, h, col = (1.1*helix_a4_pt1[[0,2]].max()) + 2, 2, 'k'
    plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
    plt.text((x1+x2)*0.5, y+h, "***" , ha='center', va='bottom', color=col)
if p_a4_alt < 0.05 and p_a4_alt > 0.01:
    x1, x2 = 2, 3 #Columns for AD crystal and BBR
    y, h, col = (1.1*helix_a4_pt1[[1,2]].max()) + 2, 2, 'k'
    plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
    plt.text((x1+x2)*0.5, y+h, "*" , ha='center', va='bottom', color=col)
if p_a4_alt < 0.01 and p_a4_alt > 0.001:
    x1, x2 = 2, 3 #Columns for AD crystal and BBR
    y, h, col = (1.1*helix_a4_pt1[[1,2]].max()) + 2, 2, 'k'
    plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
    plt.text((x1+x2)*0.5, y+h, "**" , ha='center', va='bottom', color=col)
if p_a4_alt < 0.001:
    x1, x2 = 2, 3 #Columns for AD crystal and BBR
    y, h, col = (1.1*helix_a4_pt1[[1,2]].max()) + 2, 2, 'k'
    plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
    plt.text((x1+x2)*0.5, y+h, "***" , ha='center', va='bottom', color=col)
fig.savefig('a4_pt1_helix_inter_cmpr.png')
plt.close(fig)

#Plot Bar graph comparing averages for a4 helix
fig = plt.figure()
plt.title('Interactions with Second Half of a4 Helix')    
plt.ylabel('Mean Number of Interactions')
plt.bar(num, helix_a4_pt2[:], color = ['green', 'green', 'blue'], width=0.8)
plt.errorbar(num, helix_a4_pt2[:], yerr=err_a4_pt2[:], fmt='o', color='black')
plt.xticks(num, label, fontsize=8)
if p_a4_pt2 < 0.05 and p_a4_pt2 > 0.01:
    x1, x2 = 1, 3 #Columns for AD crystal and BBR
    y, h, col = (1.1*helix_a4_pt2[[0,2]].max()) + 2, 2, 'k'
    plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
    plt.text((x1+x2)*0.5, y+h, "*" , ha='center', va='bottom', color=col)
if p_a4_pt2 < 0.01 and p_a4_pt2 > 0.001:
    x1, x2 = 1, 3 #Columns for AD crystal and BBR
    y, h, col = (1.1*helix_a4_pt2[[0,2]].max()) + 2, 2, 'k'
    plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
    plt.text((x1+x2)*0.5, y+h, "**" , ha='center', va='bottom', color=col)
if p_a4_pt2 < 0.001:
    x1, x2 = 1, 3 #Columns for AD crystal and BBR
    y, h, col = (1.1*helix_a4_pt2[[0,2]].max()) + 2, 2, 'k'
    plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
    plt.text((x1+x2)*0.5, y+h, "***" , ha='center', va='bottom', color=col)
if p_a4_pt2_alt < 0.05 and p_a4_pt2_alt > 0.01:
    x1, x2 = 2, 3 #Columns for AD crystal and BBR
    y, h, col = (1.1*helix_a4_pt2[[1,2]].max()) + 2, 2, 'k'
    plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
    plt.text((x1+x2)*0.5, y+h, "*" , ha='center', va='bottom', color=col)
if p_a4_pt2_alt < 0.01 and p_a4_pt2_alt > 0.001:
    x1, x2 = 2, 3 #Columns for AD crystal and BBR
    y, h, col = (1.1*helix_a4_pt2[[1,2]].max()) + 2, 2, 'k'
    plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
    plt.text((x1+x2)*0.5, y+h, "**" , ha='center', va='bottom', color=col)
if p_a4_pt2_alt < 0.001:
    x1, x2 = 2, 3 #Columns for AD crystal and BBR
    y, h, col = (1.1*helix_a4_pt2[[1,2]].max()) + 2, 2, 'k'
    plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
    plt.text((x1+x2)*0.5, y+h, "***" , ha='center', va='bottom', color=col)
fig.savefig('a4_pt2_helix_inter_cmpr.png')
plt.close(fig)

#Plot Bar graph comparing averages for a6 helix
fig = plt.figure()
plt.title('Interactions with First Helf of a6 Helix')  
plt.ylabel('Mean Number of Interactions')
plt.bar(num, helix_a6_pt1[:], color = ['green', 'green', 'blue'], width=0.8)
plt.errorbar(num, helix_a6_pt1[:], yerr=err_a6_pt1[:], fmt='o', color='black')
plt.xticks(num, label, fontsize=8)
if p_a6_pt1 < 0.05 and p_a6_pt1 > 0.01:
    x1, x2 = 1, 3 #Columns for AD crystal and BBR
    y, h, col = (1.1*helix_a6_pt1[[0,2]].max()) + 2, 2, 'k'
    plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
    plt.text((x1+x2)*0.5, y+h, "*" , ha='center', va='bottom', color=col)
if p_a6_pt1 < 0.01 and p_a6_pt1 > 0.001:
    x1, x2 = 1, 3 #Columns for AD crystal and BBR
    y, h, col = (1.1*helix_a6_pt1[[0,2]].max()) + 2, 2, 'k'
    plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
    plt.text((x1+x2)*0.5, y+h, "**" , ha='center', va='bottom', color=col)
if p_a6_pt1 < 0.001:
    x1, x2 = 1, 3 #Columns for AD crystal and BBR
    y, h, col = (1.1*helix_a6_pt1[[0,2]].max()) + 2, 2, 'k'
    plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
    plt.text((x1+x2)*0.5, y+h, "***" , ha='center', va='bottom', color=col)
if p_a6_alt < 0.05 and p_a6_alt > 0.01:
    x1, x2 = 2, 3 #Columns for AD crystal and BBR
    y, h, col = (1.1*helix_a6_pt1[[1,2]].max()) + 2, 2, 'k'
    plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
    plt.text((x1+x2)*0.5, y+h, "*" , ha='center', va='bottom', color=col)
if p_a6_alt < 0.01 and p_a6_alt > 0.001:
    x1, x2 = 2, 3 #Columns for AD crystal and BBR
    y, h, col = (1.1*helix_a6_pt1[[1,2]].max()) + 2, 2, 'k'
    plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
    plt.text((x1+x2)*0.5, y+h, "**" , ha='center', va='bottom', color=col)
if p_a6_alt < 0.001:
    x1, x2 = 2, 3 #Columns for AD crystal and BBR
    y, h, col = (1.1*helix_a6_pt1[[1,2]].max()) + 2, 2, 'k'
    plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
    plt.text((x1+x2)*0.5, y+h, "***" , ha='center', va='bottom', color=col)
fig.savefig('a6_pt1_helix_inter_cmpr.png')
plt.close(fig)

#Plot Bar graph comparing averages for a6 helix
fig = plt.figure()
plt.title('Interactions with Second Helf of a6 Helix')  
plt.ylabel('Mean Number of Interactions')
plt.bar(num, helix_a6_pt2[:], color = ['green', 'green', 'blue'], width=0.8)
plt.errorbar(num, helix_a6_pt2[:], yerr=err_a6_pt2[:], fmt='o', color='black')
plt.xticks(num, label, fontsize=8)
if p_a6_pt2 < 0.05 and p_a6_pt2 > 0.01:
    x1, x2 = 1, 3 #Columns for AD crystal and BBR
    y, h, col = (1.1*helix_a6_pt2[[0,2]].max()) + 2, 2, 'k'
    plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
    plt.text((x1+x2)*0.5, y+h, "*" , ha='center', va='bottom', color=col)
if p_a6_pt2 < 0.01 and p_a6_pt2 > 0.001:
    x1, x2 = 1, 3 #Columns for AD crystal and BBR
    y, h, col = (1.1*helix_a6_pt2[[0,2]].max()) + 2, 2, 'k'
    plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
    plt.text((x1+x2)*0.5, y+h, "**" , ha='center', va='bottom', color=col)
if p_a6_pt2 < 0.001:
    x1, x2 = 1, 3 #Columns for AD crystal and BBR
    y, h, col = (1.1*helix_a6_pt2[[0,2]].max()) + 2, 2, 'k'
    plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
    plt.text((x1+x2)*0.5, y+h, "***" , ha='center', va='bottom', color=col)
if p_a6_pt2_alt < 0.05 and p_a6_pt2_alt > 0.01:
    x1, x2 = 2, 3 #Columns for AD crystal and BBR
    y, h, col = (1.1*helix_a6_pt2[[1,2]].max()) + 2, 2, 'k'
    plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
    plt.text((x1+x2)*0.5, y+h, "*" , ha='center', va='bottom', color=col)
if p_a6_pt2_alt < 0.01 and p_a6_pt2_alt > 0.001:
    x1, x2 = 2, 3 #Columns for AD crystal and BBR
    y, h, col = (1.1*helix_a6_pt2[[1,2]].max()) + 2, 2, 'k'
    plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
    plt.text((x1+x2)*0.5, y+h, "**" , ha='center', va='bottom', color=col)
if p_a6_pt2_alt < 0.001:
    x1, x2 = 2, 3 #Columns for AD crystal and BBR
    y, h, col = (1.1*helix_a6_pt2[[1,2]].max()) + 2, 2, 'k'
    plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
    plt.text((x1+x2)*0.5, y+h, "***" , ha='center', va='bottom', color=col)
fig.savefig('a6_pt2_helix_inter_cmpr.png')
plt.close(fig)

#Plot Bar graph comparing averages for a7 helix
fig = plt.figure()
plt.title('Interactions with First Half of a7 Helix')    
plt.ylabel('Mean Number of Interactions')
plt.bar(num, helix_a7_pt1[:], color = ['green', 'green', 'blue'], width=0.8)
plt.errorbar(num, helix_a7_pt1[:], yerr=err_a7_pt1[:], fmt='o', color='black')
plt.xticks(num, label, fontsize=8)
if p_a7_pt1 < 0.05 and p_a7_pt1 > 0.01:
    x1, x2 = 1, 3 #Columns for AD crystal and BBR
    y, h, col = (1.1*helix_a7_pt1[[0,2]].max()) + 2, 2, 'k'
    plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
    plt.text((x1+x2)*0.5, y+h, "*" , ha='center', va='bottom', color=col)
if p_a7_pt1 < 0.01 and p_a7_pt1 > 0.001:
    x1, x2 = 1, 3 #Columns for AD crystal and BBR
    y, h, col = (1.1*helix_pt1[[0,2]].max()) + 2, 2, 'k'
    plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
    plt.text((x1+x2)*0.5, y+h, "**" , ha='center', va='bottom', color=col)
if p_a7_pt1 < 0.001:
    x1, x2 = 1, 3 #Columns for AD crystal and BBR
    y, h, col = (1.1*helix_a7_pt1[[0,2]].max()) + 2, 2, 'k'
    plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
    plt.text((x1+x2)*0.5, y+h, "***" , ha='center', va='bottom', color=col)
if p_a7_alt < 0.05 and p_a7_alt > 0.01:
    x1, x2 = 2, 3 #Columns for AD crystal and BBR
    y, h, col = (1.1*helix_a7_pt1[[1,2]].max()) + 2, 2, 'k'
    plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
    plt.text((x1+x2)*0.5, y+h, "*" , ha='center', va='bottom', color=col)
if p_a7_alt < 0.01 and p_a7_alt > 0.001:
    x1, x2 = 2, 3 #Columns for AD crystal and BBR
    y, h, col = (1.1*helix_a7_pt1[[1,2]].max()) + 2, 2, 'k'
    plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
    plt.text((x1+x2)*0.5, y+h, "**" , ha='center', va='bottom', color=col)
if p_a7_alt < 0.001:
    x1, x2 = 2, 3 #Columns for AD crystal and BBR
    y, h, col = (1.1*helix_a7_pt1[[1,2]].max()) + 2, 2, 'k'
    plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
    plt.text((x1+x2)*0.5, y+h, "***" , ha='center', va='bottom', color=col)
fig.savefig('a7_pt1_helix_inter_cmpr.png')
plt.close(fig)

#Plot Bar graph comparing averages for a7 helix
fig = plt.figure()
plt.title('Interactions with Second Half of a7 Helix')    
plt.ylabel('Mean Number of Interactions')
plt.bar(num, helix_a7_pt2[:], color = ['green', 'green', 'blue'], width=0.8)
plt.errorbar(num, helix_a7_pt2[:], yerr=err_a7_pt2[:], fmt='o', color='black')
plt.xticks(num, label, fontsize=8)
if p_a7_pt2 < 0.05 and p_a7_pt2 > 0.01:
    x1, x2 = 1, 3 #Columns for AD crystal and BBR
    y, h, col = (1.1*helix_a7_pt2[[0,2]].max()) + 2, 2, 'k'
    plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
    plt.text((x1+x2)*0.5, y+h, "*" , ha='center', va='bottom', color=col)
if p_a7_pt2 < 0.01 and p_a7_pt2 > 0.001:
    x1, x2 = 1, 3 #Columns for AD crystal and BBR
    y, h, col = (1.1*helix_a7_pt2[[0,2]].max()) + 2, 2, 'k'
    plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
    plt.text((x1+x2)*0.5, y+h, "**" , ha='center', va='bottom', color=col)
if p_a7_pt2 < 0.001:
    x1, x2 = 1, 3 #Columns for AD crystal and BBR
    y, h, col = (1.1*helix_a7_pt2[[0,2]].max()) + 2, 2, 'k'
    plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
    plt.text((x1+x2)*0.5, y+h, "***" , ha='center', va='bottom', color=col)
if p_a7_pt2_alt < 0.05 and p_a7_pt2_alt > 0.01:
    x1, x2 = 2, 3 #Columns for AD crystal and BBR
    y, h, col = (1.1*helix_a7_pt2[[1,2]].max()) + 2, 2, 'k'
    plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
    plt.text((x1+x2)*0.5, y+h, "*" , ha='center', va='bottom', color=col)
if p_a7_pt2_alt < 0.01 and p_a7_pt2_alt > 0.001:
    x1, x2 = 2, 3 #Columns for AD crystal and BBR
    y, h, col = (1.1*helix_a7_pt2[[1,2]].max()) + 2, 2, 'k'
    plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
    plt.text((x1+x2)*0.5, y+h, "**" , ha='center', va='bottom', color=col)
if p_a7_pt2_alt < 0.001:
    x1, x2 = 2, 3 #Columns for AD crystal and BBR
    y, h, col = (1.1*helix_a7_pt2[[1,2]].max()) + 2, 2, 'k'
    plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
    plt.text((x1+x2)*0.5, y+h, "***" , ha='center', va='bottom', color=col)
fig.savefig('a7_pt2_helix_inter_cmpr.png')
plt.close(fig)

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

