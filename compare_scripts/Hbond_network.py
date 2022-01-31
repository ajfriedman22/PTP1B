#!/ usr / bin / env python

import mdtraj as md
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt

#Open files for percent of time each h-bond is formed
File_a7 = open('../../rebuild_a7/Hbond_per_single.txt', 'r').readlines()
File_a7_AD = open('../../AD_rebuild_a7/Hbond_per_single.txt', 'r').readlines()
File_apo = open('../../Apo/Hbond_per_single.txt', 'r').readlines()
File_AD = open('../../AD/Hbond_per_single.txt', 'r').readlines()
File_dis7 = open('../../Apo_dis/config7/Hbond_per_single.txt', 'r').readlines()
File_dis7_AD = open('../../AD_dis/config7/Hbond_per_single.txt', 'r').readlines()
File_dis9 = open('../../Apo_dis/config9/Hbond_per_single.txt', 'r').readlines()
File_dis9_AD = open('../../AD_dis/config9/Hbond_per_single.txt', 'r').readlines()
File_dis11 = open('../../Apo_dis/config11/Hbond_per_single.txt', 'r').readlines()
File_dis11_AD = open('../../AD_dis/config11/Hbond_per_single.txt', 'r').readlines()
File_1sug = open('../../1sug/Hbond_per_single.txt', 'r').readlines()
File_1sug_AD = open('../../1sug_AD/Hbond_per_single.txt', 'r').readlines()
File_1sug_alt_AD = open('../../1sug_AD_dis_alt/run_1/Hbond_per_single.txt', 'r').readlines()
File_1sug_alt2_AD = open('../../1sug_AD_dis_alt/run_2/Hbond_per_single.txt', 'r').readlines()
File_1sug_na7 = open('../../1sug_no_a7/Hbond_per_single.txt', 'r').readlines()
File_1sug_na7_AD = open('../../1sug_no_a7_AD/Hbond_per_single.txt', 'r').readlines()
File_1sug_dis7 = open('../../1sug_dis/config7/Hbond_per_single.txt', 'r').readlines()
File_1sug_dis7_AD = open('../../1sug_dis_AD/config7/Hbond_per_single.txt', 'r').readlines()
File_1sug_dis9 = open('../../1sug_dis/config9/Hbond_per_single.txt', 'r').readlines()
File_1sug_dis9_AD = open('../../1sug_dis_AD/config9/Hbond_per_single.txt', 'r').readlines()
File_1sug_dis11 = open('../../1sug_dis/config11/Hbond_per_single.txt', 'r').readlines()
File_1sug_dis11_AD = open('../../1sug_dis_AD/config11/Hbond_per_single.txt', 'r').readlines()
File_BBR_a7 = open('../../BBR_a7/Hbond_per_single.txt', 'r').readlines()
File_BBR_dis = open('../../BBR_dis/config9/Hbond_per_single.txt', 'r').readlines()

#Open File for bond names
Hbond_name = open('Hbond_uncommon.txt', 'r').readlines()

#Output Files for significant h-bonds
file_close = open('./Network/closed/Hbonds_list.txt', 'w')
file_open = open('./Network/open/Hbonds_list.txt', 'w')

#Output file for significant h-bonds
file_p = open('./Network/Hbond_p.txt', 'w')
file_p2 = open('./Network/Hbond_AD_p.txt', 'w')

#Make arrays of percentages the bonds are formed
per_a7, per_a7_AD, per_apo, per_AD, per_1sug, per_1sug_AD, per_1sug_na7, per_1sug_na7_AD, per_1sug_alt_AD, per_1sug_alt2_AD = [],[],[],[],[],[],[],[],[],[]
per_1sug_dis7, per_1sug_dis7_AD, per_1sug_dis9, per_1sug_dis9_AD, per_1sug_dis11, per_1sug_dis11_AD = [],[],[],[],[],[]
per_dis7, per_dis7_AD, per_dis9, per_dis9_AD, per_dis11, per_dis11_AD = [],[],[],[],[],[]
per_BBR_a7, per_BBR_dis = [],[]

for i in range(len(File_a7)):
    if i%2 != 0:
        per_a7.append(float(File_a7[i]))
        per_a7_AD.append(float(File_a7_AD[i]))
        per_apo.append(float(File_apo[i]))
        per_AD.append(float(File_AD[i]))
        per_dis7.append(float(File_dis7[i]))
        per_dis7_AD.append(float(File_dis7_AD[i]))
        per_dis9.append(float(File_dis9[i]))
        per_dis9_AD.append(float(File_dis9_AD[i]))
        per_dis11.append(float(File_dis11[i]))
        per_dis11_AD.append(float(File_dis11_AD[i]))
        per_1sug.append(float(File_1sug[i]))
        per_1sug_AD.append(float(File_1sug_AD[i]))
        per_1sug_alt_AD.append(float(File_1sug_alt_AD[i]))
        per_1sug_alt2_AD.append(float(File_1sug_alt2_AD[i]))
        per_1sug_na7.append(float(File_1sug_na7[i]))
        per_1sug_na7_AD.append(float(File_1sug_na7_AD[i]))
        per_1sug_dis7.append(float(File_1sug_dis7[i]))
        per_1sug_dis7_AD.append(float(File_1sug_dis7_AD[i]))
        per_1sug_dis9.append(float(File_1sug_dis9[i]))
        per_1sug_dis9_AD.append(float(File_1sug_dis9_AD[i]))
        per_1sug_dis11.append(float(File_1sug_dis11[i]))
        per_1sug_dis11_AD.append(float(File_1sug_dis11_AD[i]))
        per_BBR_a7.append(float(File_BBR_a7[i]))
        per_BBR_dis.append(float(File_BBR_dis[i]))

#Make Array for percentages of each bond
avg = np.zeros([len(per_apo), 9])
err = np.zeros([len(per_apo), 9])
for i in range(len(per_apo)):
    #seperate the open and closed trajectories
    open_na7 = [per_apo[i], per_AD[i]]
    open_a7 = [per_a7[i], per_1sug_dis7[i], per_1sug_dis9[i], per_a7_AD[i], per_1sug_AD[i], per_1sug_dis9_AD[i], per_dis9[i]]
    open_apo_a7 = [per_a7[i], per_1sug_dis7[i], per_1sug_dis9[i], per_dis9[i]]
    closed_na7 = [per_1sug_na7_AD[i]]
    closed_a7 = [per_1sug[i], per_1sug_dis11[i], per_dis11[i], per_dis7_AD[i], per_dis9_AD[i], per_dis11_AD[i]]
    closed_apo_a7 = [per_1sug[i], per_1sug_dis11[i], per_dis11[i]]
    bw = [per_1sug_na7[i], per_1sug_dis7_AD[i], per_dis7[i]]
    open_AD_crys = [per_1sug_dis11_AD[i], per_1sug_alt_AD[i], per_1sug_alt2_AD[i]]
    open_BBR_crys = [per_BBR_a7[i], per_BBR_dis[i]]

    #Determine average for open vs closed structures
    avg[i][0] = np.mean(open_na7) #Open + no a7
    avg[i][1] = np.mean(closed_na7) #Closed + no a7
    avg[i][2] = np.mean(open_a7) #Open + a7
    avg[i][3] = np.mean(closed_a7) #closed + a7
    avg[i][4] = np.mean(open_apo_a7) # open + apo + a7
    avg[i][5] = np.mean(closed_apo_a7) #closed + apo + a7
    avg[i][6] = np.mean(bw) #b/w open and closed
    avg[i][7] = np.mean(open_AD_crys) #Open + AD bound crys or alt 1
    avg[i][8] = np.mean(open_BBR_crys) #Open + BBR bound crys or alt 1
    
    #Error for the averages
    err[i][0] = stats.sem(open_na7, nan_policy = 'raise') #Open + no a7
    err[i][1] = 0 #Closed + no a7
    err[i][2] = stats.sem(open_a7, nan_policy = 'raise') #Open + a7
    err[i][3] = stats.sem(closed_a7, nan_policy = 'raise') #closed + a7
    err[i][4] = stats.sem(open_apo_a7, nan_policy = 'raise') #Open + apo + a7
    err[i][5] = stats.sem(closed_apo_a7, nan_policy = 'raise') #closed + apo + a7
    err[i][6] = stats.sem(bw, nan_policy = 'raise') #b/w open and closed
    err[i][7] = stats.sem(open_AD_crys) #Open + AD bound crys or alt 1
    err[i][8] = stats.sem(open_BBR_crys) #Open + BBR bound crys or alt 1

    #Calculate the difference b/w scenerios
    diff = avg[i][5] - avg[i][4] #Closed apo - Open apo
   
    #Determine the p-value for open vs closed with a7
    st, p = stats.ttest_ind(open_apo_a7, closed_apo_a7, equal_var = False) #Welch's t-test b/w open apo and closed apo
    st2, p2 = stats.ttest_ind(open_a7, open_AD_crys, equal_var = False) #Welch's t-test b/w open and AD
    st3, p3 = stats.ttest_ind(closed_a7, open_AD_crys, equal_var = False) #Welch's t-test b/w closed and AD
    st4, p4 = stats.ttest_ind(open_a7, open_BBR_crys, equal_var = False) #Welch's t-test b/w open and BBR
    st5, p5 = stats.ttest_ind(closed_a7, open_BBR_crys, equal_var = False) #Welch's t-test b/w closed and BBR
    st6, p6 = stats.ttest_ind(open_BBR_crys, open_AD_crys, equal_var = False) #Welch's t-test b/w BBR and AD

    #Only print graphs with significant differences
    num = [5, 10, 15]
    Method = ['Open Apo + a7', 'Closed Apo + a7', 'b/w Open + Closed']
    if diff > 0 and p <= 0.05:
        #Write the hbonds to file
        file_close.write(str(Hbond_name[i]))
        file_p.write(str(Hbond_name[i]) + str(p) + '\n')
        #Plot Bar graph comparing averages
        fig = plt.figure()
        plt.title('Time Hbond ' + str(Hbond_name[i]) + ' is formed')    
        plt.ylabel('% Time Hbond Formed')
        plt.bar(num,avg[i][[4, 5, 6]], color = ['blue', 'red', 'purple'], width=4.5)
        plt.errorbar(num, avg[i][[4, 5, 6]], yerr=err[i][[4, 5, 6]], fmt='o', color='black')
        plt.xticks(num, Method, fontsize=8)
        if p < 0.05 and p > 0.01:
            x1, x2 = 5, 10 #Columns for Open Apo + a7 and Closed Apo + a7
            y, h, col = (1.1*avg[i][[4, 5]].max()) + 2, 2, 'k'
            plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
            plt.text((x1+x2)*0.5, y+h, "*" , ha='center', va='bottom', color=col)
        if p < 0.01 and p > 0.001:
            x1, x2 = 5, 10 #Columns for Open + a7 and Closed + a7
            y, h, col = (1.1*avg[i][[4, 5]].max()) + 2, 2, 'k'
            plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
            plt.text((x1+x2)*0.5, y+h, "**" , ha='center', va='bottom', color=col)
        if p < 0.001:
            x1, x2 = 5, 10 #Columns for Open + a7 and Closed + a7
            y, h, col = (1.1*avg[i][[4, 5]].max()) + 2, 2, 'k'
            plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
            plt.text((x1+x2)*0.5, y+h, "***" , ha='center', va='bottom', color=col)
        fig.savefig('./Network/closed/' + str(Hbond_name[i]) + '_group_percent.png')
        plt.close(fig)
    if diff < 0 and p <= 0.05:
        #Write the hbonds to file
        file_open.write(str(Hbond_name[i]))
        file_p.write(str(Hbond_name[i]) + str(p) + '\n')
        #Plot Bar graph comparing averages
        fig = plt.figure()
        plt.title('Time Hbond ' + str(Hbond_name[i]) + ' is formed')    
        plt.ylabel('% Time Hbond Formed')
        plt.bar(num,avg[i][[4, 5, 6]], color = ['blue', 'red', 'purple'], width=4.5)
        plt.errorbar(num, avg[i][[4, 5, 6]], yerr=err[i][[4, 5, 6]], fmt='o', color='black')
        plt.xticks(num, Method, fontsize=8)
        if p < 0.05 and p > 0.01:
            x1, x2 = 5, 10 #Columns for Open + a7 and Closed + a7
            y, h, col = (1.1*avg[i][[4, 5]].max()) + 2, 2, 'k'
            plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
            plt.text((x1+x2)*0.5, y+h, "*" , ha='center', va='bottom', color=col)
        if p < 0.01 and p > 0.001:
            x1, x2 = 5, 10 #Columns for Open + a7 and Closed + a7
            y, h, col = (1.1*avg[i][[4, 5]].max()) + 2, 2, 'k'
            plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
            plt.text((x1+x2)*0.5, y+h, "**" , ha='center', va='bottom', color=col)
        if p < 0.001:
            x1, x2 = 5, 10 #Columns for Open + a7 and Closed + a7
            y, h, col = (1.1*avg[i][[4, 5]].max()) + 2, 2, 'k'
            plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
            plt.text((x1+x2)*0.5, y+h, "***" , ha='center', va='bottom', color=col)
        fig.savefig('./Network/open/' + str(Hbond_name[i]) + '_group_percent.png')
        plt.close(fig)
    num = [5, 10, 15]
    Method = ['Open', 'Closed', 'Open + AD Crys']
    if p3 <= 0.05 and p2 <= 0.05:
        #Write the hbonds to file
        file_p2.write(str(Hbond_name[i]) + str(p2) + ' ' + str(p3) + '\n')
        #Plot Bar graph comparing averages
        fig = plt.figure()
        plt.title('Time Hbond ' + str(Hbond_name[i]) + ' is formed')    
        plt.ylabel('% Time Hbond Formed')
        plt.bar(num,avg[i][[2,3,7]], color = ['blue', 'red', 'blue'], width=4.5)
        plt.errorbar(num, avg[i][[2,3,7]], yerr=err[i][[2,3,7]], fmt='o', color='black')
        plt.xticks(num, Method, fontsize=8)
        if p2 < 0.05 and p2 > 0.01:
            x1, x2 = 5, 15 #Columns for Open + a7 and AD
            y, h, col = (1.1*avg[i][[2,7]].max()) + 2, 2, 'k'
            plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
            plt.text((x1+x2)*0.5, y+h, "*" , ha='center', va='bottom', color=col)
        if p2 < 0.01 and p2 > 0.001:
            x1, x2 = 5, 15 #Columns for Open + a7 and Closed + a7
            y, h, col = (1.1*avg[i][[2,7]].max()) + 2, 2, 'k'
            plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
            plt.text((x1+x2)*0.5, y+h, "**" , ha='center', va='bottom', color=col)
        if p2 < 0.001:
            x1, x2 = 5, 15 #Columns for Open + a7 and Closed + a7
            y, h, col = (1.1*avg[i][[2,7]].max()) + 2, 2, 'k'
            plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
            plt.text((x1+x2)*0.5, y+h, "***" , ha='center', va='bottom', color=col)
        if p3 < 0.05 and p3 > 0.01:
            x1, x2 = 10, 15 #Columns for Open + a7 and Closed + a7
            y, h, col = (1.1*avg[i][[3,7]].max()) + 2, 2, 'k'
            plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
            plt.text((x1+x2)*0.5, y+h, "*" , ha='center', va='bottom', color=col)
        if p3 < 0.01 and p3 > 0.001:
            x1, x2 = 10, 15 #Columns for Open + a7 and Closed + a7
            y, h, col = (1.1*avg[i][[3,7]].max()) + 2, 2, 'k'
            plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
            plt.text((x1+x2)*0.5, y+h, "**" , ha='center', va='bottom', color=col)
        if p3 < 0.001:
            x1, x2 = 10, 15 #Columns for Open + a7 and Closed + a7
            y, h, col = (1.1*avg[i][[3,7]].max()) + 2, 2, 'k'
            plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
            plt.text((x1+x2)*0.5, y+h, "***" , ha='center', va='bottom', color=col)
        fig.savefig('./Network/AD/' + str(Hbond_name[i]) + '_group_percent.png')
        plt.close(fig)
    num = [5, 10, 15]
    Method = ['Open', 'Closed', 'Open + BBR Crys']
    if p4 <= 0.05 and p5 <= 0.05:
        #Write the hbonds to file
        file_p2.write(str(Hbond_name[i]) + str(p4) + ' ' + str(p5) + '\n')
        #Plot Bar graph comparing averages
        fig = plt.figure()
        plt.title('Time Hbond ' + str(Hbond_name[i]) + ' is formed')    
        plt.ylabel('% Time Hbond Formed')
        plt.bar(num,avg[i][[2,3,8]], color = ['blue', 'red', 'blue'], width=4.5)
        plt.errorbar(num, avg[i][[2,3,8]], yerr=err[i][[2,3,8]], fmt='o', color='black')
        plt.xticks(num, Method, fontsize=8)
        if p4 < 0.05 and p4 > 0.01:
            x1, x2 = 5, 15 #Columns for Open + a7 and BBR
            y, h, col = (1.1*avg[i][[2,8]].max()) + 2, 2, 'k'
            plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
            plt.text((x1+x2)*0.5, y+h, "*" , ha='center', va='bottom', color=col)
        if p4 < 0.01 and p4 > 0.001:
            x1, x2 = 5, 15 #Columns for Open + a7 and BBR
            y, h, col = (1.1*avg[i][[2,8]].max()) + 2, 2, 'k'
            plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
            plt.text((x1+x2)*0.5, y+h, "**" , ha='center', va='bottom', color=col)
        if p4 < 0.001:
            x1, x2 = 5, 15 #Columns for Open + a7 and BBR
            y, h, col = (1.1*avg[i][[2,8]].max()) + 2, 2, 'k'
            plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
            plt.text((x1+x2)*0.5, y+h, "***" , ha='center', va='bottom', color=col)
        if p5 < 0.05 and p5 > 0.01:
            x1, x2 = 10, 15 #Columns for Closed + a7 and BBR
            y, h, col = (1.1*avg[i][[3,8]].max()) + 2, 2, 'k'
            plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
            plt.text((x1+x2)*0.5, y+h, "*" , ha='center', va='bottom', color=col)
        if p5 < 0.01 and p5 > 0.001:
            x1, x2 = 10, 15 #Columns for Closed + a7 and BBR
            y, h, col = (1.1*avg[i][[3,8]].max()) + 2, 2, 'k'
            plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
            plt.text((x1+x2)*0.5, y+h, "**" , ha='center', va='bottom', color=col)
        if p5 < 0.001:
            x1, x2 = 10, 15 #Columns for Closed + a7 and BBR
            y, h, col = (1.1*avg[i][[3,8]].max()) + 2, 2, 'k'
            plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
            plt.text((x1+x2)*0.5, y+h, "***" , ha='center', va='bottom', color=col)
        fig.savefig('./Network/BBR/' + str(Hbond_name[i]) + '_group_percent.png')
        plt.close(fig)
    num = [5, 10, 15, 20]
    Method = ['Open', 'Closed', 'AD', 'BBR']
    if (p2 <= 0.05 and p3 <= 0.05) or (p4 <= 0.05 and p5 <= 0.05) and p6 <= 0.05:
        #Write the hbonds to file
        file_p2.write(str(Hbond_name[i]) + str(p2) + ' ' + str(p3) + ' ' + str(p4) + ' ' + str(p5) + ' ' + str(p6) + '\n')
        #Plot Bar graph comparing averages
        fig = plt.figure()
        plt.title('Time Hbond ' + str(Hbond_name[i]) + ' is formed')
        plt.ylabel('% Time Hbond Formed')
        plt.bar(num,avg[i][[4, 5, 7, 8]], color = ['blue', 'blue'], width=4.5)
        plt.errorbar(num, avg[i][[4, 5, 7, 8]], yerr=err[i][[4, 5, 7, 8]], fmt='o', color='black')
        plt.xticks(num, Method, fontsize=8)
        if p6 < 0.05 and p6 > 0.01:
            x1, x2 = 15, 20 #Columns for AD and BBR
            y, h, col = (1.1*avg[i][[7, 8]].max()) + 2, 2, 'k'
            plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
            plt.text((x1+x2)*0.5, y+h, "*" , ha='center', va='bottom', color=col)
        if p6 < 0.01 and p6 > 0.001:
            x1, x2 = 15, 20 #Columns for AD and BBR
            y, h, col = (1.1*avg[i][[7, 8]].max()) + 2, 2, 'k'
            plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
            plt.text((x1+x2)*0.5, y+h, "**" , ha='center', va='bottom', color=col)
        if p6 < 0.001:
            x1, x2 = 15, 20 #Columns for AD and BBR
            y, h, col = (1.1*avg[i][[7, 8]].max()) + 2, 2, 'k'
            plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
            plt.text((x1+x2)*0.5, y+h, "***" , ha='center', va='bottom', color=col)
        fig.savefig('./Network/AD_BBR_diff/' + str(Hbond_name[i]) + '_group_percent.png')
        plt.close(fig)

#Open files for percent of time each h-bond found in literature is formed
File_a7 = open('../../rebuild_a7/Hbond_lit_per.txt', 'r').readlines()
File_a7_AD = open('../../AD_rebuild_a7/Hbond_lit_per.txt', 'r').readlines()
File_dis7 = open('../../Apo_dis/config7/Hbond_lit_per.txt', 'r').readlines()
File_dis7_AD = open('../../AD_dis/config7/Hbond_lit_per.txt', 'r').readlines()
File_dis9 = open('../../Apo_dis/config9/Hbond_lit_per.txt', 'r').readlines()
File_dis9_AD = open('../../AD_dis/config9/Hbond_lit_per.txt', 'r').readlines()
File_dis11 = open('../../Apo_dis/config11/Hbond_lit_per.txt', 'r').readlines()
File_dis11_AD = open('../../AD_dis/config11/Hbond_lit_per.txt', 'r').readlines()
File_1sug = open('../../1sug/Hbond_lit_per.txt', 'r').readlines()
File_1sug_AD = open('../../1sug_AD/Hbond_lit_per.txt', 'r').readlines()
File_1sug_alt_AD = open('../../1sug_AD_dis_alt/run_1/Hbond_lit_per.txt', 'r').readlines()
File_1sug_alt2_AD = open('../../1sug_AD_dis_alt/run_2/Hbond_lit_per.txt', 'r').readlines()
File_1sug_dis7 = open('../../1sug_dis/config7/Hbond_lit_per.txt', 'r').readlines()
File_1sug_dis7_AD = open('../../1sug_dis_AD/config7/Hbond_lit_per.txt', 'r').readlines()
File_1sug_dis9 = open('../../1sug_dis/config9/Hbond_lit_per.txt', 'r').readlines()
File_1sug_dis9_AD = open('../../1sug_dis_AD/config9/Hbond_lit_per.txt', 'r').readlines()
File_1sug_dis11 = open('../../1sug_dis/config11/Hbond_lit_per.txt', 'r').readlines()
File_1sug_dis11_AD = open('../../1sug_dis_AD/config11/Hbond_lit_per.txt', 'r').readlines()
File_BBR_a7 = open('../../BBR_a7/Hbond_lit_per.txt', 'r').readlines()
File_BBR_dis = open('../../BBR_dis/config9/Hbond_lit_per.txt', 'r').readlines()

Hbond_name_lit = open('Hbond_lit.txt', 'r').readlines()

per_a7, per_a7_AD, per_apo, per_AD, per_1sug, per_1sug_AD, per_1sug_na7, per_1sug_na7_AD, per_1sug_alt_AD, per_1sug_alt2_AD = [],[],[],[],[],[],[],[],[],[]
per_1sug_dis7, per_1sug_dis7_AD, per_1sug_dis9, per_1sug_dis9_AD, per_1sug_dis11, per_1sug_dis11_AD = [],[],[],[],[],[]
per_dis7, per_dis7_AD, per_dis9, per_dis9_AD, per_dis11, per_dis11_AD = [],[],[],[],[],[]
per_BBR_a7, per_BBR_dis = [],[]

for i in range(len(Hbond_name_lit)):
    per_a7.append(float(File_a7[i]))
    per_a7_AD.append(float(File_a7_AD[i]))
    per_dis7.append(float(File_dis7[i]))
    per_dis7_AD.append(float(File_dis7_AD[i]))
    per_dis9.append(float(File_dis9[i]))
    per_dis9_AD.append(float(File_dis9_AD[i]))
    per_dis11.append(float(File_dis11[i]))
    per_dis11_AD.append(float(File_dis11_AD[i]))
    per_1sug.append(float(File_1sug[i]))
    per_1sug_AD.append(float(File_1sug_AD[i]))
    per_1sug_alt_AD.append(float(File_1sug_alt_AD[i]))
    per_1sug_alt2_AD.append(float(File_1sug_alt2_AD[i]))
    per_1sug_dis7.append(float(File_1sug_dis7[i]))
    per_1sug_dis7_AD.append(float(File_1sug_dis7_AD[i]))
    per_1sug_dis9.append(float(File_1sug_dis9[i]))
    per_1sug_dis9_AD.append(float(File_1sug_dis9_AD[i]))
    per_1sug_dis11.append(float(File_1sug_dis11[i]))
    per_1sug_dis11_AD.append(float(File_1sug_dis11_AD[i]))
    per_BBR_a7.append(float(File_BBR_a7[i]))
    per_BBR_dis.append(float(File_BBR_dis[i]))

avg = np.zeros([len(Hbond_name_lit), 7])
err = np.zeros([len(Hbond_name_lit), 7])
for i in range(len(Hbond_name_lit)):
    #seperate the open and closed trajectories
    open_a7 = [per_a7[i], per_1sug_dis7[i], per_1sug_dis9[i], per_a7_AD[i], per_1sug_AD[i], per_1sug_dis9_AD[i], per_dis9[i]]
    open_apo_a7 = [per_a7[i], per_1sug_dis7[i], per_1sug_dis9[i], per_dis9[i]]
    closed_a7 = [per_1sug[i], per_1sug_dis11[i], per_dis11[i], per_dis7_AD[i], per_dis9_AD[i], per_dis11_AD[i]]
    closed_apo_a7 = [per_1sug[i], per_1sug_dis11[i], per_dis11[i]]
    bw = [per_1sug_dis7_AD[i], per_dis7[i]]
    open_AD_crys = [per_1sug_dis11_AD[i], per_1sug_alt_AD[i], per_1sug_alt2_AD[i]]
    open_BBR_crys = [per_BBR_a7[i], per_BBR_dis[i]]

    #Determine average for open vs closed structures
    avg[i][0] = np.mean(open_a7) #Open + a7
    avg[i][1] = np.mean(closed_a7) #closed + a7
    avg[i][2] = np.mean(open_apo_a7) # open + apo + a7
    avg[i][3] = np.mean(closed_apo_a7) #closed + apo + a7
    avg[i][4] = np.mean(bw) #b/w open and closed
    avg[i][5] = np.mean(open_AD_crys) #Open + AD bound crys or alt 1
    avg[i][6] = np.mean(open_BBR_crys) #Open + BBR bound crys or alt 1
    
    #Error for the averages
    err[i][0] = stats.sem(open_a7, nan_policy = 'raise') #Open + a7
    err[i][1] = stats.sem(closed_a7, nan_policy = 'raise') #closed + a7
    err[i][2] = stats.sem(open_apo_a7, nan_policy = 'raise') #Open + apo + a7
    err[i][3] = stats.sem(closed_apo_a7, nan_policy = 'raise') #closed + apo + a7
    err[i][4] = stats.sem(bw, nan_policy = 'raise') #b/w open and closed
    err[i][5] = stats.sem(open_AD_crys) #Open + AD bound crys or alt 1
    err[i][6] = stats.sem(open_BBR_crys) #Open + BBR bound crys or alt 1

   
    #Determine the p-value for open vs closed with a7
    st, p = stats.ttest_ind(open_apo_a7, closed_apo_a7, equal_var = False) #Welch's t-test b/w open apo and closed apo
    st2, p2 = stats.ttest_ind(open_apo_a7, open_AD_crys, equal_var = False) #Welch's t-test b/w open and AD
    st3, p3 = stats.ttest_ind(closed_apo_a7, open_AD_crys, equal_var = False) #Welch's t-test b/w closed and AD
    st4, p4 = stats.ttest_ind(open_apo_a7, open_BBR_crys, equal_var = False) #Welch's t-test b/w open and BBR
    st5, p5 = stats.ttest_ind(closed_apo_a7, open_BBR_crys, equal_var = False) #Welch's t-test b/w closed and BBR
    st6, p6 = stats.ttest_ind(open_BBR_crys, open_AD_crys, equal_var = False) #Welch's t-test b/w BBR and AD

    #Print all lit graphs
    num = [5, 10, 15]
    Method = ['Open Apo', 'Closed Apo', 'b/w Open + Closed']
    #Plot Bar graph comparing averages
    fig = plt.figure()
    plt.title('Time Hbond ' + str(Hbond_name_lit[i]) + ' is formed')    
    plt.ylabel('% Time Hbond Formed')
    plt.bar(num,avg[i][[2, 3, 4]], color = ['blue', 'red', 'purple'], width=4.5)
    plt.errorbar(num, avg[i][[2, 3, 4]], yerr=err[i][[2, 3, 4]], fmt='o', color='black')
    plt.xticks(num, Method, fontsize=8)
    if p < 0.05 and p > 0.01:
        x1, x2 = 5, 10 #Columns for Open Apo + a7 and Closed Apo + a7
        y, h, col = (1.1*avg[i][[2, 3]].max()) + 2, 2, 'k'
        plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
        plt.text((x1+x2)*0.5, y+h, "*" , ha='center', va='bottom', color=col)
    if p < 0.01 and p > 0.001:
        x1, x2 = 5, 10 #Columns for Open + a7 and Closed + a7
        y, h, col = (1.1*avg[i][[2, 3]].max()) + 2, 2, 'k'
        plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
        plt.text((x1+x2)*0.5, y+h, "**" , ha='center', va='bottom', color=col)
    if p < 0.001:
        x1, x2 = 5, 10 #Columns for Open + a7 and Closed + a7
        y, h, col = (1.1*avg[i][[2, 3]].max()) + 2, 2, 'k'
        plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
        plt.text((x1+x2)*0.5, y+h, "***" , ha='center', va='bottom', color=col)
    if p < 0.05:
        fig.savefig('./Network/Lit/Significant/' + str(Hbond_name_lit[i]) + '_apo_percent.png')
    else:
        fig.savefig('./Network/Lit/Not_Significant/' + str(Hbond_name_lit[i]) + '_apo_percent.png')
    plt.close(fig)

    num = [5, 10, 15, 20]
    Method = ['Open Apo', 'Closed Apo', 'AD Open', 'BBR Open']
    #Plot Bar graph comparing averages
    fig = plt.figure()
    plt.title('Time Hbond ' + str(Hbond_name_lit[i]) + ' is formed')    
    plt.ylabel('% Time Hbond Formed')
    plt.bar(num,avg[i][[2, 3, 5, 6]], color = ['blue', 'red', 'blue', 'blue'], width=4.5)
    plt.errorbar(num, avg[i][[2, 3, 5, 6]], yerr=err[i][[2, 3, 5, 6]], fmt='o', color='black')
    plt.xticks(num, Method, fontsize=8)
    if p < 0.05 and p > 0.01:
        x1, x2 = 5, 10 #Columns for Open Apo + a7 and Closed Apo + a7
        y, h, col = (1.1*avg[i][[2, 3]].max()) + 2, 2, 'k'
        plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
        plt.text((x1+x2)*0.5, y+h, "*" , ha='center', va='bottom', color=col)
    if p < 0.01 and p > 0.001:
        x1, x2 = 5, 10 #Columns for Open + a7 and Closed + a7
        y, h, col = (1.1*avg[i][[2, 3]].max()) + 2, 2, 'k'
        plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
        plt.text((x1+x2)*0.5, y+h, "**" , ha='center', va='bottom', color=col)
    if p < 0.001:
        x1, x2 = 5, 10 #Columns for Open + a7 and Closed + a7
        y, h, col = (1.1*avg[i][[2, 3]].max()) + 2, 2, 'k'
        plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
        plt.text((x1+x2)*0.5, y+h, "***" , ha='center', va='bottom', color=col)
    if p3 < 0.05 and p3 > 0.01:
        x1, x2 = 10, 15 #Columns for Open Apo + a7 and Closed Apo + a7
        y, h, col = (1.1*avg[i][[3, 5]].max()) + 2, 2, 'k'
        plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
        plt.text((x1+x2)*0.5, y+h, "*" , ha='center', va='bottom', color=col)
    if p3 < 0.01 and p3 > 0.001:
        x1, x2 = 10, 15 #Columns for Open + a7 and Closed + a7
        y, h, col = (1.1*avg[i][[3, 5]].max()) + 2, 2, 'k'
        plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
        plt.text((x1+x2)*0.5, y+h, "**" , ha='center', va='bottom', color=col)
    if p3 < 0.001:
        x1, x2 = 5, 10 #Columns for Open + a7 and Closed + a7
        y, h, col = (1.1*avg[i][[3, 5]].max()) + 2, 2, 'k'
        plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
        plt.text((x1+x2)*0.5, y+h, "***" , ha='center', va='bottom', color=col)
    if p5 < 0.05 and p5 > 0.01:
        x1, x2 = 10, 20 #Columns for Open Apo + a7 and Closed Apo + a7
        y, h, col = (1.1*avg[i][[3, 6]].max()) + 2, 2, 'k'
        plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
        plt.text((x1+x2)*0.5, y+h, "*" , ha='center', va='bottom', color=col)
    if p5 < 0.01 and p5 > 0.001:
        x1, x2 = 10, 20 #Columns for Open + a7 and Closed + a7
        y, h, col = (1.1*avg[i][[3, 6]].max()) + 2, 2, 'k'
        plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
        plt.text((x1+x2)*0.5, y+h, "**" , ha='center', va='bottom', color=col)
    if p5 < 0.001:
        x1, x2 = 10, 20 #Columns for Open + a7 and Closed + a7
        y, h, col = (1.1*avg[i][[3, 6]].max()) + 2, 2, 'k'
        plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
        plt.text((x1+x2)*0.5, y+h, "***" , ha='center', va='bottom', color=col)
    if p < 0.05:
        fig.savefig('./Network/Lit/Significant/' + str(Hbond_name_lit[i]) + '_percent.png')
    else:
        fig.savefig('./Network/Lit/Not_Significant/' + str(Hbond_name_lit[i]) + '_percent.png')
    plt.close(fig)

