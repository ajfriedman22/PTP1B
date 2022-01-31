import numpy as np
from matplotlib import pyplot as plt
from scipy import stats

#Input Apo files
beta_Apo_open, L11_Apo_open, a3_Apo_open, a6_Apo_open, a7_Apo_open = [],[],[],[],[]
beta_Apo_close, L11_Apo_close, a3_Apo_close, a6_Apo_close, a7_Apo_close = [],[],[],[],[]

for i in open('beta_1sug2_rmsd.txt', 'r').readlines():
    beta_Apo_close.append(float(i)*10)
for i in open('beta_1sug3_rmsd.txt', 'r').readlines():
    beta_Apo_open.append(float(i)*10)
for i in open('beta_1sug_dis7_rmsd.txt', 'r').readlines():
    beta_Apo_open.append(float(i)*10)
for i in open('beta_1sug_dis9_rmsd.txt', 'r').readlines():
    beta_Apo_open.append(float(i)*10)
for i in open('beta_1sug_dis11_rmsd.txt', 'r').readlines():
    beta_Apo_close.append(float(i)*10)
for i in open('beta_a7_rmsd.txt', 'r').readlines():
    beta_Apo_open.append(float(i)*10)
for i in open('beta_dis9_rmsd.txt', 'r').readlines():
    beta_Apo_open.append(float(i)*10)
for i in open('beta_dis11_rmsd.txt', 'r').readlines():
    beta_Apo_close.append(float(i)*10)

for i in open('L11_1sug2_rmsd.txt', 'r').readlines():
    L11_Apo_close.append(float(i)*10)
for i in open('L11_1sug3_rmsd.txt', 'r').readlines():
    L11_Apo_open.append(float(i)*10)
for i in open('L11_1sug_dis7_rmsd.txt', 'r').readlines():
    L11_Apo_open.append(float(i)*10)
for i in open('L11_1sug_dis9_rmsd.txt', 'r').readlines():
    L11_Apo_open.append(float(i)*10)
for i in open('L11_1sug_dis11_rmsd.txt', 'r').readlines():
    L11_Apo_close.append(float(i)*10)
for i in open('L11_a7_rmsd.txt', 'r').readlines():
    L11_Apo_open.append(float(i)*10)
for i in open('L11_dis9_rmsd.txt', 'r').readlines():
    L11_Apo_open.append(float(i)*10)
for i in open('L11_dis11_rmsd.txt', 'r').readlines():
    L11_Apo_close.append(float(i)*10)

for i in open('a3_1sug2_rmsd.txt', 'r').readlines():
    a3_Apo_close.append(float(i)*10)
for i in open('a3_1sug3_rmsd.txt', 'r').readlines():
    a3_Apo_open.append(float(i)*10)
for i in open('a3_1sug_dis7_rmsd.txt', 'r').readlines():
    a3_Apo_open.append(float(i)*10)
for i in open('a3_1sug_dis9_rmsd.txt', 'r').readlines():
    a3_Apo_open.append(float(i)*10)
for i in open('a3_1sug_dis11_rmsd.txt', 'r').readlines():
    a3_Apo_close.append(float(i)*10)
for i in open('a3_a7_rmsd.txt', 'r').readlines():
    a3_Apo_open.append(float(i)*10)
for i in open('a3_dis9_rmsd.txt', 'r').readlines():
    a3_Apo_open.append(float(i)*10)
for i in open('a3_dis11_rmsd.txt', 'r').readlines():
    a3_Apo_close.append(float(i)*10)

for i in open('a6_1sug2_rmsd.txt', 'r').readlines():
    a6_Apo_close.append(float(i)*10)
for i in open('a6_1sug3_rmsd.txt', 'r').readlines():
    a6_Apo_open.append(float(i)*10)
for i in open('a6_1sug_dis7_rmsd.txt', 'r').readlines():
    a6_Apo_open.append(float(i)*10)
for i in open('a6_1sug_dis9_rmsd.txt', 'r').readlines():
    a6_Apo_open.append(float(i)*10)
for i in open('a6_1sug_dis11_rmsd.txt', 'r').readlines():
    a6_Apo_close.append(float(i)*10)
for i in open('a6_a7_rmsd.txt', 'r').readlines():
    a6_Apo_open.append(float(i)*10)
for i in open('a6_dis9_rmsd.txt', 'r').readlines():
    a6_Apo_open.append(float(i)*10)
for i in open('a6_dis11_rmsd.txt', 'r').readlines():
    a6_Apo_close.append(float(i)*10)

for i in open('a7_1sug2_rmsd.txt', 'r').readlines():
    a7_Apo_close.append(float(i)*10)
for i in open('a7_1sug3_rmsd.txt', 'r').readlines():
    a7_Apo_open.append(float(i)*10)
for i in open('a7_1sug_dis7_rmsd.txt', 'r').readlines():
    a7_Apo_open.append(float(i)*10)
for i in open('a7_1sug_dis9_rmsd.txt', 'r').readlines():
    a7_Apo_open.append(float(i)*10)
for i in open('a7_1sug_dis11_rmsd.txt', 'r').readlines():
    a7_Apo_close.append(float(i)*10)
for i in open('a7_a7_rmsd.txt', 'r').readlines():
    a7_Apo_open.append(float(i)*10)
for i in open('a7_dis9_rmsd.txt', 'r').readlines():
    a7_Apo_open.append(float(i)*10)
for i in open('a7_dis11_rmsd.txt', 'r').readlines():
    a7_Apo_close.append(float(i)*10)

#Input AD bound files
beta_AD_open, L11_AD_open, a3_AD_open, a6_AD_open, a7_AD_open = [],[],[],[],[]
beta_AD_close, L11_AD_close, a3_AD_close, a6_AD_close, a7_AD_close = [],[],[],[],[]

for i in open('beta_1sug_AD_dis11_rmsd.txt', 'r').readlines():
    beta_AD_open.append(float(i)*10)
for i in open('beta_1sug_AD_dis11_2_rmsd.txt', 'r').readlines():
    beta_AD_close.append(float(i)*10)
for i in open('beta_1sug_AD_alt1_rmsd.txt', 'r').readlines():
    beta_AD_open.append(float(i)*10)
for i in open('beta_1sug_AD_alt2_rmsd.txt', 'r').readlines():
    beta_AD_open.append(float(i)*10)

for i in open('L11_1sug_AD_dis11_rmsd.txt', 'r').readlines():
    L11_AD_open.append(float(i)*10)
for i in open('L11_1sug_AD_dis11_2_rmsd.txt', 'r').readlines():
    L11_AD_close.append(float(i)*10)
for i in open('L11_1sug_AD_alt1_rmsd.txt', 'r').readlines():
    L11_AD_open.append(float(i)*10)
for i in open('L11_1sug_AD_alt2_rmsd.txt', 'r').readlines():
    L11_AD_open.append(float(i)*10)

for i in open('a3_1sug_AD_dis11_rmsd.txt', 'r').readlines():
    a3_AD_open.append(float(i)*10)
for i in open('a3_1sug_AD_dis11_2_rmsd.txt', 'r').readlines():
    a3_AD_close.append(float(i)*10)
for i in open('a3_1sug_AD_alt1_rmsd.txt', 'r').readlines():
    a3_AD_open.append(float(i)*10)
for i in open('a3_1sug_AD_alt2_rmsd.txt', 'r').readlines():
    a3_AD_open.append(float(i)*10)

for i in open('a6_1sug_AD_dis11_rmsd.txt', 'r').readlines():
    a6_AD_open.append(float(i)*10)
for i in open('a6_1sug_AD_dis11_2_rmsd.txt', 'r').readlines():
    a6_AD_close.append(float(i)*10)
for i in open('a6_1sug_AD_alt1_rmsd.txt', 'r').readlines():
    a6_AD_open.append(float(i)*10)
for i in open('a6_1sug_AD_alt2_rmsd.txt', 'r').readlines():
    a6_AD_open.append(float(i)*10)

for i in open('a7_1sug_AD_dis11_rmsd.txt', 'r').readlines():
    a7_AD_open.append(float(i)*10)
for i in open('a7_1sug_AD_dis11_2_rmsd.txt', 'r').readlines():
    a7_AD_close.append(float(i)*10)
for i in open('a7_1sug_AD_alt1_rmsd.txt', 'r').readlines():
    a7_AD_open.append(float(i)*10)
for i in open('a7_1sug_AD_alt2_rmsd.txt', 'r').readlines():
    a7_AD_open.append(float(i)*10)

#Input BBR bound files
beta_BBR_open, L11_BBR_open, a3_BBR_open, a6_BBR_open, a7_BBR_open = [],[],[],[],[]
beta_BBR_close, L11_BBR_close, a3_BBR_close, a6_BBR_close, a7_BBR_close = [],[],[],[],[]

for i in open('beta_BBR_rmsd.txt', 'r').readlines():
    beta_BBR_open.append(float(i)*10)
for i in open('beta_1sug_BBR_rmsd.txt', 'r').readlines():
    beta_BBR_close.append(float(i)*10)
for i in open('beta_BBR_dis9_rmsd.txt', 'r').readlines():
    beta_BBR_open.append(float(i)*10)
for i in open('beta_1sug_BBR_dis7_rmsd.txt', 'r').readlines():
    beta_BBR_close.append(float(i)*10)

for i in open('L11_BBR_rmsd.txt', 'r').readlines():
    L11_BBR_open.append(float(i)*10)
for i in open('L11_1sug_BBR_rmsd.txt', 'r').readlines():
    L11_BBR_close.append(float(i)*10)
for i in open('L11_BBR_dis9_rmsd.txt', 'r').readlines():
    L11_BBR_open.append(float(i)*10)
for i in open('L11_1sug_BBR_dis7_rmsd.txt', 'r').readlines():
    L11_BBR_close.append(float(i)*10)

for i in open('a3_BBR_rmsd.txt', 'r').readlines():
    a3_BBR_open.append(float(i)*10)
for i in open('a3_1sug_BBR_rmsd.txt', 'r').readlines():
    a3_BBR_close.append(float(i)*10)
for i in open('a3_BBR_dis9_rmsd.txt', 'r').readlines():
    a3_BBR_open.append(float(i)*10)
for i in open('a3_1sug_BBR_dis7_rmsd.txt', 'r').readlines():
    a3_BBR_close.append(float(i)*10)

for i in open('a6_BBR_rmsd.txt', 'r').readlines():
    a6_BBR_open.append(float(i)*10)
for i in open('a6_1sug_BBR_rmsd.txt', 'r').readlines():
    a6_BBR_close.append(float(i)*10)
for i in open('a6_BBR_dis9_rmsd.txt', 'r').readlines():
    a6_BBR_open.append(float(i)*10)
for i in open('a6_1sug_BBR_dis7_rmsd.txt', 'r').readlines():
    a6_BBR_close.append(float(i)*10)

for i in open('a7_BBR_rmsd.txt', 'r').readlines():
    a7_BBR_open.append(float(i)*10)
for i in open('a7_1sug_BBR_rmsd.txt', 'r').readlines():
    a7_BBR_close.append(float(i)*10)
for i in open('a7_BBR_dis9_rmsd.txt', 'r').readlines():
    a7_BBR_open.append(float(i)*10)
for i in open('a7_1sug_BBR_dis7_rmsd.txt', 'r').readlines():
    a7_BBR_close.append(float(i)*10)

#Calculate Mean and SEN for each group
beta = np.array([np.mean(beta_Apo_open), np.mean(beta_Apo_close), np.mean(beta_AD_open), np.mean(beta_AD_close), np.mean(beta_BBR_open), 
    np.mean(beta_BBR_close)])
L11 = np.array([np.mean(L11_Apo_open), np.mean(L11_Apo_close), np.mean(L11_AD_open), np.mean(L11_AD_close), np.mean(L11_BBR_open), 
    np.mean(L11_BBR_close)])
a3 = np.array([np.mean(a3_Apo_open), np.mean(a3_Apo_close), np.mean(a3_AD_open), np.mean(a3_AD_close), np.mean(a3_BBR_open), np.mean(a3_BBR_close)])
a6 = np.array([np.mean(a6_Apo_open), np.mean(a6_Apo_close), np.mean(a6_AD_open), np.mean(a6_AD_close), np.mean(a6_BBR_open), np.mean(a6_BBR_close)])
a7 = np.array([np.mean(a7_Apo_open), np.mean(a7_Apo_close), np.mean(a7_AD_open), np.mean(a7_AD_close), np.mean(a7_BBR_open), np.mean(a7_BBR_close)])

beta_sem = np.array([stats.sem(beta_Apo_open), stats.sem(beta_Apo_close), stats.sem(beta_AD_open), stats.sem(beta_AD_close), stats.sem(beta_BBR_open), 
    stats.sem(beta_BBR_close)])
L11_sem = np.array([stats.sem(L11_Apo_open), stats.sem(L11_Apo_close), stats.sem(L11_AD_open), stats.sem(L11_AD_close), stats.sem(L11_BBR_open), 
    stats.sem(L11_BBR_close)])
a3_sem = np.array([stats.sem(a3_Apo_open), stats.sem(a3_Apo_close), stats.sem(a3_AD_open), stats.sem(a3_AD_close), stats.sem(a3_BBR_open), 
    stats.sem(a3_BBR_close)])
a6_sem = np.array([stats.sem(a6_Apo_open), stats.sem(a6_Apo_close), stats.sem(a6_AD_open), stats.sem(a6_AD_close), stats.sem(a6_BBR_open), 
    stats.sem(a6_BBR_close)])
a7_sem = np.array([stats.sem(a7_Apo_open), stats.sem(a7_Apo_close), stats.sem(a7_AD_open), stats.sem(a7_AD_close), stats.sem(a7_BBR_open), 
    stats.sem(a7_BBR_close)])

#Run t-test on Subset of groups for beta sheet region
st_beta1, p_beta1 = stats.ttest_ind(beta_Apo_open, beta_AD_open, equal_var = False) #Apo Open vs AD Open
st_beta2, p_beta2 = stats.ttest_ind(beta_Apo_close, beta_AD_open, equal_var = False) #Apo Closed vs AD Open
st_beta3, p_beta3 = stats.ttest_ind(beta_AD_close, beta_AD_open, equal_var = False) #AD Open vs AD Closed
st_beta4, p_beta4 = stats.ttest_ind(beta_Apo_open, beta_BBR_open, equal_var = False) #Apo Open vs BBR Open
st_beta5, p_beta5 = stats.ttest_ind(beta_Apo_close, beta_BBR_open, equal_var = False) #Apo Closed vs BBR Open
st_beta6, p_beta6 = stats.ttest_ind(beta_BBR_open, beta_BBR_close, equal_var = False) #BBR Open vs BBR Closed
st_beta7, p_beta7 = stats.ttest_ind(beta_BBR_open, beta_AD_open, equal_var = False) #AD Open vs BBR Open

#Print p-values to file
output = open('p_values.txt', 'w')
output.write('Beta Sheet RMSD\n')
output.write('Apo Open vs AD Open: ' + str(p_beta1))
output.write('Apo Close vs AD Open: ' + str(p_beta2))
output.write('Apo Close vs AD Open: ' + str(p_beta3))
output.write('Apo Open vs BBR Open: ' + str(p_beta4))
output.write('Apo Close vs BBR Open: ' + str(p_beta5))
output.write('BBR Open vs BBR_close: ' + str(p_beta6))
output.write('BBR Open vs AD Open: ' + str(p_beta7))

#Run t-test on Subset of groups for beta sheet region
st_L111, p_L111 = stats.ttest_ind(L11_Apo_open, L11_AD_open, equal_var = False) #Apo Open vs AD Open
st_L112, p_L112 = stats.ttest_ind(L11_Apo_close, L11_AD_open, equal_var = False) #Apo Closed vs AD Open
st_L113, p_L113 = stats.ttest_ind(L11_AD_close, L11_AD_open, equal_var = False) #AD Open vs AD Closed
st_L114, p_L114 = stats.ttest_ind(L11_Apo_open, L11_BBR_open, equal_var = False) #Apo Open vs BBR Open
st_L115, p_L115 = stats.ttest_ind(L11_Apo_close, L11_BBR_open, equal_var = False) #Apo Closed vs BBR Open
st_L116, p_L116 = stats.ttest_ind(L11_BBR_open, L11_BBR_close, equal_var = False) #BBR Open vs BBR Closed
st_L117, p_L117 = stats.ttest_ind(L11_BBR_open, L11_AD_open, equal_var = False) #AD Open vs BBR Open

#Print p-values to file
output.write('L11 loop RMSD\n')
output.write('Apo Open vs AD Open: ' + str(p_L111))
output.write('Apo Close vs AD Open: ' + str(p_L112))
output.write('Apo Close vs AD Open: ' + str(p_L113))
output.write('Apo Open vs BBR Open: ' + str(p_L114))
output.write('Apo Close vs BBR Open: ' + str(p_L115))
output.write('BBR Open vs BBR_close: ' + str(p_L116))
output.write('BBR Open vs AD Open: ' + str(p_L117))

#Run t-test on Subset of groups for alpha 3 helix
st_a31, p_a31 = stats.ttest_ind(a3_Apo_open, a3_AD_open, equal_var = False) #Apo Open vs AD Open
st_a32, p_a32 = stats.ttest_ind(a3_Apo_close, a3_AD_open, equal_var = False) #Apo Closed vs AD Open
st_a33, p_a33 = stats.ttest_ind(a3_AD_close, a3_AD_open, equal_var = False) #AD Open vs AD Closed
st_a34, p_a34 = stats.ttest_ind(a3_Apo_open, a3_BBR_open, equal_var = False) #Apo Open vs BBR Open
st_a35, p_a35 = stats.ttest_ind(a3_Apo_close, a3_BBR_open, equal_var = False) #Apo Closed vs BBR Open
st_a36, p_a36 = stats.ttest_ind(a3_BBR_open, a3_BBR_close, equal_var = False) #BBR Open vs BBR Closed
st_a37, p_a37 = stats.ttest_ind(a3_BBR_open, a3_AD_open, equal_var = False) #AD Open vs BBR Open

#Print p-values to file
output = open('p_values.txt', 'w')
output.write('a3 Helix RMSD\n')
output.write('Apo Open vs AD Open: ' + str(p_a31))
output.write('Apo Close vs AD Open: ' + str(p_a32))
output.write('Apo Close vs AD Open: ' + str(p_a33))
output.write('Apo Open vs BBR Open: ' + str(p_a34))
output.write('Apo Close vs BBR Open: ' + str(p_a35))
output.write('BBR Open vs BBR_close: ' + str(p_a36))
output.write('BBR Open vs AD Open: ' + str(p_a37))

#Run t-test on Subset of groups for alpha 3 helix
st_a61, p_a61 = stats.ttest_ind(a6_Apo_open, a6_AD_open, equal_var = False) #Apo Open vs AD Open
st_a62, p_a62 = stats.ttest_ind(a6_Apo_close, a6_AD_open, equal_var = False) #Apo Closed vs AD Open
st_a63, p_a63 = stats.ttest_ind(a6_AD_close, a6_AD_open, equal_var = False) #AD Open vs AD Closed
st_a64, p_a64 = stats.ttest_ind(a6_Apo_open, a6_BBR_open, equal_var = False) #Apo Open vs BBR Open
st_a65, p_a65 = stats.ttest_ind(a6_Apo_close, a6_BBR_open, equal_var = False) #Apo Closed vs BBR Open
st_a66, p_a66 = stats.ttest_ind(a6_BBR_open, a6_BBR_close, equal_var = False) #BBR Open vs BBR Closed
st_a67, p_a67 = stats.ttest_ind(a6_BBR_open, a6_AD_open, equal_var = False) #AD Open vs BBR Open

#Print p-values to file
output = open('p_values.txt', 'w')
output.write('a3 Helix RMSD\n')
output.write('Apo Open vs AD Open: ' + str(p_a61))
output.write('Apo Close vs AD Open: ' + str(p_a62))
output.write('Apo Close vs AD Open: ' + str(p_a63))
output.write('Apo Open vs BBR Open: ' + str(p_a64))
output.write('Apo Close vs BBR Open: ' + str(p_a65))
output.write('BBR Open vs BBR_close: ' + str(p_a66))
output.write('BBR Open vs AD Open: ' + str(p_a67))

#Run t-test on Subset of groups for alpha 3 helix
st_a71, p_a71 = stats.ttest_ind(a7_Apo_open, a7_AD_open, equal_var = False) #Apo Open vs AD Open
st_a72, p_a72 = stats.ttest_ind(a7_Apo_close, a7_AD_open, equal_var = False) #Apo Closed vs AD Open
st_a73, p_a73 = stats.ttest_ind(a7_AD_close, a7_AD_open, equal_var = False) #AD Open vs AD Closed
st_a74, p_a74 = stats.ttest_ind(a7_Apo_open, a7_BBR_open, equal_var = False) #Apo Open vs BBR Open
st_a75, p_a75 = stats.ttest_ind(a7_Apo_close, a7_BBR_open, equal_var = False) #Apo Closed vs BBR Open
st_a76, p_a76 = stats.ttest_ind(a7_BBR_open, a7_BBR_close, equal_var = False) #BBR Open vs BBR Closed
st_a77, p_a77 = stats.ttest_ind(a7_BBR_open, a7_AD_open, equal_var = False) #AD Open vs BBR Open

#Print p-values to file
output = open('p_values.txt', 'w')
output.write('a3 Helix RMSD\n')
output.write('Apo Open vs AD Open: ' + str(p_a71))
output.write('Apo Close vs AD Open: ' + str(p_a72))
output.write('Apo Close vs AD Open: ' + str(p_a73))
output.write('Apo Open vs BBR Open: ' + str(p_a74))
output.write('Apo Close vs BBR Open: ' + str(p_a75))
output.write('BBR Open vs BBR_close: ' + str(p_a76))
output.write('BBR Open vs AD Open: ' + str(p_a77))

#Plot Apo and AD
num = [5, 10, 15, 20]
Method = ['Apo Open', 'Apo Closed', 'AD Open', 'AD Closed']
fig = plt.figure()
plt.title('RMSD of Beta Sheet Relative to Individual Apo Closed')
plt.ylabel('RMSD (A)')
plt.bar(num, beta[[0, 1, 2, 3]], color = ['blue', 'red', 'blue', 'red'], width=4.5)
plt.errorbar(num, beta[[0, 1, 2, 3]], yerr=beta_sem[[0, 1, 2, 3]], fmt='o', color='black')
plt.xticks(num, Method, fontsize=8)
if p_beta1 <= 0.05 and p_beta1 > 0.01:
    x1, x2 = 5, 15 #Columns for Open Apo + a7 and Closed Apo + a7
    y, h, col = (1.1*beta[[0, 2]].max()), 0.2, 'b'
    plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
    plt.text((x1+x2)*0.5, y+h, "*" , ha='center', va='bottom', color=col)
if p_beta1 <= 0.01 and p_beta1 > 0.001:
    x1, x2 = 5, 15 #Columns for Open Apo + a7 and Closed Apo + a7
    y, h, col = (1.1*beta[[0, 2]].max()), 0.2, 'b'
    plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
    plt.text((x1+x2)*0.5, y+h, "**" , ha='center', va='bottom', color=col)
if p_beta1 <= 0.001:
    x1, x2 = 5, 15 #Columns for Open Apo + a7 and Closed Apo + a7
    y, h, col = (1.1*beta[[0, 2]].max()), 0.2, 'b'
    plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
    plt.text((x1+x2)*0.5, y+h, "***" , ha='center', va='bottom', color=col)
if p_beta2 <= 0.05 and p_beta2 > 0.01:
    x1, x2 = 10, 15 #Columns for Open Apo + a7 and Closed Apo + a7
    y, h, col = (1.1*beta[[1, 2]].max()), 0.1, 'r'
    plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
    plt.text((x1+x2)*0.5, y+h, "*" , ha='center', va='bottom', color=col)
if p_beta2 <= 0.02 and p_beta2 > 0.001:
    x1, x2 = 10, 15 #Columns for Open Apo + a7 and Closed Apo + a7
    y, h, col = (1.1*beta[[1, 2]].max()), 0.1, 'r'
    plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
    plt.text((x1+x2)*0.5, y+h, "**" , ha='center', va='bottom', color=col)
if p_beta2 <= 0.001:
    x1, x2 = 10, 15 #Columns for Open Apo + a7 and Closed Apo + a7
    y, h, col = (1.1*beta[[1, 2]].max()), 0.1, 'r'
    plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
    plt.text((x1+x2)*0.5, y+h, "***" , ha='center', va='bottom', color=col)
fig.savefig('Beta_rmsd.png')
plt.close(fig)

fig = plt.figure()
plt.title('RMSD of L11 Loop Relative to Individual Apo Closed')
plt.ylabel('RMSD (A)')
plt.bar(num, L11[[0, 1, 2, 3]], color = ['blue', 'red', 'blue', 'red'], width=4.5)
plt.errorbar(num, L11[[0, 1, 2, 3]], yerr=L11_sem[[0, 1, 2, 3]], fmt='o', color='black')
plt.xticks(num, Method, fontsize=8)
if p_L111 <= 0.05 and p_L111 > 0.01:
    x1, x2 = 5, 15 #Columns for Open Apo + a7 and Closed Apo + a7
    y, h, col = (1.1*L11[[0, 2]].max()), 0.2, 'b'
    plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
    plt.text((x1+x2)*0.5, y+h, "*" , ha='center', va='bottom', color=col)
if p_L111 <= 0.01 and p_L111 > 0.001:
    x1, x2 = 5, 15 #Columns for Open Apo + a7 and Closed Apo + a7
    y, h, col = (1.1*L11[[0, 2]].max()), 0.2, 'b'
    plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
    plt.text((x1+x2)*0.5, y+h, "**" , ha='center', va='bottom', color=col)
if p_L111 <= 0.001:
    x1, x2 = 5, 15 #Columns for Open Apo + a7 and Closed Apo + a7
    y, h, col = (1.1*L11[[0, 2]].max()), 0.2, 'b'
    plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
    plt.text((x1+x2)*0.5, y+h, "***" , ha='center', va='bottom', color=col)
if p_L112 <= 0.05 and p_L112 > 0.01:
    x1, x2 = 10, 15 #Columns for Open Apo + a7 and Closed Apo + a7
    y, h, col = (1.1*L11[[1, 2]].max()), 0.1, 'r'
    plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
    plt.text((x1+x2)*0.5, y+h, "*" , ha='center', va='bottom', color=col)
if p_L112 <= 0.02 and p_L112 > 0.001:
    x1, x2 = 10, 15 #Columns for Open Apo + a7 and Closed Apo + a7
    y, h, col = (1.1*L11[[1, 2]].max()), 0.1, 'r'
    plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
    plt.text((x1+x2)*0.5, y+h, "**" , ha='center', va='bottom', color=col)
if p_L112 <= 0.001:
    x1, x2 = 10, 15 #Columns for Open Apo + a7 and Closed Apo + a7
    y, h, col = (1.1*L11[[1, 2]].max()), 0.1, 'r'
    plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
    plt.text((x1+x2)*0.5, y+h, "***" , ha='center', va='bottom', color=col)
fig.savefig('L11_rmsd.png')
plt.close(fig)

fig = plt.figure()
plt.title('RMSD of a3 Helix Relative to Individual Apo Closed')
plt.ylabel('RMSD (A)')
plt.bar(num, a3[[0, 1, 2, 3]], color = ['blue', 'red', 'blue', 'red'], width=4.5)
plt.errorbar(num, a3[[0, 1, 2, 3]], yerr=a3_sem[[0, 1, 2, 3]], fmt='o', color='black')
plt.xticks(num, Method, fontsize=8)
if p_a31 <= 0.05 and p_a31 > 0.01:
    x1, x2 = 5, 15 #Columns for Open Apo + a7 and Closed Apo + a7
    y, h, col = (1.1*a3[[0, 2]].max()), 0.2, 'b'
    plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
    plt.text((x1+x2)*0.5, y+h, "*" , ha='center', va='bottom', color=col)
if p_a31 <= 0.01 and p_a31 > 0.001:
    x1, x2 = 5, 15 #Columns for Open Apo + a7 and Closed Apo + a7
    y, h, col = (1.1*a3[[0, 2]].max()), 0.2, 'b'
    plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
    plt.text((x1+x2)*0.5, y+h, "**" , ha='center', va='bottom', color=col)
if p_a31 <= 0.001:
    x1, x2 = 5, 15 #Columns for Open Apo + a7 and Closed Apo + a7
    y, h, col = (1.1*a3[[0, 2]].max()), 0.2, 'b'
    plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
    plt.text((x1+x2)*0.5, y+h, "***" , ha='center', va='bottom', color=col)
if p_a32 <= 0.05 and p_a32 > 0.01:
    x1, x2 = 10, 15 #Columns for Open Apo + a7 and Closed Apo + a7
    y, h, col = (1.1*a3[[1, 2]].max()), 0.1, 'r'
    plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
    plt.text((x1+x2)*0.5, y+h, "*" , ha='center', va='bottom', color=col)
if p_a32 <= 0.02 and p_a32 > 0.001:
    x1, x2 = 10, 15 #Columns for Open Apo + a7 and Closed Apo + a7
    y, h, col = (1.1*a3[[1, 2]].max()), 0.1, 'r'
    plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
    plt.text((x1+x2)*0.5, y+h, "**" , ha='center', va='bottom', color=col)
if p_a32 <= 0.001:
    x1, x2 = 10, 15 #Columns for Open Apo + a7 and Closed Apo + a7
    y, h, col = (1.1*a3[[1, 2]].max()), 0.1, 'r'
    plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
    plt.text((x1+x2)*0.5, y+h, "***" , ha='center', va='bottom', color=col)
fig.savefig('a3_rmsd.png')
plt.close(fig)

fig = plt.figure()
plt.title('RMSD of a6 Helix Relative to Individual Apo Closed')
plt.ylabel('RMSD (A)')
plt.bar(num, a6[[0, 1, 2, 3]], color = ['blue', 'red', 'blue', 'red'], width=4.5)
plt.errorbar(num, a6[[0, 1, 2, 3]], yerr=a6_sem[[0, 1, 2, 3]], fmt='o', color='black')
plt.xticks(num, Method, fontsize=8)
if p_a61 <= 0.05 and p_a61 > 0.01:
    x1, x2 = 5, 15 #Columns for Open Apo + a7 and Closed Apo + a7
    y, h, col = (1.1*a6[[0, 2]].max()), 0.2, 'b'
    plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
    plt.text((x1+x2)*0.5, y+h, "*" , ha='center', va='bottom', color=col)
if p_a61 <= 0.01 and p_a61 > 0.001:
    x1, x2 = 5, 15 #Columns for Open Apo + a7 and Closed Apo + a7
    y, h, col = (1.1*a6[[0, 2]].max()), 0.2, 'b'
    plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
    plt.text((x1+x2)*0.5, y+h, "**" , ha='center', va='bottom', color=col)
if p_a61 <= 0.001:
    x1, x2 = 5, 15 #Columns for Open Apo + a7 and Closed Apo + a7
    y, h, col = (1.1*a6[[0, 2]].max()), 0.2, 'b'
    plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
    plt.text((x1+x2)*0.5, y+h, "***" , ha='center', va='bottom', color=col)
if p_a62 <= 0.05 and p_a62 > 0.01:
    x1, x2 = 10, 15 #Columns for Open Apo + a7 and Closed Apo + a7
    y, h, col = (1.1*a6[[1, 2]].max()), 0.1, 'r'
    plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
    plt.text((x1+x2)*0.5, y+h, "*" , ha='center', va='bottom', color=col)
if p_a62 <= 0.02 and p_a62 > 0.001:
    x1, x2 = 10, 15 #Columns for Open Apo + a7 and Closed Apo + a7
    y, h, col = (1.1*a6[[1, 2]].max()), 0.1, 'r'
    plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
    plt.text((x1+x2)*0.5, y+h, "**" , ha='center', va='bottom', color=col)
if p_a62 <= 0.001:
    x1, x2 = 10, 15 #Columns for Open Apo + a7 and Closed Apo + a7
    y, h, col = (1.1*a6[[1, 2]].max()), 0.1, 'r'
    plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
    plt.text((x1+x2)*0.5, y+h, "***" , ha='center', va='bottom', color=col)
fig.savefig('a6_rmsd.png')
plt.close(fig)

fig = plt.figure()
plt.title('RMSD of a7 Helix Relative to Individual Apo Closed')
plt.ylabel('RMSD (A)')
plt.bar(num, a7[[0, 1, 2, 3]], color = ['blue', 'red', 'blue', 'red'], width=4.5)
plt.errorbar(num, a7[[0, 1, 2, 3]], yerr=a7_sem[[0, 1, 2, 3]], fmt='o', color='black')
plt.xticks(num, Method, fontsize=8)
if p_a71 <= 0.05 and p_a71 > 0.01:
    x1, x2 = 5, 15 #Columns for Open Apo + a7 and Closed Apo + a7
    y, h, col = (1.1*a7[[0, 2]].max()), 0.2, 'b'
    plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
    plt.text((x1+x2)*0.5, y+h, "*" , ha='center', va='bottom', color=col)
if p_a71 <= 0.01 and p_a71 > 0.001:
    x1, x2 = 5, 15 #Columns for Open Apo + a7 and Closed Apo + a7
    y, h, col = (1.1*a7[[0, 2]].max()), 0.2, 'b'
    plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
    plt.text((x1+x2)*0.5, y+h, "**" , ha='center', va='bottom', color=col)
if p_a71 <= 0.001:
    x1, x2 = 5, 15 #Columns for Open Apo + a7 and Closed Apo + a7
    y, h, col = (1.1*a7[[0, 2]].max()), 0.2, 'b'
    plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
    plt.text((x1+x2)*0.5, y+h, "***" , ha='center', va='bottom', color=col)
if p_a72 <= 0.05 and p_a72 > 0.01:
    x1, x2 = 10, 15 #Columns for Open Apo + a7 and Closed Apo + a7
    y, h, col = (1.1*a7[[1, 2]].max()), 0.1, 'r'
    plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
    plt.text((x1+x2)*0.5, y+h, "*" , ha='center', va='bottom', color=col)
if p_a72 <= 0.02 and p_a72 > 0.001:
    x1, x2 = 10, 15 #Columns for Open Apo + a7 and Closed Apo + a7
    y, h, col = (1.1*a7[[1, 2]].max()), 0.1, 'r'
    plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
    plt.text((x1+x2)*0.5, y+h, "**" , ha='center', va='bottom', color=col)
if p_a72 <= 0.001:
    x1, x2 = 10, 15 #Columns for Open Apo + a7 and Closed Apo + a7
    y, h, col = (1.1*a7[[1, 2]].max()), 0.1, 'r'
    plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
    plt.text((x1+x2)*0.5, y+h, "***" , ha='center', va='bottom', color=col)
fig.savefig('a7_rmsd.png')
plt.close(fig)

#Plot AD and BBR
num = [5, 10, 15]
Method = ['Apo Open', 'Apo Closed', 'BBR Open']
fig = plt.figure()
plt.title('RMSD of Beta Sheet Relative to Individual Apo Closed')
plt.ylabel('RMSD (A)')
plt.bar(num, beta[[0, 1, 4]], color = ['blue', 'red', 'blue'], width=4.5)
plt.errorbar(num, beta[[0, 1, 4]], yerr=beta_sem[[0, 1, 4]], fmt='o', color='black')
plt.xticks(num, Method, fontsize=8)
if p_beta4 <= 0.05 and p_beta4 > 0.01:
    x1, x2 = 5, 15 #Columns for Open Apo + a7 and Closed Apo + a7
    y, h, col = (1.1*beta[[0, 4]].max()), 0.2, 'b'
    plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
    plt.text((x1+x2)*0.5, y+h, "*" , ha='center', va='bottom', color=col)
if p_beta4 <= 0.01 and p_beta4 > 0.001:
    x1, x2 = 5, 15 #Columns for Open Apo + a7 and Closed Apo + a7
    y, h, col = (1.1*beta[[0, 4]].max()), 0.2, 'b'
    plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
    plt.text((x1+x2)*0.5, y+h, "**" , ha='center', va='bottom', color=col)
if p_beta4 <= 0.001:
    x1, x2 = 5, 15 #Columns for Open Apo + a7 and Closed Apo + a7
    y, h, col = (1.1*beta[[0, 4]].max()), 0.2, 'b'
    plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
    plt.text((x1+x2)*0.5, y+h, "***" , ha='center', va='bottom', color=col)
if p_beta5 <= 0.05 and p_beta5 > 0.01:
    x1, x2 = 10, 15 #Columns for Open Apo + a7 and Closed Apo + a7
    y, h, col = (1.1*beta[[1, 4]].max()), 0.1, 'r'
    plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
    plt.text((x1+x2)*0.5, y+h, "*" , ha='center', va='bottom', color=col)
if p_beta5 <= 0.02 and p_beta5 > 0.001:
    x1, x2 = 10, 15 #Columns for Open Apo + a7 and Closed Apo + a7
    y, h, col = (1.1*beta[[1, 4]].max()), 0.1, 'r'
    plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
    plt.text((x1+x2)*0.5, y+h, "**" , ha='center', va='bottom', color=col)
if p_beta5 <= 0.001:
    x1, x2 = 10, 15 #Columns for Open Apo + a7 and Closed Apo + a7
    y, h, col = (1.1*beta[[1, 4]].max()), 0.1, 'r'
    plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
    plt.text((x1+x2)*0.5, y+h, "***" , ha='center', va='bottom', color=col)
fig.savefig('Beta_BBR_rmsd.png')
plt.close(fig)

fig = plt.figure()
plt.title('RMSD of L11 Relative to Individual Apo Closed')
plt.ylabel('RMSD (A)')
plt.bar(num, L11[[0, 1, 4]], color = ['blue', 'red', 'blue'], width=4.5)
plt.errorbar(num, L11[[0, 1, 4]], yerr=L11_sem[[0, 1, 4]], fmt='o', color='black')
plt.xticks(num, Method, fontsize=8)
if p_L114 <= 0.05 and p_L114 > 0.01:
    x1, x2 = 5, 15 #Columns for Open Apo + a7 and Closed Apo + a7
    y, h, col = (1.1*L11[[0, 4]].max()), 0.2, 'b'
    plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
    plt.text((x1+x2)*0.5, y+h, "*" , ha='center', va='bottom', color=col)
if p_L114 <= 0.01 and p_L114 > 0.001:
    x1, x2 = 5, 15 #Columns for Open Apo + a7 and Closed Apo + a7
    y, h, col = (1.1*L11[[0, 4]].max()), 0.2, 'b'
    plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
    plt.text((x1+x2)*0.5, y+h, "**" , ha='center', va='bottom', color=col)
if p_L114 <= 0.001:
    x1, x2 = 5, 15 #Columns for Open Apo + a7 and Closed Apo + a7
    y, h, col = (1.1*L11[[0, 4]].max()), 0.2, 'b'
    plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
    plt.text((x1+x2)*0.5, y+h, "***" , ha='center', va='bottom', color=col)
if p_L115 <= 0.05 and p_L115 > 0.01:
    x1, x2 = 10, 15 #Columns for Open Apo + a7 and Closed Apo + a7
    y, h, col = (1.1*L11[[1, 4]].max()), 0.1, 'r'
    plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
    plt.text((x1+x2)*0.5, y+h, "*" , ha='center', va='bottom', color=col)
if p_L115 <= 0.02 and p_L115 > 0.001:
    x1, x2 = 10, 15 #Columns for Open Apo + a7 and Closed Apo + a7
    y, h, col = (1.1*L11[[1, 4]].max()), 0.1, 'r'
    plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
    plt.text((x1+x2)*0.5, y+h, "**" , ha='center', va='bottom', color=col)
if p_L115 <= 0.001:
    x1, x2 = 10, 15 #Columns for Open Apo + a7 and Closed Apo + a7
    y, h, col = (1.1*L11[[1, 4]].max()), 0.1, 'r'
    plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
    plt.text((x1+x2)*0.5, y+h, "***" , ha='center', va='bottom', color=col)
fig.savefig('L11_BBR_rmsd.png')
plt.close(fig)

fig = plt.figure()
plt.title('RMSD of a3 Relative to Individual Apo Closed')
plt.ylabel('RMSD (A)')
plt.bar(num, a3[[0, 1, 4]], color = ['blue', 'red', 'blue'], width=4.5)
plt.errorbar(num, a3[[0, 1, 4]], yerr=a3_sem[[0, 1, 4]], fmt='o', color='black')
plt.xticks(num, Method, fontsize=8)
if p_a34 <= 0.05 and p_a34 > 0.01:
    x1, x2 = 5, 15 #Columns for Open Apo + a7 and Closed Apo + a7
    y, h, col = (1.1*a3[[0, 4]].max()), 0.2, 'b'
    plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
    plt.text((x1+x2)*0.5, y+h, "*" , ha='center', va='bottom', color=col)
if p_a34 <= 0.01 and p_a34 > 0.001:
    x1, x2 = 5, 15 #Columns for Open Apo + a7 and Closed Apo + a7
    y, h, col = (1.1*a3[[0, 4]].max()), 0.2, 'b'
    plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
    plt.text((x1+x2)*0.5, y+h, "**" , ha='center', va='bottom', color=col)
if p_a34 <= 0.001:
    x1, x2 = 5, 15 #Columns for Open Apo + a7 and Closed Apo + a7
    y, h, col = (1.1*a3[[0, 4]].max()), 0.2, 'b'
    plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
    plt.text((x1+x2)*0.5, y+h, "***" , ha='center', va='bottom', color=col)
if p_a35 <= 0.05 and p_a35 > 0.01:
    x1, x2 = 10, 15 #Columns for Open Apo + a7 and Closed Apo + a7
    y, h, col = (1.1*a3[[1, 4]].max()), 0.1, 'r'
    plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
    plt.text((x1+x2)*0.5, y+h, "*" , ha='center', va='bottom', color=col)
if p_a35 <= 0.02 and p_a35 > 0.001:
    x1, x2 = 10, 15 #Columns for Open Apo + a7 and Closed Apo + a7
    y, h, col = (1.1*a3[[1, 4]].max()), 0.1, 'r'
    plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
    plt.text((x1+x2)*0.5, y+h, "**" , ha='center', va='bottom', color=col)
if p_a35 <= 0.001:
    x1, x2 = 10, 15 #Columns for Open Apo + a7 and Closed Apo + a7
    y, h, col = (1.1*a3[[1, 4]].max()), 0.1, 'r'
    plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
    plt.text((x1+x2)*0.5, y+h, "***" , ha='center', va='bottom', color=col)
fig.savefig('a3_BBR_rmsd.png')
plt.close(fig)

fig = plt.figure()
plt.title('RMSD of a6 Relative to Individual Apo Closed')
plt.ylabel('RMSD (A)')
plt.bar(num, a6[[0, 1, 4]], color = ['blue', 'red', 'blue'], width=4.5)
plt.errorbar(num, a6[[0, 1, 4]], yerr=a6_sem[[0, 1, 4]], fmt='o', color='black')
plt.xticks(num, Method, fontsize=8)
if p_a64 <= 0.05 and p_a64 > 0.01:
    x1, x2 = 5, 15 #Columns for Open Apo + a7 and Closed Apo + a7
    y, h, col = (1.1*a6[[0, 4]].max()), 0.2, 'b'
    plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
    plt.text((x1+x2)*0.5, y+h, "*" , ha='center', va='bottom', color=col)
if p_a64 <= 0.01 and p_a64 > 0.001:
    x1, x2 = 5, 15 #Columns for Open Apo + a7 and Closed Apo + a7
    y, h, col = (1.1*a6[[0, 4]].max()), 0.2, 'b'
    plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
    plt.text((x1+x2)*0.5, y+h, "**" , ha='center', va='bottom', color=col)
if p_a64 <= 0.001:
    x1, x2 = 5, 15 #Columns for Open Apo + a7 and Closed Apo + a7
    y, h, col = (1.1*a6[[0, 4]].max()), 0.2, 'b'
    plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
    plt.text((x1+x2)*0.5, y+h, "***" , ha='center', va='bottom', color=col)
if p_a65 <= 0.05 and p_a65 > 0.01:
    x1, x2 = 10, 15 #Columns for Open Apo + a7 and Closed Apo + a7
    y, h, col = (1.1*a6[[1, 4]].max()), 0.1, 'r'
    plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
    plt.text((x1+x2)*0.5, y+h, "*" , ha='center', va='bottom', color=col)
if p_a65 <= 0.02 and p_a65 > 0.001:
    x1, x2 = 10, 15 #Columns for Open Apo + a7 and Closed Apo + a7
    y, h, col = (1.1*a6[[1, 4]].max()), 0.1, 'r'
    plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
    plt.text((x1+x2)*0.5, y+h, "**" , ha='center', va='bottom', color=col)
if p_a65 <= 0.001:
    x1, x2 = 10, 15 #Columns for Open Apo + a7 and Closed Apo + a7
    y, h, col = (1.1*a6[[1, 4]].max()), 0.1, 'r'
    plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
    plt.text((x1+x2)*0.5, y+h, "***" , ha='center', va='bottom', color=col)
fig.savefig('a6_BBR_rmsd.png')
plt.close(fig)

fig = plt.figure()
plt.title('RMSD of a7 Relative to Individual Apo Closed')
plt.ylabel('RMSD (A)')
plt.bar(num, a7[[0, 1, 4]], color = ['blue', 'red', 'blue'], width=4.5)
plt.errorbar(num, a7[[0, 1, 4]], yerr=a7_sem[[0, 1, 4]], fmt='o', color='black')
plt.xticks(num, Method, fontsize=8)
if p_a74 <= 0.05 and p_a74 > 0.01:
    x1, x2 = 5, 15 #Columns for Open Apo + a7 and Closed Apo + a7
    y, h, col = (1.1*a7[[0, 4]].max()), 0.2, 'b'
    plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
    plt.text((x1+x2)*0.5, y+h, "*" , ha='center', va='bottom', color=col)
if p_a74 <= 0.01 and p_a74 > 0.001:
    x1, x2 = 5, 15 #Columns for Open Apo + a7 and Closed Apo + a7
    y, h, col = (1.1*a7[[0, 4]].max()), 0.2, 'b'
    plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
    plt.text((x1+x2)*0.5, y+h, "**" , ha='center', va='bottom', color=col)
if p_a74 <= 0.001:
    x1, x2 = 5, 15 #Columns for Open Apo + a7 and Closed Apo + a7
    y, h, col = (1.1*a7[[0, 4]].max()), 0.2, 'b'
    plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
    plt.text((x1+x2)*0.5, y+h, "***" , ha='center', va='bottom', color=col)
if p_a75 <= 0.05 and p_a75 > 0.01:
    x1, x2 = 10, 15 #Columns for Open Apo + a7 and Closed Apo + a7
    y, h, col = (1.1*a7[[1, 4]].max()), 0.1, 'r'
    plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
    plt.text((x1+x2)*0.5, y+h, "*" , ha='center', va='bottom', color=col)
if p_a75 <= 0.02 and p_a75 > 0.001:
    x1, x2 = 10, 15 #Columns for Open Apo + a7 and Closed Apo + a7
    y, h, col = (1.1*a7[[1, 4]].max()), 0.1, 'r'
    plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
    plt.text((x1+x2)*0.5, y+h, "**" , ha='center', va='bottom', color=col)
if p_a75 <= 0.001:
    x1, x2 = 10, 15 #Columns for Open Apo + a7 and Closed Apo + a7
    y, h, col = (1.1*a7[[1, 4]].max()), 0.1, 'r'
    plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
    plt.text((x1+x2)*0.5, y+h, "***" , ha='center', va='bottom', color=col)
fig.savefig('a7_BBR_rmsd.png')
plt.close(fig)

