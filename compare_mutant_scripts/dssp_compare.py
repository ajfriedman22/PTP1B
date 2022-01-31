#!/ usr / bin / env python

#Import packages
from matplotlib import pyplot as plt
import numpy as np

#Function for determining helix orientation
def determine_helix(list_dssp):
    num_struct = np.zeros(10)
    num_alpha_helix = np.zeros(10)
    per_struct = np.zeros(8)
    per = np.zeros(8)
    for char in list_dssp:
        if char[0]=='H':
            num_alpha_helix[0] +=1
        if char[0] == 'T' or char[0] == 'G' or char[0] == 'I':
            num_struct[0] += 1
        if char[2]=='H':
            num_alpha_helix[1] +=1
        if char[2] == 'T' or char[2] == 'G' or char[2] == 'I':
            num_struct[1] += 1
        if char[4]=='H':
            num_alpha_helix[2] += 1
        if char[4] == 'T' or char[4] == 'G' or char[4] == 'I':
            num_struct[2] += 1
        if char[6]=='H':
            num_alpha_helix[3] +=1
        if char[6] == 'T' or char[6] == 'G' or char[6] == 'I':
            num_struct[3] += 1
        if char[8]=='H':
            num_alpha_helix[4] += 1
        if char[8] == 'T' or char[8] == 'G' or char[8] == 'I':
            num_struct[4] += 1
        if char[10]=='H':
            num_alpha_helix[5] += 1
        if char[10] == 'T' or char[10] == 'G' or char[10] == 'I':
            num_struct[5] += 1
        if char[12]=='H':
            num_alpha_helix[6] += 1
        if char[12] == 'T' or char[12] == 'G' or char[12] == 'I':
            num_struct[6] += 1
        if char[14]=='H':
            num_alpha_helix[7] +=1
        if char[14] == 'T' or char[14] == 'G' or char[14] == 'I':
            num_struct[7] += 1
        if char[16]=='H':
            num_alpha_helix[8] +=1
        if char[16] == 'T' or char[16] == 'G' or char[16] == 'I':
            num_struct[8] += 1
        if char[18]=='H':
            num_alpha_helix[9] +=1
        if char[18] == 'T' or char[18] == 'G' or char[18] == 'I':
            num_struct[9] += 1
    list_len = len(list_dssp)
    
    per[0] = round(100*num_alpha_helix[2]/list_len)
    per[1] = round(100*num_alpha_helix[3]/list_len)
    per[2] = round(100*num_alpha_helix[4]/list_len) 
    per[3] = round(100*num_alpha_helix[5]/list_len) 
    per[4] = round(100*num_alpha_helix[6]/list_len) 
    per[5] = round(100*num_alpha_helix[7]/list_len) 
    per[6] = round(100*num_alpha_helix[8]/list_len) 
    per[7] = round(100*num_alpha_helix[9]/list_len)

    per_struct[0] = round(100*num_struct[2]/list_len)
    per_struct[1] = round(100*num_struct[3]/list_len)
    per_struct[2] = round(100*num_struct[4]/list_len) 
    per_struct[3] = round(100*num_struct[5]/list_len) 
    per_struct[4] = round(100*num_struct[6]/list_len) 
    per_struct[5] = round(100*num_struct[7]/list_len) 
    per_struct[6] = round(100*num_struct[8]/list_len) 
    per_struct[7] = round(100*num_struct[9]/list_len)
    return num_alpha_helix, num_struct, per, per_struct

#Read in data from input files
WT_AD = open('../../WT/AD/DSSP_WT.txt', 'r').readlines()
WT_BBR = open('../../WT/BBR/DSSP_WT_BBR.txt', 'r').readlines()
F196A_AD = open('../../F196A/AD/DSSP_F196A.txt', 'r').readlines()
F196A_BBR = open('../../F196A/BBR/DSSP_F196A_BBR.txt', 'r').readlines()
L192A_AD = open('../../L192A/AD/DSSP_L192A.txt', 'r').readlines()
L192A_BBR = open('../../L192A/BBR/DSSP_L192A_BBR.txt', 'r').readlines()
L192F_AD = open('../../L192F/AD/DSSP_L192F.txt', 'r').readlines()
L192F_BBR = open('../../L192F/BBR/DSSP_L192F_BBR.txt', 'r').readlines()
L192N_AD = open('../../L192N/AD/DSSP_L192N.txt', 'r').readlines()
L192N_BBR = open('../../L192N/BBR/DSSP_L192N_BBR.txt', 'r').readlines()
L195A_AD = open('../../L195A/AD/DSSP_L195A.txt', 'r').readlines()
L195A_BBR = open('../../L195A/BBR/DSSP_L195A_BBR.txt', 'r').readlines()
L195F_AD = open('../../L195F/AD/DSSP_L195F.txt', 'r').readlines()
L195F_BBR = open('../../L195F/BBR/DSSP_L195F_BBR.txt', 'r').readlines()
L195N_AD = open('../../L195N/AD/DSSP_L195N.txt', 'r').readlines()
L195N_BBR = open('../../L195N/BBR/DSSP_L195N_BBR.txt', 'r').readlines()
S286A_AD = open('../../S286A/AD/DSSP_S286A_AD.txt', 'r').readlines()
S286A_BBR = open('../../S286A/BBR/DSSP_S286A_BBR.txt', 'r').readlines()
F280Y_AD = open('../../F280Y/AD/DSSP_F280Y_AD.txt', 'r').readlines()
F280Y_BBR = open('../../F280Y/BBR/DSSP_F280Y_BBR.txt', 'r').readlines()
E276L_AD = open('../../E276L/AD/DSSP_E276L_AD.txt', 'r').readlines()
E276L_BBR = open('../../E276L/BBR/DSSP_E276L_BBR.txt', 'r').readlines()
E276F_AD = open('../../E276F/AD/DSSP_E276F_AD.txt', 'r').readlines()
E276F_BBR = open('../../E276F/BBR/DSSP_E276F_BBR.txt', 'r').readlines()
K279M_AD = open('../../K279M/AD/DSSP_K279M_AD.txt', 'r').readlines()
K279M_BBR = open('../../K279M/BBR/DSSP_K279M_BBR.txt', 'r').readlines()
K279W_AD = open('../../K279W/AD/DSSP_K279W_AD.txt', 'r').readlines()
K279W_BBR = open('../../K279W/BBR/DSSP_K279W_BBR.txt', 'r').readlines()
V287T_AD = open('../../V287T/AD/DSSP_V287T_AD.txt', 'r').readlines()
V287T_BBR = open('../../V287T/BBR/DSSP_V287T_BBR.txt', 'r').readlines()

#Declare array of zeros for num of times the residues are alpha helical and helix in general
num_struct = np.zeros([30, 10])
num_alpha_helix = np.zeros([30, 10])
per_struct = np.zeros([30, 8])
per_alpha_helix = np.zeros([30, 8])

#Seperate Characters in the string and record number that are in an alpha helix
num_alpha_helix[0][:], num_struct[0][:], per_struct[0][:], per_alpha_helix[0][:] = determine_helix(WT_AD)
num_alpha_helix[1][:], num_struct[1][:], per_struct[1][:], per_alpha_helix[1][:] = determine_helix(WT_BBR)
num_alpha_helix[2][:], num_struct[2][:], per_struct[2][:], per_alpha_helix[2][:] = determine_helix(F196A_AD)
num_alpha_helix[3][:], num_struct[3][:], per_struct[3][:], per_alpha_helix[3][:] = determine_helix(F196A_BBR)
num_alpha_helix[4][:], num_struct[4][:], per_struct[4][:], per_alpha_helix[4][:] = determine_helix(L192A_AD)
num_alpha_helix[5][:], num_struct[5][:], per_struct[5][:], per_alpha_helix[5][:] = determine_helix(L192A_BBR)
num_alpha_helix[6][:], num_struct[6][:], per_struct[6][:], per_alpha_helix[6][:] = determine_helix(L192F_AD)
num_alpha_helix[7][:], num_struct[7][:], per_struct[7][:], per_alpha_helix[7][:] = determine_helix(L192F_BBR)
num_alpha_helix[8][:], num_struct[8][:], per_struct[8][:], per_alpha_helix[8][:] = determine_helix(L192N_AD)
num_alpha_helix[9][:], num_struct[9][:], per_struct[9][:], per_alpha_helix[9][:] = determine_helix(L192N_BBR)
num_alpha_helix[10][:], num_struct[10][:], per_struct[10][:], per_alpha_helix[10][:] = determine_helix(L195A_AD)
num_alpha_helix[11][:], num_struct[11][:], per_struct[11][:], per_alpha_helix[11][:] = determine_helix(L195A_BBR)
num_alpha_helix[12][:], num_struct[12][:], per_struct[12][:], per_alpha_helix[12][:] = determine_helix(L195F_AD)
num_alpha_helix[13][:], num_struct[13][:], per_struct[13][:], per_alpha_helix[13][:] = determine_helix(L195F_BBR)
num_alpha_helix[14][:], num_struct[14][:], per_struct[14][:], per_alpha_helix[14][:] = determine_helix(L195N_AD)
num_alpha_helix[15][:], num_struct[15][:], per_struct[15][:], per_alpha_helix[15][:] = determine_helix(L195N_BBR)
num_alpha_helix[16][:], num_struct[16][:], per_struct[16][:], per_alpha_helix[16][:] = determine_helix(S286A_AD)
num_alpha_helix[17][:], num_struct[17][:], per_struct[17][:], per_alpha_helix[17][:] = determine_helix(S286A_BBR)
num_alpha_helix[18][:], num_struct[18][:], per_struct[18][:], per_alpha_helix[18][:] = determine_helix(F280Y_AD)
num_alpha_helix[19][:], num_struct[19][:], per_struct[19][:], per_alpha_helix[19][:] = determine_helix(F280Y_BBR)
num_alpha_helix[20][:], num_struct[20][:], per_struct[20][:], per_alpha_helix[20][:] = determine_helix(E276L_AD)
num_alpha_helix[21][:], num_struct[21][:], per_struct[21][:], per_alpha_helix[21][:] = determine_helix(E276L_BBR)
num_alpha_helix[22][:], num_struct[22][:], per_struct[22][:], per_alpha_helix[22][:] = determine_helix(E276F_AD)
num_alpha_helix[23][:], num_struct[23][:], per_struct[23][:], per_alpha_helix[23][:] = determine_helix(E276F_BBR)
num_alpha_helix[24][:], num_struct[24][:], per_struct[24][:], per_alpha_helix[24][:] = determine_helix(K279M_AD)
num_alpha_helix[25][:], num_struct[25][:], per_struct[25][:], per_alpha_helix[25][:] = determine_helix(K279M_BBR)
num_alpha_helix[26][:], num_struct[26][:], per_struct[26][:], per_alpha_helix[26][:] = determine_helix(K279W_AD)
num_alpha_helix[27][:], num_struct[27][:], per_struct[27][:], per_alpha_helix[27][:] = determine_helix(K279W_BBR)
num_alpha_helix[28][:], num_struct[28][:], per_struct[28][:], per_alpha_helix[28][:] = determine_helix(V287T_AD)
num_alpha_helix[29][:], num_struct[29][:], per_struct[29][:], per_alpha_helix[29][:] = determine_helix(V287T_BBR)

#Compare Ligand Bound Disordered Structures
num = [1, 2, 3, 4, 5, 6, 7, 8]
Method = ['Res 287', 'Res 288', 'Res 289', 'Res 290', 'Res 291', 'Res 292', 'Res 293', 'Res 294']
fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.set_title('Degree of Structure in AD Bound Disordered Conformations')
ax1.set_ylabel('% Residue was in alpha helix')
ax1.set_ylim(0,100)
ax1.plot(num,per_alpha_helix[0][:], color = 'red', label='WT AD')
plt.xticks(num, Method)
leg = ax1.legend()
fig.savefig('Helicity_cmpr_AD.png')
plt.close(fig)

