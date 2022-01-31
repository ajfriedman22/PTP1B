from matplotlib import pyplot as plt
import numpy as np

#Import files for state occupation
loc1_1sug_AD_dis11, loc2_1sug_AD_dis11 = [],[]
for i in open('../../1sug_dis_AD/config11/AD_bound_crys_v_time_strict.txt', 'r').readlines():
    loc1_1sug_AD_dis11.append(float(i))
for i in open('../../1sug_dis_AD/config11/AD_bound_alt3_v_time_strict.txt', 'r').readlines():
    loc2_1sug_AD_dis11.append(float(i))
loc1_1sug_AD_alt, loc2_1sug_AD_alt = [],[]
for i in open('../../1sug_AD_dis_alt/run_2/AD_bound_crys_v_time_strict.txt', 'r').readlines():
    loc1_1sug_AD_alt.append(float(i))
for i in open('../../1sug_AD_dis_alt/run_2/AD_bound_alt3_v_time_strict.txt', 'r').readlines():
    loc2_1sug_AD_alt.append(float(i))

#Determine the time that transitions occur
time_loc1_loc2_dis11, time_loc2_loc1_dis11 = [],[]
time_loc1_loc2_alt, time_loc2_loc1_alt = [],[]
time_convert = len(loc1_1sug_AD_dis11)/(300-5) #convert time to ns
time_convert_alt = len(loc1_1sug_AD_alt)/(300-60) #convert time to ns
for i in range(len(loc1_1sug_AD_dis11) - 1):
    if loc1_1sug_AD_dis11[i] == 1 and loc1_1sug_AD_dis11[i+1] == 0 and loc2_1sug_AD_dis11[i] == 0 and loc2_1sug_AD_dis11[i+1] == 1:
        time_loc1_loc2_dis11.append(i/time_convert)
    if loc1_1sug_AD_dis11[i] == 0 and loc1_1sug_AD_dis11[i+1] == 1 and loc2_1sug_AD_dis11[i] == 1 and loc2_1sug_AD_dis11[i+1] == 0:
        time_loc2_loc1_dis11.append(i/time_convert)
for i in range(len(loc1_1sug_AD_alt) - 1):
    if loc1_1sug_AD_alt[i] == 1 and loc1_1sug_AD_alt[i+1] == 0 and loc2_1sug_AD_alt[i] == 0 and loc2_1sug_AD_alt[i+1] == 1:
        time_loc1_loc2_alt.append(i/time_convert_alt)
    if loc1_1sug_AD_alt[i] == 0 and loc1_1sug_AD_alt[i+1] == 1 and loc2_1sug_AD_alt[i] == 1 and loc2_1sug_AD_alt[i+1] == 0:
        time_loc2_loc1_alt.append(i/time_convert_alt)

#Determine minimum time between transitions
trans_time_dis11, trans_time_alt = [],[]
for i in range(len(time_loc1_loc2_dis11)):
    trans_time_dis11.append(abs(time_loc1_loc2_dis11[i] - time_loc2_loc1_dis11[i]))
for i in range(len(time_loc1_loc2_alt)):
    trans_time_alt.append(abs(time_loc1_loc2_alt[i] - time_loc2_loc1_alt[i]))

#Print the number of transitions and minimum time between them
file_output = open('Transitions.txt', 'w')
file_output.write('Dis11\n')
file_output.write('Number of transitions loc1 to loc2: ' + str(len(time_loc1_loc2_dis11)) + '\n')
file_output.write('Number of transitions loc2 to loc1: ' + str(len(time_loc2_loc1_dis11)) + '\n')
file_output.write('Minimum time b/w transitions: ' + str(min(trans_time_dis11)) + ' ns \n')
file_output.write('Mean time b/w transitions: ' + str(np.mean(trans_time_dis11)) + ' ns \n')
file_output.write('Alt\n')
file_output.write('Number of transitions loc1 to loc2: ' + str(len(time_loc1_loc2_alt)) + '\n')
file_output.write('Number of transitions loc2 to loc1: ' + str(len(time_loc2_loc1_alt)) + '\n')
file_output.write('Minimum time b/w transitions: ' + str(min(trans_time_alt)) + ' ns\n')
file_output.write('Mean time b/w transitions: ' + str(np.mean(trans_time_alt)) + ' ns')

