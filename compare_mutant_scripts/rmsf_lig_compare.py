#Import Necessary Packages
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

def load_data(file_path):
    with open('../../' + file_path + '/analysis/rmsf_lig_' + file_name + '.xvg') as f:
        for _ in range(17):
            next(f)
        for line in f:
            cols = line.split()
            if len(cols) == 2:
                t.append(float(cols[0]))
                r.append(float(cols[1]))
    return t, r

def adjacent_values(vals, q1, q3):
    upper_adjacent_value = q3 + (q3 - q1) * 1.5
    upper_adjacent_value = np.clip(upper_adjacent_value, q3, vals[-1])

    lower_adjacent_value = q1 - (q3 - q1) * 1.5
    lower_adjacent_value = np.clip(lower_adjacent_value, vals[0], q1)
    return lower_adjacent_value, upper_adjacent_value

def set_axis_style(ax, labels):
    ax.xaxis.set_tick_params(direction='out')
    ax.xaxis.set_ticks_position('bottom')
    ax.set_xticks(np.arange(1, len(labels) + 1))
    ax.set_xticklabels(labels)
    ax.set_xlim(0.25, len(labels) + 0.75)
    ax.set_xlabel('Sample name')

#Load data
File_paths_AD = ['../AD_dis/analysis/config11', 'F196A/AD']

fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.set_title("Comparison of a7 helix RMSF")    
ax1.set_xlabel('Residues')
ax1.set_ylabel('RMSF(nm)')
ax1.plot(res_a7, a7_Apo_open, label='Apo Open', color = 'blue')
ax1.plot(res_a7, a7_AD_open, label='AD', color = 'blue', linestyle = 'dashed')
ax1.plot(res_a7, a7_BBR_open, label='BBR', color = 'blue', linestyle = 'dotted')
ax1.plot(res_a7, a7_Apo_close, label='Apo Closed', color = 'red')
ax1.legend(loc='best')
fig.savefig('a7_RMSF_compare.png') 

