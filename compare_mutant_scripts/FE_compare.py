#Import Necessary Packages
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as lines
from scipy import stats
from itertools import product
import seaborn as sns

#Record mutant FE differences
Mutants = ['F196A', 'L192F', 'L195F', 'E276F', 'V287T', 'F280Y']
AD_RFE_TI = [-0.473, 0.375, 0.044, 0, -0.056, -0.634]
AD_RFE_TI_err = [0.133, 0.103, 0.086, 0, 0.186, 0.134]
AD_RFE_MBAR = [-0.381, 0.297, 0.019, 0, 0.026, -0.236]
AD_RFE_MBAR_err = [0.109, 0.095, 0.074, 0, 0.106, 0.072]

AD_exp_low = [-0.050030578, -0.00971087, -0.125098767, -0.097286637, -0.114041226, 0.153088661]
AD_exp_high = [-0.204972546, -0.035556769, -0.2663284, -0.286222165, -0.331721332, 0.131296446]
AD_exp_low_err = [0.000981403, 0.00078874, 0.000724636, 0.000701781, 0.000815448, 0.000586185]
AD_exp_high_err = [0.000911951, 0.001351002, 0.000930425, 0.001033858, 0.001059996, 0.000928346]

num = np.linspace(0, 5, num = 6)
num2 = np.linspace(0.25, 5.25, num = 6)
num3 = np.linspace(0.5, 5.5, num = 6)
num4 = np.linspace(0.75, 5.75, num = 6)

fig,ax1 = plt.subplots()
ax1.set_title('Comparison of PTP1B Inhibition and RFE') 
ax1.set_ylabel('Relative Binding FE (kcal/mol)')
ax1.set_ylim(-0.7, 0.7)
ax1.bar(num, AD_RFE_TI, width = 0.22, color = 'blue')
#ax1.errorbar(num, AD_RFE_TI, yerr = AD_RFE_MBAR_err, color = 'blue')
ax1.bar(num2, AD_RFE_MBAR, width = 0.22, color = 'purple')
#ax1.errorbar(num2, AD_RFE_MBAR, yerr = AD_RFE_MBAR_err, color = 'purple')
ax2=ax1.twinx()
ax2.set_ylim(-0.4, 0.4)
ax2.bar(num3, AD_exp_low, width = 0.22, color = 'red')
#ax2.errorbar(num3, AD_exp_low, yerr = AD_exp_low_err, color = 'red')
ax2.bar(num4, AD_exp_high, width = 0.22, color = 'orange')
#ax2.errorbar(num4, AD_exp_high, yerr = AD_exp_high_err, color = 'orange')
ax2.set_ylabel('Relative PTP1B Inhibiton (F)')
plt.xticks(num3, Mutants, fontsize=8)
fig.savefig('cmpr_FE_exp_inhibition_AD_alt.png')
plt.close(fig)

#Alt style
num_mut = np.linspace(0.5, 5.5, num = 6)

Label = ['TI Estimate', 'MBAR Estimate']
Label_exp = [r'EXP 77 $\mu$m', r'EXP 115.5 $\mu$m']
colors = ['blue', 'purple']
colors_exp = ['red', 'orange']

fig,ax1 = plt.subplots()
ax1.set_title('Comparison of PTP1B Inhibition and RFE') 
ax1.set_ylabel('Relative Binding FE (kcal/mol)')
ax1.set_ylim(-0.7, 0.7)
ax1.bar([0, 0.25], [AD_RFE_TI[0], AD_RFE_MBAR[0]], width = 0.22, color = colors, label = Label)
plt.errorbar([0, 0.25], [AD_RFE_TI[0], AD_RFE_MBAR[0]], yerr = [AD_RFE_TI_err[0], AD_RFE_MBAR_err[0]], color = colors)
ax1.bar([1, 1.25], [AD_RFE_TI[1], AD_RFE_MBAR[1]], width = 0.22, color = colors)
plt.errorbar([1, 1.25], [AD_RFE_TI[1], AD_RFE_MBAR[1]], yerr = [AD_RFE_TI_err[1], AD_RFE_MBAR_err[1]], color = colors)
ax1.bar([2, 2.25], [AD_RFE_TI[2], AD_RFE_MBAR[2]], width = 0.22, color = colors)
plt.errorbar([2, 2.25], [AD_RFE_TI[2], AD_RFE_MBAR[2]], yerr = [AD_RFE_TI_err[2], AD_RFE_MBAR_err[2]], color = colors)
ax1.bar([3, 3.25], [AD_RFE_TI[3], AD_RFE_MBAR[3]], width = 0.22, color = colors)
plt.errorbar([3, 3.25], [AD_RFE_TI[3], AD_RFE_MBAR[3]], yerr = [AD_RFE_TI_err[3], AD_RFE_MBAR_err[3]], color = colors)
ax1.bar([4, 4.25], [AD_RFE_TI[4], AD_RFE_MBAR[4]], width = 0.22, color = colors)
plt.errorbar([4, 4.25], [AD_RFE_TI[4], AD_RFE_MBAR[4]], yerr = [AD_RFE_TI_err[4], AD_RFE_MBAR_err[4]], color = colors)
ax1.bar([5, 5.25], [AD_RFE_TI[5], AD_RFE_MBAR[5]], width = 0.22, color = colors)
plt.errorbar([5, 5.25], [AD_RFE_TI[5], AD_RFE_MBAR[5]], yerr = [AD_RFE_TI_err[5], AD_RFE_MBAR_err[5]], color = colors)

ax2.set_ylim(-0.4, 0.4)
ax2.bar([0.5, 0.75], [AD_exp_low[0], AD_exp_high[0]], width = 0.22, color = colors_exp, label = Label_exp)
ax2.errorbar([0.5, 0.75], [AD_exp_low[0], AD_exp_high[0]], yerr = [AD_exp_low_err[0], AD_exp_high_err[0]], color = colors_exp)
ax2.bar([1.5, 1.75], [AD_exp_low[1], AD_exp_high[1]], width = 0.22, color = colors_exp)
ax2.errorbar([1.5, 1.75], [AD_exp_low[1], AD_exp_high[1]], yerr = [AD_exp_low_err[1], AD_exp_high_err[1]], color = colors_exp)
ax2.bar([2.5, 2.75], [AD_exp_low[2], AD_exp_high[2]], width = 0.22, color = colors_exp)
ax2.errorbar([2.5, 2.75], [AD_exp_low[2], AD_exp_high[2]], yerr = [AD_exp_low_err[2], AD_exp_high_err[2]], color = colors_exp)
ax2.bar([3.5, 3.75], [AD_exp_low[3], AD_exp_high[3]], width = 0.22, color = colors_exp)
ax2.errorbar([3.5, 3.75], [AD_exp_low[3], AD_exp_high[3]], yerr = [AD_exp_low_err[3], AD_exp_high_err[3]], color = colors_exp)
ax2.bar([4.5, 4.75], [AD_exp_low[4], AD_exp_high[4]], width = 0.22, color = colors_exp)
ax2.errorbar([4.5, 4.75], [AD_exp_low[4], AD_exp_high[4]], yerr = [AD_exp_low_err[4], AD_exp_high_err[4]], color = colors_exp)
ax2.bar([5.5, 5.75], [AD_exp_low[5], AD_exp_high[5]], width = 0.22, color = colors_exp)
ax2.errorbar([5.5, 5.75], [AD_exp_low[5], AD_exp_high[5]], yerr = [AD_exp_low_err[5], AD_exp_high_err[5]], color = colors_exp)
ax2.set_ylabel('Relative PTP1B Inhibiton (F)')
plt.xticks(num_mut, Mutants, fontsize=8)
#plt.legend()
fig.savefig('cmpr_FE_exp_inhibition_AD.png')
plt.close(fig)

