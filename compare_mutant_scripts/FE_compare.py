#Import Necessary Packages
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as lines
from scipy import stats
from itertools import product
import seaborn as sns

#Record mutant FE differences
Mutants = ['F196A', 'L192F', 'L195F', 'E276F', 'V287T', 'F280Y']
AD_RFE_TI = [-0.473, 0.375, 0.044, ?, -0.056, -0.634]
AD_RFE_TI_err = [0.133, 0.103, 0.086, ?, 0.186, 0.134]
AD_RFE_MBAR = [-0.381, 0.297, 0.019, ?, 0.026, -0.236]
AD_RFE_MBAR_err = [0.109, 0.095, 0.074, ?, 0.106, 0.072]

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
ax1.bar(num, AD_RFE_TI, width = 0.22)
plt.errorbar(num, AD_RFE_TI, yerr = AD_RFE_MBAR_err)
ax1.bar(num2, AD_RFE_MBAR, width = 0.22)
plt.errorbar(num2, AD_RFE_MBAR, yerr = AD_RFE_MBAR_err)
ax2=ax1.twinx()
ax2.bar(num3, AD_exp_low, width = 0.22)
ax2.errorbar(num3, AD_exp_low, yerr = AD_exp_low_err)
ax2.bar(num4, AD_exp_high, width = 0.22)
ax2.errorbar(num4, AD_exp_high, yerr = AD_exp_high_err)
ax2.set_ylabel('Relative PTP1B Inhibiton (F)')
plt.xticks(num, Mutants, fontsize=8)
fig.savefig('cmpr_FE_exp_inhibition_AD.png')
plt.close(fig)
