#Import Necessary Packages
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from scipy.stats import norm
from scipy.optimize import curve_fit

def load_rmsf(path):
    raw = open('../../' + path + '/analysis/rmsf_ref_self.txt').readlines()
    r = np.zeros(len(raw))
    for i in range(len(raw)):
        r[i] = float(raw[i])*10
    return r

def plot_two(RMSF_all, ind, Label):
    #Seperate Data for Plotting
    rmsf_0 = RMSF_all[ind[0],:]
    rmsf_1 = RMSF_all[ind[1],:]
    label_0 = Label[ind[0]]
    label_1 = Label[ind[1]]
    
    #Establish atom numbers
    res = np.linspace(1, 299, num=len(rmsf_0))
    
    #Plot
    fig = plt.figure(figsize=(18,10))
    ax1 = fig.add_subplot(111)
    ax1.set_title("Comparison of RMSF between " + str(label_0) + " and " + str(label_1))
    ax1.set_ylabel('RMSF(nm)')
    ax1.set_xlabel('Residue Number')
    ax1.plot(res, rmsf_0, label=label_0)
    ax1.plot(res, rmsf_1, label=label_1)
    ax1.legend(loc='upper right')
    fig.savefig('RMSF_compare_' + label_0 + '_' + label_1 + '_self.png') 

def plot_three(RMSF_all, ind, Label):
    #Seperate Data for Plotting
    rmsf_0 = RMSF_all[ind[0],:]
    rmsf_1 = RMSF_all[ind[1],:]
    rmsf_2 = RMSF_all[ind[2],:]
    label_0 = Label[ind[0]]
    label_1 = Label[ind[1]]
    label_2 = Label[ind[2]]

    #Establish atom numbers
    res = np.linspace(1, 299, num=len(rmsf_0))
    
    #Plot
    fig = plt.figure(figsize=(18,10))
    ax1 = fig.add_subplot(111)
    ax1.set_title("Comparison of RMSF between " + str(label_0) + " and " + str(label_1) + " and " + str(label_2))
    ax1.set_ylabel('RMSF(nm)')
    ax1.set_xlabel('Residue Number')
    ax1.plot(res, rmsf_0, label=label_0)
    ax1.plot(res, rmsf_1, label=label_1)
    ax1.plot(res, rmsf_2, label=label_2)
    ax1.legend(loc='upper right')
    fig.savefig('RMSF_compare_' + label_0 + '_' + label_1 + '_' + label_2 + '_self.png') 

def plot_three_zoom(RMSF_all, ind, Label, res_range):
    #Seperate Data for Plotting
    rmsf_0all = RMSF_all[ind[0],:]
    rmsf_1all = RMSF_all[ind[1],:]
    rmsf_2all = RMSF_all[ind[2],:]
    label_0 = Label[ind[0]]
    label_1 = Label[ind[1]]
    label_2 = Label[ind[2]]

    #Establish atom numbers
    res_all = np.linspace(1, 299, num=len(rmsf_0all))
    
    #Limit to residue range
    res, rmsf_0, rmsf_1, rmsf_2 = [],[],[],[]
    for i in range(len(res_all)):
        if res_all[i] >= res_range[0] and res_all[i] <= res_range[1]:
            rmsf_0.append(rmsf_0all[i])
            rmsf_1.append(rmsf_1all[i])
            rmsf_2.append(rmsf_2all[i])
            res.append(res_all[i])
    #Plot
    fig = plt.figure(figsize=(18,10))
    ax1 = fig.add_subplot(111)
    ax1.set_title("Comparison of RMSF between " + str(label_0) + " and " + str(label_1) + " and " + str(label_2))
    ax1.set_ylabel('RMSF(nm)')
    ax1.set_xlabel('Residue Number')
    ax1.plot(res, rmsf_0, label=label_0)
    ax1.plot(res, rmsf_1, label=label_1)
    ax1.plot(res, rmsf_2, label=label_2)
    ax1.legend(loc='upper right')
    fig.savefig('RMSF_compare_' + label_0 + '_' + label_1 + '_' + label_2 + '_' + str(res_range[0]) + '_' + str(res_range[1]) + '_self.png') 

#File paths for all input files
file_path = ['../Apo_dis', 'F196A/Apo', 'L192F/Apo', 'L195F/Apo', 'F280Y/Apo', 'E276F/Apo', 'V287T/Apo']

#Load WT
rmsf = load_rmsf(file_path[0])

#Set empty array sizes
RMSF = np.zeros((len(file_path), len(rmsf)))
RMSF[0][:] = rmsf

for i in range(1, len(file_path)):
    #Load Data
    RMSF[i][:] = load_rmsf(file_path[i])

#Name Labels
Label = ['WT', 'F196A', 'L192F', 'L195F', 'F280Y', 'E276F', 'V287T']
Label_lim = ['WT', 'F196A', 'L192F', 'L195F', 'F280Y', 'V287T'] #Remove E276F charge changing mutation

#Plot all compared to WT
for i in range(1, len(Label)):
    plot_two(RMSF, [0,i], Label)

plot_three(RMSF, [0, 6, 4], Label)
plot_three(RMSF, [0, 2, 6], Label)
plot_three(RMSF, [0, 1, 6], Label)

#Plot only the substrate binding loop
plot_three_zoom(RMSF, [0, 6, 4], Label, [113, 120])
plot_three_zoom(RMSF, [0, 2, 6], Label, [113, 120])
plot_three_zoom(RMSF, [0, 1, 6], Label, [113, 120])

#Plot only the residues 290-293
plot_three_zoom(RMSF, [0, 6, 4], Label, [289, 292])
plot_three_zoom(RMSF, [0, 2, 6], Label, [289, 292])
plot_three_zoom(RMSF, [0, 1, 6], Label, [289, 292])

#Mean RMSF values
rmsf_mean = np.zeros(len(Label))
rmsf_err = np.zeros(len(Label))
rmsf_mean_lim = np.zeros(len(Label)-1)
rmsf_err_lim = np.zeros(len(Label)-1)
n = 0
for i in range(len(Label)):
    rmsf_mean[i] = np.mean(RMSF[i][:])
    rmsf_err[i] = stats.sem(RMSF[i][:])
    if i != 5:
        rmsf_mean_lim[n] = rmsf_mean[i]
        rmsf_err_lim[n] = rmsf_err[i]
        n += 1

#Plot helicity vs dG for AD
dG = [0.00, -0.437, 0.375, 0.044, -0.634, 1.526, -0.056]
dG_MBAR = [0.00, -0.381, 0.297, 0.019, -0.236, 0.373, 0.026]
dG_MBAR_err = [0.00, 0.109, 0.095, 0.074, 0.072, 0.694, 0.106]

fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.set_title(r'Comparison of RMSF of the $\alpha$-7 Helix')    
ax1.set_ylabel(r'$\Delta\Delta$G from WT')
ax1.set_xlabel(r'RMSF($\AA$)')
ax1.scatter(rmsf_mean, dG_MBAR, label = 'AD', color = 'blue')
plt.errorbar(rmsf_mean, dG_MBAR, xerr= rmsf_err, yerr = dG_MBAR_err, fmt='o', color='blue')
fig.savefig('RMSF_v_dG.png')
plt.close(fig)

#Remove E276F and add gausian approximation line
dG = [0.00, -0.437, 0.375, 0.044, -0.634, -0.056]
dG_MBAR = [0.00, -0.381, 0.297, 0.019, -0.236, 0.026]
dG_MBAR_err = [0.00, 0.109, 0.095, 0.074, 0.072, 0.106]

# Define the Gaussian function
def Gauss(x, a, x0, sigma, b):
    return a*np.exp(-(x-x0)**2/(2*sigma**2)) + b
parameters, covariance = curve_fit(Gauss, rmsf_mean_lim, dG_MBAR, p0 = [0.4, 0.87, 0.1, 0])
  
[fit_A, fit_B, fit_C, fit_D] = parameters

x = np.linspace(0, 2, num=100)

#Plot
fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.set_title(r'Comparison of RMSF of the $\alpha$-7 Helix and Binding Affinity') 
ax1.set_ylabel(r'$\Delta\Delta$G from WT')
ax1.set_xlabel(r'RMSF($\AA$)')
ax1.scatter(rmsf_mean_lim, dG_MBAR, color = 'blue', label = 'Simulation Data')
plt.errorbar(rmsf_mean_lim, dG_MBAR, xerr= rmsf_err_lim, yerr = dG_MBAR_err, fmt='o', color='blue')
plt.plot(x, Gauss(x, fit_A, fit_B, fit_C, fit_D), label = 'Gaussian Fit')
# Loop for annotation of all points
for i in range(len(Label_lim)):
    if i == 0:
        plt.annotate(Label_lim[i], (rmsf_mean_lim[i], dG_MBAR[i] + 0.05))
    elif i == 5:
        plt.annotate(Label_lim[i], (rmsf_mean_lim[i]-0.1, dG_MBAR[i] - 0.1))
    else:
        plt.annotate(Label_lim[i], (rmsf_mean_lim[i], dG_MBAR[i] - 0.1))
plt.legend(loc='best')
fig.savefig('RMSF_v_dG_lim.png')
plt.close(fig)

#Plot DSSP against relative inhibition
F_BBR = [0.00, -0.259789014, -0.303991355, -0.488133372, -0.447825789, -0.243071805, 0.164257023]
F_BBR_se = [0.00, 0.001039749, 0.001238509, 0.001227868, 0.001397659, 0.001188003, 0.000777901]
F_AD = [0.00, -0.204972546, -0.035556769, -0.2663284, -0.286222165, -0.331721332, 0.131296446]
F_AD_se = [0.00, 0.000911951, 0.001351002, 0.000930425, 0.001033858, 0.001059996, 0.000928346]

#Make plot
fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.set_title(r'Comparison of RMSF of the $\alpha$-7')    
ax1.set_xlabel(r'RMSF of the $\alpha$-7 Helix')
ax1.set_ylabel('Relative Inhibition (F)')
ax1.scatter(rmsf_mean, F_AD, label = 'AD', color = 'blue')
plt.errorbar(rmsf_mean, F_AD, xerr= rmsf_err, yerr = F_AD, fmt='o', color='blue')
#ax1.scatter(rmsf_mean, F_BBR, label = 'BBR', color = 'purple')
#plt.errorbar(rmsf_mean, F_BBR, xerr= rmsf_err, yerr = F_BBR, fmt='o', color='purple')
fig.savefig('rmsf_v_F.png')
plt.close(fig)


