import mdtraj as md
import numpy as np
import argparse

#Declare arguments
parser = argparse.ArgumentParser(description = 'Determination of RMSD of the beta sheet, a3, a6, and a7 helices')
parser.add_argument('-t', required=True, help='File name for input trajectory')
parser.add_argument('-g', required=True, help= 'File name for input topology (gro format)')
#parser.add_argument('-r', required=True, help='Reference Structure')
parser.add_argument('-f', required=True, help='File name base')
parser.add_argument('-s', required=False, default=False , type = bool, help='True is the base is 1sug')

#Import Arguments
args = parser.parse_args()
File_traj = args.t + '.xtc'
File_gro = args.g + '.gro'
#File_ref = args.r + '.pdb'
File_ref = '../../1sug/Cluster_helix/1sug_cluster_PTP1B.pdb'
File_base = args.f
check_s = args.s

#Load trajectory and reference
traj = md.load(File_traj, top=File_gro)
top = traj.topology
traj_bb = traj.atom_slice(top.select('protein and backbone'))
ref_full = md.load(File_ref)
top_ref = ref_full.topology
ref = ref_full.atom_slice(top_ref.select('protein and backbone'))

#Allign to reference structure
if check_s == True:
    traj_bb.superpose(ref)
else:
    all_but1 = np.linspace(4, 1195, num=1192, dtype=int)
    all_ = np.linspace(0, 1191, num=1192, dtype=int)

    traj_bb.superpose(ref, atom_indices = all_but1, ref_atom_indices=all_)

#Determine atom indices for structure and reference
atom_ref_beta = np.linspace(548, 708, num=160) #residue 139 to 178
atom_ref_L11 = np.linspace(593, 608, num=16)#residue 150 to 153)
atom_ref_a3 = np.linspace(737, 796, num=59) #residue 186 to 200
atom_ref_a6 = np.linspace(1049, 1120, num=71) #residue 264 to 281
atom_ref_a7 = np.linspace(1132, 1188, num=61) #residue 285 to 297
#Arrays for uncorrelated samples
beta_uncorr, L11_uncorr, a3_uncorr, a6_uncorr, a7_uncorr = [],[],[],[],[]
if check_s == True:
    #Compute RMSD for B-sheet
    rmsd_beta = md.rmsd(traj_bb, ref, atom_indices=atom_ref_beta)
    #Determine Maximum and Minimum indices
    max_ind, min_ind, diff = [],[],[]
    for i in range(len(rmsd_beta)-1):
        if rmsd_beta[i+1] > rmsd_beta[i] and rmsd_beta[i-1] > rmsd_beta[i]:
            min_ind.append(i)
        if rmsd_beta[i+1] < rmsd_beta[i] and rmsd_beta[i-1] < rmsd_beta[i]:
            max_ind.append(i)
    #Compute Approx Correlation Time for B-sheet
    for i in range(min([len(min_ind), len(max_ind)])):
        diff.append(abs(max_ind[i]-min_ind[i]))
    corr_time = max(diff) * 2
    for i in range(len(rmsd_beta)):
        if i%corr_time == 0:
            beta_uncorr.append(rmsd_beta[i])
    np.savetxt('beta_' + File_base + '_rmsd.txt', beta_uncorr)

    #Compute RMSD for L11 loop
    rmsd_L11 = md.rmsd(traj_bb, ref, atom_indices=atom_ref_L11)
    #Determine Maximum and Minimum indices
    max_ind, min_ind, diff = [],[],[]
    for i in range(len(rmsd_L11)-1):
        if rmsd_L11[i+1] > rmsd_L11[i] and rmsd_L11[i-1] > rmsd_L11[i]:
            min_ind.append(i)
        if rmsd_L11[i+1] < rmsd_L11[i] and rmsd_L11[i-1] < rmsd_L11[i]:
            max_ind.append(i)
    #Compute Approx Correlation Time for B-sheet
    for i in range(min([len(min_ind), len(max_ind)])):
        diff.append(abs(max_ind[i]-min_ind[i]))
    corr_time = max(diff) * 2
    for i in range(len(rmsd_L11)):
        if i%corr_time == 0:
            L11_uncorr.append(rmsd_L11[i])
    np.savetxt('L11_' + File_base + '_rmsd.txt', L11_uncorr)

    #Compute RMSD for a3 helix
    rmsd_a3 = md.rmsd(traj_bb, ref, atom_indices=atom_ref_a3)
    #Determine Maximum and Minimum indices
    max_ind, min_ind, diff = [],[],[]
    for i in range(len(rmsd_a3)-1):
        if rmsd_a3[i+1] > rmsd_a3[i] and rmsd_a3[i-1] > rmsd_a3[i]:
            min_ind.append(i)
        if rmsd_a3[i+1] < rmsd_a3[i] and rmsd_a3[i-1] < rmsd_a3[i]:
            max_ind.append(i)
    #Compute Approx Correlation Time for B-sheet
    for i in range(min([len(min_ind), len(max_ind)])):
        diff.append(abs(max_ind[i]-min_ind[i]))
    corr_time = max(diff) * 2
    for i in range(len(rmsd_a3)):
        if i%corr_time == 0:
            a3_uncorr.append(rmsd_a3[i])
    np.savetxt('a3_' + File_base + '_rmsd.txt', a3_uncorr)

    #Compute RMSD for a6 helix
    rmsd_a6 = md.rmsd(traj_bb, ref, atom_indices=atom_ref_a6)
    #Determine Maximum and Minimum indices
    max_ind, min_ind, diff = [],[],[]
    for i in range(len(rmsd_a6)-1):
        if rmsd_a6[i+1] > rmsd_a6[i] and rmsd_a6[i-1] > rmsd_a6[i]:
            min_ind.append(i)
        if rmsd_a6[i+1] < rmsd_a6[i] and rmsd_a6[i-1] < rmsd_a6[i]:
            max_ind.append(i)
    #Compute Approx Correlation Time for B-sheet
    for i in range(min([len(min_ind), len(max_ind)])):
        diff.append(abs(max_ind[i]-min_ind[i]))
    corr_time = max(diff) * 2
    for i in range(len(rmsd_a6)):
        if i%corr_time == 0:
            a6_uncorr.append(rmsd_a6[i])
    np.savetxt('a6_' + File_base + '_rmsd.txt', a6_uncorr)
   
    #Compute RMSD for a3 helix
    rmsd_a7 = md.rmsd(traj_bb, ref, atom_indices=atom_ref_a7)
    #Determine Maximum and Minimum indices
    max_ind, min_ind, diff = [],[],[]
    for i in range(len(rmsd_a7)-1):
        if rmsd_a7[i+1] > rmsd_a7[i] and rmsd_a7[i-1] > rmsd_a7[i]:
            min_ind.append(i)
        if rmsd_a7[i+1] < rmsd_a7[i] and rmsd_a7[i-1] < rmsd_a7[i]:
            max_ind.append(i)
    #Compute Approx Correlation Time for B-sheet
    for i in range(min([len(min_ind), len(max_ind)])):
        diff.append(abs(max_ind[i]-min_ind[i]))
    corr_time = max(diff) * 2
    for i in range(len(rmsd_a7)):
        if i%corr_time == 0:
            a7_uncorr.append(rmsd_a7[i])
    np.savetxt('a7_' + File_base + '_rmsd.txt', a7_uncorr)

if check_s == False:
    atom_beta = np.linspace(552, 712, num=160) #residue 139 to 178
    atom_L11 = np.linspace(597, 612, num=16)#residue 150 to 153)
    atom_a3 = np.linspace(741, 780, num=59) #residue 186 to 200
    atom_a6 = np.linspace(1053, 1124, num=71) #residue 264 to 281
    atom_a7 = np.linspace(1136, 1192, num=61) #residue 285 to 297

    #Compute RMSD for B-sheet
    rmsd_beta = md.rmsd(traj_bb, ref, atom_indices = atom_beta, ref_atom_indices=atom_ref_beta)
    #Determine Maximum and Minimum indices
    max_ind, min_ind, diff = [],[],[]
    for i in range(len(rmsd_beta)-1):
        if rmsd_beta[i+1] > rmsd_beta[i] and rmsd_beta[i-1] > rmsd_beta[i]:
            min_ind.append(i)
        if rmsd_beta[i+1] < rmsd_beta[i] and rmsd_beta[i-1] < rmsd_beta[i]:
            max_ind.append(i)
    #Compute Approx Correlation Time for B-sheet
    for i in range(min([len(min_ind), len(max_ind)])):
        diff.append(abs(max_ind[i]-min_ind[i]))
    corr_time = max(diff) * 2
    for i in range(len(rmsd_beta)):
        if i%corr_time == 0:
            beta_uncorr.append(rmsd_beta[i])
    np.savetxt('beta_' + File_base + '_rmsd.txt', beta_uncorr)
    
    #Compute RMSD for L11 loop
    rmsd_L11 = md.rmsd(traj_bb, ref, atom_indices = atom_L11, ref_atom_indices=atom_ref_L11)
    #Determine Maximum and Minimum indices
    max_ind, min_ind, diff = [],[],[]
    for i in range(len(rmsd_L11)-1):
        if rmsd_L11[i+1] > rmsd_L11[i] and rmsd_L11[i-1] > rmsd_L11[i]:
            min_ind.append(i)
        if rmsd_L11[i+1] < rmsd_L11[i] and rmsd_L11[i-1] < rmsd_L11[i]:
            max_ind.append(i)
    #Compute Approx Correlation Time for B-sheet
    for i in range(min([len(min_ind), len(max_ind)])):
        diff.append(abs(max_ind[i]-min_ind[i]))
    corr_time = max(diff) * 2
    for i in range(len(rmsd_L11)):
        if i%corr_time == 0:
            L11_uncorr.append(rmsd_L11[i])
    np.savetxt('L11_' + File_base + '_rmsd.txt', L11_uncorr)

    #Compute RMSD for a3 helix
    rmsd_a3 = md.rmsd(traj_bb, ref, atom_indices = atom_a3, ref_atom_indices=atom_ref_a3)
    #Determine Maximum and Minimum indices
    max_ind, min_ind, diff = [],[],[]
    for i in range(len(rmsd_a3)-1):
        if rmsd_a3[i+1] > rmsd_a3[i] and rmsd_a3[i-1] > rmsd_a3[i]:
            min_ind.append(i)
        if rmsd_a3[i+1] < rmsd_a3[i] and rmsd_a3[i-1] < rmsd_a3[i]:
            max_ind.append(i)
    #Compute Approx Correlation Time for B-sheet
    for i in range(min([len(min_ind), len(max_ind)])):
        diff.append(abs(max_ind[i]-min_ind[i]))
    corr_time = max(diff) * 2
    for i in range(len(rmsd_a3)):
        if i%corr_time == 0:
            a3_uncorr.append(rmsd_a3[i])
    np.savetxt('a3_' + File_base + '_rmsd.txt', a3_uncorr)

    #Compute RMSD for a6 helix
    rmsd_a6 = md.rmsd(traj_bb, ref, atom_indices = atom_a6, ref_atom_indices=atom_ref_a6)
    #Determine Maximum and Minimum indices
    max_ind, min_ind, diff = [],[],[]
    for i in range(len(rmsd_a6)-1):
        if rmsd_a6[i+1] > rmsd_a6[i] and rmsd_a6[i-1] > rmsd_a6[i]:
            min_ind.append(i)
        if rmsd_a6[i+1] < rmsd_a6[i] and rmsd_a6[i-1] < rmsd_a6[i]:
            max_ind.append(i)
    #Compute Approx Correlation Time for B-sheet
    for i in range(min([len(min_ind), len(max_ind)])):
        diff.append(abs(max_ind[i]-min_ind[i]))
    corr_time = max(diff) * 2
    for i in range(len(rmsd_a6)):
        if i%corr_time == 0:
            a6_uncorr.append(rmsd_a6[i])
    np.savetxt('a6_' + File_base + '_rmsd.txt', a6_uncorr)
   
    #Compute RMSD for a3 helix
    rmsd_a7 = md.rmsd(traj_bb, ref, atom_indices = atom_a7, ref_atom_indices=atom_ref_a7)
    #Determine Maximum and Minimum indices
    max_ind, min_ind, diff = [],[],[]
    for i in range(len(rmsd_a7)-1):
        if rmsd_a7[i+1] > rmsd_a7[i] and rmsd_a7[i-1] > rmsd_a7[i]:
            min_ind.append(i)
        if rmsd_a7[i+1] < rmsd_a7[i] and rmsd_a7[i-1] < rmsd_a7[i]:
            max_ind.append(i)
    #Compute Approx Correlation Time for B-sheet
    for i in range(min([len(min_ind), len(max_ind)])):
        diff.append(abs(max_ind[i]-min_ind[i]))
    corr_time = max(diff) * 2
    for i in range(len(rmsd_a7)):
        if i%corr_time == 0:
            a7_uncorr.append(rmsd_a7[i])
    np.savetxt('a7_' + File_base + '_rmsd.txt', a7_uncorr)


