#!/ usr / bin / env python

import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
from sklearn.decomposition import PCA
from itertools import combinations
import argparse
from itertools import product
from statistics import stdev
import sys

#Import custom modules
sys.path.insert(1,'/ocean/projects/cts160011p/afriedma/code/PTP1B/util')
import mdfunc
import uncorr

#Declare arguments
parser = argparse.ArgumentParser(description = 'Determination of DSSP, H-bonds, Ligand Contacts, Helical interactions and PCA for GROMACS Trajectory of PTP1B')
parser.add_argument('-t', required=True, help='File name for input trajectory')
parser.add_argument('-g', required=True, help= 'File name for input topology (gro format)')
parser.add_argument('-f', required=True, help='Base for all output file names')
parser.add_argument('-r', required=False, help='Should the reference structure for RMSD and RMSF be Apo Open or Closed or other?')
parser.add_argument('-rn', required=False, help='Reference name for RMSD')
parser.add_argument('-l', required=False, default='none', help='If ligand analysis should be performed, which ligand is present?')
parser.add_argument('-a', required=False, default=False, type=bool, help='Should DSSP be calculated?')
parser.add_argument('-b', required=False, default=False, type=bool, help='Should Hbond Analysis be preformed?')
parser.add_argument('-i', required=False, default=False, type=bool, help='Should helical interactions be measured?')
parser.add_argument('-p', required=False, default=False, type=bool, help='Should PCA be preformed?')
parser.add_argument('-rms', required=False, default=False, type=bool, help='Should RMSF and RMSD analysis be computed?')

#Import Arguments
args = parser.parse_args()
File_traj = args.t + '.xtc'
File_gro = args.g + '.gro'
File_base = args.f
ref_type = args.r
dssp_check = args.a
hbond_check = args.b
check_hel = args.i
pca_ck = args.p
rms_chk = args.rms
lig = args.l
if lig != 'none':
    lig_check = True
else:
    lig_check = False

#Load trajectories
traj, traj_bb, traj_prot, traj_ns = mdfunc.mdtraj_load(File_traj, File_gro) 

#Load reference trajectory
if rms_chk == True: 
    if ref_type == 'Apo_open':
        ref = '/ocean/projects/cts160011p/afriedma/PTP1B/Apo_dis/analysis/bb_cluster/Apo_dis_bb_cluster.pdb'
        ref_name = 'open'
    elif ref_type == 'Apo_closed':
        ref = '/ocean/projects/cts160011p/afriedma/PTP1B/Apo_1SUG/analysis/1sug/bb_cluster/1sug_bb_cluster.pdb'
        ref_name = 'closed'
    else:
        ref = ref_type
        ref_name = args.rn
        print(ref)
    ref_pdb = md.load_pdb(ref)
    top_ref = ref_pdb.topology
    ref_bb = ref_pdb.atom_slice(top_ref.select('backbone'))
print('Topology Loaded')

if rms_chk == True:
    if ref_type == 'Apo_open' or ref_name == 'self':
        #Calculate RMSF from reference structure
        rmsf_data = md.rmsf(traj_bb, ref_bb, parallel=True, precentered=False)

        np.savetxt('rmsf_ref_' + str(ref_name) + '.txt', rmsf_data)

        #Calculate RMSD for full protein
        rmsd_full = md.rmsd(traj_bb, ref_bb, parallel=True, precentered=False)
        t_full = uncorr.ind(rmsd_full)
        rmsd_full_uncorr = uncorr.sort(rmsd_full, t_full)

        np.savetxt('rmsd_full_ref_' + str(ref_name) + '.txt', rmsd_full_uncorr)
        
    #Set Topology
    top_bb = traj_bb.topology
    top_ref_bb = ref_bb.topology
    if ref_type != 'Apo_closed':
        ref_WPD = ref_bb.atom_slice(top_ref_bb.select('176 <= resid and resid <= 184' )) #Limit trajectory to the WPD loop of PTP1B only
        ref_WPD_a3 = ref_bb.atom_slice(top_ref_bb.select('184 <= resid and resid <= 187' )) #Limit trajectory
        ref_P = ref_bb.atom_slice(top_ref_bb.select('213 <= resid and resid <= 222' )) #Limit trajectory to the P loop of PTP1B only
        ref_CYS = ref_bb.atom_slice(top_ref_bb.select('resid == 214' )) #Limit trajectory to CYS215 of PTP1B only
        ref_SBL = ref_bb.atom_slice(top_ref_bb.select('112 <= resid and resid <= 119' )) #Limit trajectory to the S loop of PTP1B only
        ref_a3 = ref_bb.atom_slice(top_ref_bb.select('185 <= resid and resid <= 199' ))
        ref_a3_top = ref_bb.atom_slice(top_ref_bb.select('185 <= resid and resid <= 190' ))
        ref_a4 = ref_bb.atom_slice(top_ref_bb.select('220 <= resid and resid <= 237' ))
        ref_a4_P = ref_bb.atom_slice(top_ref_bb.select('218 <= resid and resid <= 227' ))
        ref_a5 = ref_bb.atom_slice(top_ref_bb.select('244 <= resid and resid <= 251' ))
        ref_a6 = ref_bb.atom_slice(top_ref_bb.select('263 <= resid and resid <= 280' ))
        ref_a6_bot = ref_bb.atom_slice(top_ref_bb.select('274 <= resid and resid <= 280' ))
        if traj_ns.n_residues > 297: #Only compute these distances if the a7 helix is present
            ref_a7 = ref_bb.atom_slice(top_ref_bb.select('286 <= resid and resid <= 294' ))
        ref_L11 = ref_bb.atom_slice(top_ref_bb.select('150 <= resid and resid <= 153' ))
        ref_Q = ref_bb.atom_slice(top_ref_bb.select('258 <= resid and resid <= 262' ))
        ref_beg = ref_bb.atom_slice(top_ref_bb.select('26 <= resid and resid <= 35' ))
    else:
        ref_WPD = ref_bb.atom_slice(top_ref_bb.select('175 <= resid and resid <= 183' )) #Limit trajectory to the WPD loop of PTP1B only
        ref_WPD_a3 = ref_bb.atom_slice(top_ref_bb.select('183 <= resid and resid <= 186' )) #Limit trajectory
        ref_P = ref_bb.atom_slice(top_ref_bb.select('212 <= resid and resid <= 221' )) #Limit trajectory to the P loop of PTP1B only
        ref_CYS = ref_bb.atom_slice(top_ref_bb.select('resid == 213' )) #Limit trajectory to CYS215 of PTP1B only
        ref_SBL = ref_bb.atom_slice(top_ref_bb.select('111 <= resid and resid <= 118' )) #Limit trajectory to the S loop of PTP1B only
        ref_a3 = ref_bb.atom_slice(top_ref_bb.select('184 <= resid and resid <= 198' ))
        ref_a3_top = ref_bb.atom_slice(top_ref_bb.select('184 <= resid and resid <= 189' ))
        ref_a4 = ref_bb.atom_slice(top_ref_bb.select('219 <= resid and resid <= 236' ))
        ref_a4_P = ref_bb.atom_slice(top_ref_bb.select('217 <= resid and resid <= 226' ))
        ref_a5 = ref_bb.atom_slice(top_ref_bb.select('243 <= resid and resid <= 250' ))
        ref_a6 = ref_bb.atom_slice(top_ref_bb.select('262 <= resid and resid <= 279' ))
        ref_a6_bot = ref_bb.atom_slice(top_ref_bb.select('273 <= resid and resid <= 279' ))
        if traj_ns.n_residues > 297: #Only compute these distances if the a7 helix is present
            ref_a7 = ref_bb.atom_slice(top_ref_bb.select('285 <= resid and resid <= 293' ))
        ref_L11 = ref_bb.atom_slice(top_ref_bb.select('149 <= resid and resid <= 152' ))
        ref_Q = ref_bb.atom_slice(top_ref_bb.select('257 <= resid and resid <= 261' ))
        ref_beg = ref_bb.atom_slice(top_ref_bb.select('25 <= resid and resid <= 34' ))

    if traj_prot.n_residues == 299 or traj_prot.n_residues == 287: #If protein contains the first residue either with or without the a7 helix residues
        print('True')
        #Caculate RMSD for the WPD loop
        traj_WPD = traj_bb.atom_slice(top_bb.select('176 <= resid and resid <= 184' )) #Limit trajectory to the WPD loop of PTP1B only
        #Caculate RMSD for the connection b/w WPD loop and a3 helix
        traj_WPD_a3 = traj_bb.atom_slice(top_bb.select('184 <= resid and resid <= 187' )) #Limit trajectory to res 185 to 188
        #Caculate RMSD for the P loop
        traj_P = traj_bb.atom_slice(top_bb.select('213 <= resid and resid <= 222' )) #Limit trajectory to the P loop of PTP1B only
        traj_CYS = traj_bb.atom_slice(top_bb.select('resid == 214' )) #Limit trajectory to CYS215 of PTP1B only
        #Caculate RMSD for the Substrate Binding Loop loop
        traj_SBL = traj_bb.atom_slice(top_bb.select('112 <= resid and resid <= 119' )) #Limit trajectory to the S loop of PTP1B only
        #Caculate RMSD for the a3 helix
        traj_a3 = traj_bb.atom_slice(top_bb.select('185 <= resid and resid <= 199' ))
        traj_a3_top = traj_bb.atom_slice(top_bb.select('185 <= resid and resid <= 190'))
        #Caculate RMSD for the a4 helix
        traj_a4 = traj_bb.atom_slice(top_bb.select('220 <= resid and resid <= 237' ))
        traj_a4_P = traj_bb.atom_slice(top_bb.select('218 <= resid and resid <= 227' ))
        #Caculate RMSD for the a5 helix
        traj_a5 = traj_bb.atom_slice(top_bb.select('244 <= resid and resid <= 251' ))
        #Caculate RMSD for the a6 helix
        traj_a6 = traj_bb.atom_slice(top_bb.select('263 <= resid and resid <= 280' ))
        traj_a6_bot = traj_bb.atom_slice(top_bb.select('274 <= resid and resid <= 280' ))

        if traj_ns.n_residues > 297: #Only compute these distances if the a7 helix is present
            #Caculate RMSD for the a7 helix
            traj_a7 = traj_bb.atom_slice(top_bb.select('286 <= resid and resid <= 294' ))
        traj_L11 = traj_bb.atom_slice(top_bb.select('150 <= resid and resid <= 153' ))
        #Caculate RMSD for the Q Loop
        traj_Q = traj_bb.atom_slice(top_bb.select('258 <= resid and resid <= 262' ))
        #Caculate RMSD for the loop b/w a2 and a3
        traj_beg = traj_bb.atom_slice(top_bb.select('26 <= resid and resid <= 35' ))

    else:
        #Caculate RMSD for the WPD loop
        traj_WPD = traj_bb.atom_slice(top_bb.select('175 <= resid and resid <= 183' )) #Limit trajectory to the WPD loop of PTP1B only
        #Caculate RMSD for the connection b/w WPD loop and a3 helix
        traj_WPD_a3 = traj_bb.atom_slice(top_bb.select('183 <= resid and resid <= 186' )) #Limit trajectory to res 185 to 188
        #Caculate RMSD for the P loop
        traj_P = traj_bb.atom_slice(top_bb.select('212 <= resid and resid <= 221' )) #Limit trajectory to the P loop of PTP1B only
        traj_CYS = traj_bb.atom_slice(top_bb.select('resid == 213' )) #Limit trajectory to CYS215 of PTP1B only
        #Caculate RMSD for the Substrate Binding Loop loop
        traj_SBL = traj_bb.atom_slice(top_bb.select('111 <= resid and resid <= 118' )) #Limit trajectory to the S loop of PTP1B only
        #Caculate RMSD for the a3 helix
        traj_a3 = traj_bb.atom_slice(top_bb.select('184 <= resid and resid <= 198' ))
        traj_a3_top = traj_bb.atom_slice(top_bb.select('184 <= resid and resid <= 189'))
        #Caculate RMSD for the a4 helix
        traj_a4 = traj_bb.atom_slice(top_bb.select('219 <= resid and resid <= 236' ))
        traj_a4_P = traj_bb.atom_slice(top_bb.select('217 <= resid and resid <= 226' ))
        #Caculate RMSD for the a5 helix
        traj_a5 = traj_bb.atom_slice(top_bb.select('243 <= resid and resid <= 250' ))
        #Caculate RMSD for the a6 helix
        traj_a6 = traj_bb.atom_slice(top_bb.select('262 <= resid and resid <= 279' ))
        traj_a6_bot = traj_bb.atom_slice(top_bb.select('273 <= resid and resid <= 279' ))
        
        if traj_ns.n_residues > 297: #Only compute these distances if the a7 helix is present
            #Caculate RMSD for the a7 helix
            traj_a7 = traj_bb.atom_slice(top_bb.select('285 <= resid and resid <= 293' ))
        traj_L11 = traj_bb.atom_slice(top_bb.select('150 <= resid and resid <= 153' ))
        traj_Q = traj_bb.atom_slice(top_bb.select('257 <= resid and resid <= 261' ))
        traj_beg = traj_bb.atom_slice(top_bb.select('25 <= resid and resid <= 34' ))

    #Caculate RMSD for the WPD loop
    rmsd_WPD = md.rmsd(traj_WPD, ref_WPD, parallel=True, precentered=False)
    
    #Caculate RMSD for the connection b/w WPD loop and a3 helix
    rmsd_WPD_a3 = md.rmsd(traj_WPD, ref_WPD, parallel=True, precentered=False)
    
    #Caculate RMSD for the P loop
    rmsd_P = md.rmsd(traj_P, ref_P, parallel=True, precentered=False)
    rmsd_CYS = md.rmsd(traj_CYS, ref_CYS, parallel=True, precentered=False)
   
    #Caculate RMSD for the Substrate Binding Loop loop
    rmsd_SBL = md.rmsd(traj_SBL, ref_SBL, parallel=True, precentered=False)
   
    #Caculate RMSD for the a3 helix
    rmsd_a3 = md.rmsd(traj_a3, ref_a3, parallel=True, precentered=False)
    rmsd_a3_top = md.rmsd(traj_a3_top, ref_a3_top, parallel=True, precentered=False)

    #Caculate RMSD for the a4 helix
    rmsd_a4 = md.rmsd(traj_a4, ref_a4, parallel=True, precentered=False)
    rmsd_a4_P = md.rmsd(traj_a4_P, ref_a4_P, parallel=True, precentered=False)
   
    #Caculate RMSD for the a5 helix
    rmsd_a5 = md.rmsd(traj_a5, ref_a5, parallel=True, precentered=False)

    #Caculate RMSD for the a6 helix
    rmsd_a6 = md.rmsd(traj_a6, ref_a6, parallel=True, precentered=False)
    rmsd_a6_bot = md.rmsd(traj_a6_bot, ref_a6_bot, parallel=True, precentered=False)
    
    if traj_ns.n_residues > 297: #Only compute these distances if the a7 helix is present
        #Caculate RMSD for the a7 helix
        rmsd_a7 = md.rmsd(traj_a7, ref_a7, parallel=True, precentered=False)
    
    #Caculate RMSD for the L11
    rmsd_L11 = md.rmsd(traj_L11, ref_L11, parallel=True, precentered=False)
       
    #Caculate RMSD for the a3 helix
    rmsd_Q = md.rmsd(traj_Q, ref_Q, parallel=True, precentered=False)

    #Caculate RMSD for the loop between a2 and a3
    rmsd_beg = md.rmsd(traj_beg, ref_beg, parallel=True, precentered=False)

    #Seperate uncorrelated samples
    t_WPD = uncorr.ind(rmsd_WPD)
    rmsd_WPD_uncorr = uncorr.sort(rmsd_WPD, t_WPD)

    t_WPD_a3 = uncorr.ind(rmsd_WPD_a3)
    rmsd_WPD_a3_uncorr = uncorr.sort(rmsd_WPD, t_WPD_a3)

    t_P = uncorr.ind(rmsd_P)
    rmsd_P_uncorr = uncorr.sort(rmsd_P, t_P)

    t_CYS = uncorr.ind(rmsd_CYS)
    rmsd_CYS_uncorr = uncorr.sort(rmsd_CYS, t_CYS)

    t_SBL = uncorr.ind(rmsd_SBL)
    rmsd_SBL_uncorr = uncorr.sort(rmsd_SBL, t_SBL)

    t_a3 = uncorr.ind(rmsd_a3)
    rmsd_a3_uncorr = uncorr.sort(rmsd_a3, t_a3)

    t_a3 = uncorr.ind(rmsd_a3)
    rmsd_a3_uncorr = uncorr.sort(rmsd_a3, t_a3)

    t_a3_top = uncorr.ind(rmsd_a3_top)
    rmsd_a3_top_uncorr = uncorr.sort(rmsd_a3_top, t_a3_top)

    t_a4 = uncorr.ind(rmsd_a4)
    rmsd_a4_uncorr = uncorr.sort(rmsd_a4, t_a4)

    t_a4_P = uncorr.ind(rmsd_a4_P)
    rmsd_a4_P_uncorr = uncorr.sort(rmsd_a4_P, t_a4_P)

    t_a5 = uncorr.ind(rmsd_a5)
    rmsd_a5_uncorr = uncorr.sort(rmsd_a5, t_a5)

    t_a6 = uncorr.ind(rmsd_a6)
    rmsd_a6_uncorr = uncorr.sort(rmsd_a6, t_a6)

    t_a6_bot = uncorr.ind(rmsd_a6_bot)
    rmsd_a6_bot_uncorr = uncorr.sort(rmsd_a6_bot, t_a6_bot)

    if traj_ns.n_residues > 297: #Only compute these distances if the a7 helix is present
        t_a7 = uncorr.ind(rmsd_a7)
        rmsd_a7_uncorr = uncorr.sort(rmsd_a7, t_a7)

    t_L11 = uncorr.ind(rmsd_L11)
    rmsd_L11_uncorr = uncorr.sort(rmsd_L11, t_L11)

    t_Q = uncorr.ind(rmsd_Q)
    rmsd_Q_uncorr = uncorr.sort(rmsd_Q, t_Q)

    t_beg = uncorr.ind(rmsd_beg)
    rmsd_beg_uncorr = uncorr.sort(rmsd_beg, t_beg)


    np.savetxt('rmsd_WPD_ref_' + str(ref_name) + '.txt', rmsd_WPD_uncorr)
    np.savetxt('rmsd_WPD_a3_ref_' + str(ref_name) + '.txt', rmsd_WPD_a3_uncorr)
    np.savetxt('rmsd_P_ref_' + str(ref_name) + '.txt', rmsd_P_uncorr)
    np.savetxt('rmsd_CYS_ref_' + str(ref_name) + '.txt', rmsd_CYS_uncorr)
    np.savetxt('rmsd_SBL_ref_' + str(ref_name) + '.txt', rmsd_SBL_uncorr)
    np.savetxt('rmsd_a3_ref_' + str(ref_name) + '.txt', rmsd_a3_uncorr)
    np.savetxt('rmsd_a3_top_ref_' + str(ref_name) + '.txt', rmsd_a3_top_uncorr)
    np.savetxt('rmsd_a4_ref_' + str(ref_name) + '.txt', rmsd_a4_uncorr)
    np.savetxt('rmsd_a4_P_ref_' + str(ref_name) + '.txt', rmsd_a4_P_uncorr)
    np.savetxt('rmsd_a5_ref_' + str(ref_name) + '.txt', rmsd_a5_uncorr)
    np.savetxt('rmsd_a6_ref_' + str(ref_name) + '.txt', rmsd_a6_uncorr)
    np.savetxt('rmsd_a6_bot_ref_' + str(ref_name) + '.txt', rmsd_a6_bot_uncorr)
    if traj_ns.n_residues > 297: #Only compute these distances if the a7 helix is present
        np.savetxt('rmsd_a7_ref_' + str(ref_name) + '.txt', rmsd_a7_uncorr)
    np.savetxt('rmsd_L11_ref_' + str(ref_name) + '.txt', rmsd_L11_uncorr)
    np.savetxt('rmsd_Q_ref_' + str(ref_name) + '.txt', rmsd_Q_uncorr)
    np.savetxt('rmsd_beg_ref_' + str(ref_name) + '.txt', rmsd_beg_uncorr)
    
    print(len(t_full))
    print(len(t_WPD))
    print('RMSD and RMSF Analysis Completed')
else:
    print('RMSD and RMSF Analysis Skipped')

#Only do DSSP if requested when running
if dssp_check == True:
    #Compute Phi and Psi angles for all residues in the a7 helix in all frames
    phi_ind, phi_angle = md.compute_phi(traj_a7, periodic = True, opt = True)
    psi_ind, psi_angle = md.compute_psi(traj_a7, periodic = True, opt = True)
    time, angles = np.shape(phi_angle)

    #Limit to uncorrelated data
    phi_uncorr = np.zeros((len(t_full), angles))
    psi_uncorr = np.zeros((len(t_full), angles))
    for i in range(angles):
        phi_uncorr[:,i] = uncorr.sort(phi_angle[:,i], t_full)
        psi_uncorr[:,i] = uncorr.sort(psi_angle[:,i], t_full)

    #Compute Secondary Structure
    dssp_list = md.compute_dssp(traj_a7, simplified=False) #Compute DSSP for all residues in the a7 helix for all trajectory frames
    file_dssp = open('DSSP_'+ File_base + '.txt','w') #Create output file for DSSP and write over if file is present

    #Limit to Uncorrelated data
    frame_max,residue=dssp_list.shape
    print(residue)
    dssp_uncorr = np.full((len(t_full) - 1, residue), None)
    for i in range(residue):
        dssp_res = dssp_list[:,i] #Seperate all time values for a single residue
        dssp_res_mod = []
        for j in dssp_res:
            if j == ' ':
                dssp_res_mod.append('L')
            else:
                dssp_res_mod.append(j)
        dssp_uncorr[:,i] = uncorr.char(dssp_res_mod, t_full)
    
    #Output DSSP to file
    for i in range(len(t_full) - 1):
        for j in range(residue):
            file_dssp.write(dssp_uncorr[i,j] + ' ')
        file_dssp.write('\n')
    file_dssp.close() #close file
    
    #Save psi and phi angles to files and overwrite if file is present
    np.savetxt('phi_' + File_base + '.txt', phi_angle)
    np.savetxt('psi_' + File_base + '.txt', psi_angle)

    #Print to screen to notify that DSSP file has been written properly 
    print('DSSP File Written')
    
#Skip DSSP Analysis if not requested
else:
    print('DSSP Skipped')

#H-bond determination
if hbond_check == True:
    #Determine list of H-bonds present in the trajectory for over 60% of the frames
    hbonds = md.baker_hubbard(traj_ns, freq=0.6, exclude_water=True, periodic=False)
    label = lambda hbond : '%s -- %s' % (traj_ns.topology.atom(hbond[0]), traj_ns.topology.atom(hbond[2])) #Extract labels for h-bonds
    np.savetxt('Hbonds_atom_' + File_base + '.txt', hbonds) #Save all atom indicies of h-bonds to file for further analysis
    
    #Write all h-bonds present for >60% of trajectory to file
    file_object = open('Hbonds_'+ File_base +'.txt', 'w') 
    for hbond in hbonds:
        file_object.write(label(hbond)) #Maintain same order as atom indicies
        file_object.write('\n')
    file_object.close() #close file

    #Determine the exact percentage of time that each h-bond present for >60% of the trajectory is formed
    per = [] #Declare empty array for percentage of time h-bond is formed
    da_distances = md.compute_distances(traj_ns, hbonds[:,[1,2]], periodic=False) #Compute distance between h-bond donor and acceptor
    da_angles = md.compute_angles(traj_ns, hbonds[:,:], periodic=False) #Compute angle between h-bond donor and acceptor
    [num_t, num_h] = np.shape(da_distances) #save values for number of frames(num_t) and number of bonds(num_b) to caculate
    for j in range(num_h): #Loop through all h-bonds
        count = 0 #Initialize count
        for i in range(num_t): #Loop through all frames
            if da_distances[i,j] <= 0.25 and da_angles[i,j] >= 2.094: #If distance between donor and acceptor is less than 2.5A and the angle is greater than 120 degrees or ~ 2.094 radians
                count +=1
        per.append(100*count/num_t) #Percentage of time each h-bond is present in trajectory
    np.savetxt('Hbonds_per_' + File_base + '.txt',per)

    #If ligand analysis is requested determine h-bonds between ligand and PTP1B residues
    if lig_check == True:
        file_lig = open('Hbonds_lig_' +  File_base + '.txt', 'w') #Create file for list of h-bonds determined to be present more than 10% of the trajectory
        if lig == 'both':            
            ligand = traj_ns.topology.select('resname AD') #Select the ligand by name (based on topology) from the trajectory
            ligand2 = traj_ns.topology.select('resname BBR') #Select the ligand by name (based on topology) from the trajectory
        else:
            ligand = traj_ns.topology.select('resname ' + str(lig)) #Select the ligand by name (based on topology) from the trajectory
        protein = traj_ns.topology.select('protein') #Select protein from the trajectory
        hbonds = md.baker_hubbard(traj_ns, freq = 0.1, exclude_water=True, periodic=True) #Find all h-bonds present >10% of the trajectory
        label = lambda hbond : '%s -- %s' % (traj.topology.atom(hbond[0]), traj.topology.atom(hbond[2])) #Seperate the names of the h-bonds based on topology nomenclature
        if lig == 'both':
            file_lig.write('AD:\n')
            for hbond in hbonds: #search all h-bonds
                if (hbond[0] in ligand) and (hbond[2] in protein) or (hbond[2] in ligand) and (hbond[0] in protein): #seperate only those h-bonds which form between the ligand and the protein
                    file_lig.write(str(label(hbond)) + '\n') #Write h-bond names to file
                    file_lig.write(str(hbond[0]) + ' ' + str(hbond[1]) + ' ' + str(hbond[2]) + '\n') #Write atoms involved in h-bond to file
            file_lig.write('BBR:\n')
            for hbond in hbonds: #search all h-bonds
                if (hbond[0] in ligand2) and (hbond[2] in protein) or (hbond[2] in ligand2) and (hbond[0] in protein): #seperate only those h-bonds which form between the ligand and the protein
                    file_lig.write(str(label(hbond)) + '\n') #Write h-bond names to file
                    file_lig.write(str(hbond[0]) + ' ' + str(hbond[1]) + ' ' + str(hbond[2]) + '\n') #Write atoms involved in h-bond to file
            file_lig.close() #close file
    print('Hbond Analysis Written')

#Skip Hbond analysis if desired
else:
    print('Hbond Analysis Skipped')

#Compute Interactions between the a3, a6, and a7 helices
if check_hel == True:
    #Set Residue Pairs
    if traj_prot.n_residues == 299 or traj_prot.n_residues == 287: #If protein contains the first residue either with or without the a7 helix residues
        group_WPD = np.linspace(176, 184, num = 9) #residues in the WPD loop
        group_3 = np.linspace(185, 199, num = 15) #residues in the a3 helix
        group_4 = np.linspace(220, 237, num = 18) #residues in the a4 helix
        group_5 = np.linspace(244, 251, num = 8) #residues in the a5 helix
        group_6 = np.linspace(263, 280, num = 18) #residues in the a6 helix
        group_bend = np.linspace(281, 285, num = 5) #residues in the bend b/w the a6 and a7 helices
        group_7 = np.linspace(286, 297, num = 12) #residues in the a7 helix
        group_L11 = np.linspace(149, 152, num = 4) #residues in the L11 loop
        pair_other = [[150, 190], [263, 184], [149, 177], [80, 198], [116, 181], [116, 216], [191, 224], [181, 219], [186, 268], [177, 189], [151, 176]] #Residue pairs for distances not included above

        #Set the index for the ligand
        if traj_ns.n_residues == 300: #w/ a7 helix
            if lig == 'both':
                group_l = [299]
                group_l2 = [300]
            else:
                group_l = [299]
        if traj_ns.n_residues == 288: #w/o a7 helix
            if lig == 'both':
                group_l = [287]
                group_l2 = [288]
            else:
                group_l = [287]
    else: #Missing first residue leads to shifted indicies
        group_WPD = np.linspace(175, 183, num = 9) #residues in the WPD loop
        group_3 = np.linspace(184, 198, num = 15) #residues in the a3 helix
        group_4 = np.linspace(219, 236, num = 18) #residues in the a4 helix
        group_5 = np.linspace(243, 250, num = 8) #residues in the a5 helix
        group_6 = np.linspace(262, 279, num = 18) #residues in the a6 helix
        group_bend = np.linspace(280, 284, num = 5) #residues in the bend b/w the a6 and a7 helices
        group_7 = np.linspace(285, 296, num = 12) #residues in the a7 helix
        group_L11 = np.linspace(148, 151, num = 4) #residues in the L11 loop
        pair_other = [[149, 189], [262, 183], [148, 176], [79, 197], [115, 180], [115, 215], [190, 223], [180, 218], [185, 267], [176, 188], [150, 175]] #Residue pairs for distances not included above

        #Set the index for the ligand
        if traj_ns.n_residues == 299: #w/ a7 helix
            if lig == 'both':
                group_l = [298]
                group_l2 = [299]
            else:
                group_l = [298]
        if traj_ns.n_residues == 287: #w/o a7 helix
            if lig == 'both':
                group_l = [286]
                group_l2 = [287]
            else:
                group_l = [286]
    
    #Set pairs to compute distances
    pair_a3_a6 = list(product(group_3, group_6))
    pair_a3_WPD = list(product(group_3, group_WPD))
    pair_a3_bend = list(product(group_3, group_bend))
    if traj_ns.n_residues > 297: #Only include these pairs if the a7 is present in the trajectory
        pair_a7_a3 = list(product(group_3, group_7))
        pair_a7_a6 = list(product(group_6, group_7))
        pair_a7_L11 = list(product(group_L11, group_7))

    #Set indicies for sections of helices as described from their standard orientation for visualization
    a3_a7_pt1_ind = [5, 6, 7, 8, 9, 10, 11, 17, 18, 19, 20, 21, 22, 23, 29, 30, 31, 32, 33, 34, 35, 42, 43, 44, 45, 46, 47, 48, 54, 55, 56, 57, 58, 59, 60, 67, 68, 69, 70, 71, 72, 73, 79, 80, 81, 82, 83, 84, 85, 91, 92, 93, 94, 
            95, 96, 97] #Upper region of both helices
    a3_a7_pt2_ind = [98, 99, 100, 101, 102, 110, 111, 112, 113, 114, 122, 123, 124, 125, 126, 134, 135, 136, 137, 138, 146, 147, 148, 149, 150, 158, 159, 160, 161, 162, 170, 171, 172, 173, 174] #Lower region of both helices
    a6_a7_pt1_ind = [5, 6, 7, 8, 9, 10, 11, 17, 18, 19, 20, 21, 22, 23, 29, 30, 31, 32, 33, 34, 35, 42, 43, 44, 45, 46, 47, 48, 54, 55, 56, 57, 58, 59, 60, 67, 68, 69, 70, 71, 72, 73, 79, 80, 81, 82, 83, 84, 85, 91, 92, 93, 94, 
            95, 96, 97, 103, 104, 105, 106, 107, 108, 109, 115, 116, 117, 118, 119, 120, 121] #Upper region of both helices
    a6_a7_pt2_ind = [122, 123, 124, 125, 126, 134, 135, 136, 137, 138, 146, 147, 148, 149, 150, 158, 159, 160, 161, 162, 170, 171, 172, 173, 174, 182, 183, 184, 185, 186, 194, 195, 196, 197, 198, 206, 207, 208, 209, 
            210] #Lower region of both Helices
    a6_a7_pt3_ind = [127, 128, 129, 130, 131, 132, 133, 139, 140, 141, 142, 143, 144, 145, 151, 152, 153, 154, 155, 156, 157, 163, 164, 165, 166, 167, 168, 169, 175, 176, 177, 178, 179, 180, 181, 187, 188, 189, 190, 191, 192, 
            193, 199, 200, 201, 202, 203, 204, 205, 211, 212, 213, 214, 215, 216, 217] #Upper region of a7 with lower region of the a6 helix

    #Compute distances for all pairs for closest atoms to compute contacts
    [dist_a3_a6_all, pairs] = md.compute_contacts(traj_ns, contacts=pair_a3_a6, scheme='closest', ignore_nonprotein = False, periodic=True, soft_min = False)
    [dist_a3_bend_all, pairs] = md.compute_contacts(traj_ns, contacts=pair_a3_bend, scheme='closest', ignore_nonprotein = False, periodic=True, soft_min = False)
    [dist_a3_WPD_all, pairs] = md.compute_contacts(traj_ns, contacts=pair_a3_WPD, scheme='closest', ignore_nonprotein = False, periodic=True, soft_min = False)
    if traj_ns.n_residues > 297: #Only compute these distances if the a7 helix is present
        [dist_a7_a3_all, pairs] = md.compute_contacts(traj_ns, contacts=pair_a7_a3, scheme='closest', ignore_nonprotein = False, periodic=True, soft_min = False)
        [dist_a7_a6_all, pairs] = md.compute_contacts(traj_ns, contacts=pair_a7_a6, scheme='closest', ignore_nonprotein = False, periodic=True, soft_min = False)
        [dist_a7_L11_all, pairs] = md.compute_contacts(traj_ns, contacts=pair_a7_L11, scheme='closest', ignore_nonprotein = False, periodic=True, soft_min = False)
    #Compute distances for all pairs for alpha carbonds to compute distances between residues
    [dist_ca_a3_WPD_all, pairs] = md.compute_contacts(traj_ns, contacts=pair_a3_WPD, scheme='ca', ignore_nonprotein = False, periodic=True, soft_min = False)
    [dist_ca_a3_bend_all, pairs] = md.compute_contacts(traj_ns, contacts=pair_a3_bend, scheme='ca', ignore_nonprotein = False, periodic=True, soft_min = False)
    [dist_ca_other_all, pairs] = md.compute_contacts(traj_ns, contacts=pair_other, scheme='ca', ignore_nonprotein = False, periodic=True, soft_min = False)
    if traj_ns.n_residues > 297: #Only compute these distances if the a7 helix is present
        [dist_ca_a7_a3_all, pairs] = md.compute_contacts(traj_ns, contacts=pair_a7_a3, scheme='ca', ignore_nonprotein = False, periodic=True, soft_min = False)
        [dist_ca_a7_a6_all, pairs] = md.compute_contacts(traj_ns, contacts=pair_a7_a6, scheme='ca', ignore_nonprotein = False, periodic=True, soft_min = False)
        [dist_ca_a7_L11_all, pairs] = md.compute_contacts(traj_ns, contacts=pair_a7_L11, scheme='ca', ignore_nonprotein = False, periodic=True, soft_min = False)

    #Compute number of time points and number of diatance pairs for following loops
    time, num_pairs_a3_a6 = np.shape(dist_a3_a6_all)
    time, num_pairs_a3_bend = np.shape(dist_a3_bend_all)
    time, num_pairs_a3_WPD = np.shape(dist_a3_WPD_all)
    time, num_pairs_ca_other = np.shape(dist_ca_other_all)
    if traj_ns.n_residues > 297:
        time, num_pairs_a7_a3 = np.shape(dist_a7_a3_all)
        time, num_pairs_a7_a6 = np.shape(dist_a7_a6_all)
        time, num_pairs_a7_L11 = np.shape(dist_a7_L11_all)
    
    #Limit to uncorrelated samples
    test_dist = dist_a3_a6_all[:,0] #Seperate a test case for determining correlated intervals
    t_dist = uncorr.ind(test_dist) #Determine indices of uncorrelated samples
    time_uncorr = len(t_dist)

    #Set open arrays for all samples
    dist_a3_a6 = np.zeros((time_uncorr, num_pairs_a3_a6))
    dist_a3_bend = np.zeros((time_uncorr, num_pairs_a3_bend))
    dist_ca_a3_bend = np.zeros((time_uncorr, num_pairs_a3_bend))
    dist_a3_WPD = np.zeros((time_uncorr, num_pairs_a3_WPD))
    dist_ca_a3_WPD = np.zeros((time_uncorr, num_pairs_a3_WPD))
    dist_ca_other = np.zeros((time_uncorr, num_pairs_ca_other))
    if traj_ns.n_residues > 297:
        dist_a7_a3 = np.zeros((time_uncorr, num_pairs_a7_a3))
        dist_ca_a7_a3 = np.zeros((time_uncorr, num_pairs_a7_a3))
        dist_a7_a6 = np.zeros((time_uncorr, num_pairs_a7_a6))
        dist_ca_a7_a6 = np.zeros((time_uncorr, num_pairs_a7_a6))
        dist_a7_L11 = np.zeros((time_uncorr, num_pairs_a7_L11))
        dist_ca_a7_L11 = np.zeros((time_uncorr, num_pairs_a7_L11))
    
    #Set new arrays with uncorrelated samples
    for i in range(num_pairs_a3_a6):
        dist = dist_a3_a6_all[:,i]
        dist_a3_a6[:,i] = uncorr.sort(dist, t_dist)
    for i in range(num_pairs_a3_bend):
        dist = dist_a3_bend_all[:,i]
        dist_a3_bend[:,i] = uncorr.sort(dist, t_dist)
        dist = dist_ca_a3_bend_all[:,i]
        dist_ca_a3_bend[:,i] = uncorr.sort(dist, t_dist)
    for i in range(num_pairs_a3_WPD):
        dist = dist_a3_WPD_all[:,i]
        dist_a3_WPD[:,i] = uncorr.sort(dist, t_dist)
        dist = dist_ca_a3_WPD_all[:,i]
        dist_ca_a3_WPD[:,i] = uncorr.sort(dist, t_dist)
    for i in range(num_pairs_ca_other):
        dist = dist_ca_other_all[:,i]
        dist_ca_other[:,i] = uncorr.sort(dist, t_dist)
    if traj_ns.n_residues > 297: #Only open these files if the a7 helix is present
        for i in range(num_pairs_a7_a3):
            dist = dist_a7_a3_all[:,i]
            dist_a7_a3[:,i] = uncorr.sort(dist, t_dist)
            dist = dist_ca_a7_a3_all[:,i]
            dist_ca_a7_a3[:,i] = uncorr.sort(dist, t_dist)
        for i in range(num_pairs_a7_a6):
            dist = dist_a7_a6_all[:,i]
            dist_a7_a6[:,i] = uncorr.sort(dist, t_dist)
            dist = dist_ca_a7_a6_all[:,i]
            dist_ca_a7_a6[:,i] = uncorr.sort(dist, t_dist)
        for i in range(num_pairs_a7_L11):
            dist = dist_a7_L11_all[:,i]
            dist_a7_L11[:,i] = uncorr.sort(dist, t_dist)
            dist = dist_ca_a7_L11_all[:,i]
            dist_ca_a7_L11[:,i] = uncorr.sort(dist, t_dist)

    #Open files for the number of interactions with each protein region at each point in time
    file_200_282 = open('200_282_dist.txt', 'w') #a3 to a6
    file_179_191 = open('179_191_dist.txt', 'w') #WPD to a3
    file_185_191 = open('185_191_dist.txt', 'w') #WPD to a3
    file_151_191 = open('151_191_dist.txt', 'w') #L11 to a3
    file_264_185 = open('264_185_dist.txt', 'w') #WPD to a6
    file_178_150 = open('178_150_dist.txt', 'w') #WPD to L11
    file_81_199 = open('81_199_dist.txt', 'w')
    file_117_182 = open('117_182_dist.txt', 'w') #WPD to SBL
    file_117_217 = open('117_217_dist.txt', 'w') #P to SBL
    file_192_225 = open('192_225_dist.txt', 'w') #a3 top to a5
    file_182_220 = open('182_220_dist.txt', 'w') #WPD tp P-loop
    file_187_269 = open('187_269_dist.txt', 'w') #a3 top to a6 top
    file_178_190 = open('178_190_dist.txt', 'w') #Start of WPD loop to side of a3
    file_152_177 = open('152_177_dist.txt', 'w') #Start of WPD loop to L11 loop
    file_mean = open('inter_mean.txt', 'w') #File for the mean number of interactions between any two given residues
    if traj_ns.n_residues > 297: #Only open these files if the a7 helix is present
        file_200_287 = open('200_287_dist.txt', 'w') #a3 to a7 (bottom)
        file_276_292 = open('276_292_dist.txt', 'w') #a6 to a7 (middle)
        file_189_295 = open('189_295_dist.txt', 'w') #a3 to a7 (top)
        file_280_287 = open('280_287_dist.txt', 'w') #a6 to a7 (bottom)
        file_152_297 = open('152_297_dist.txt', 'w') #L11 to a7

    #Declare arrays to keep track of all contacts
    a3_a6_all, a3_a7_all, a6_a7_all, a7_L11_all = [],[],[],[] #Arrays for the total contacts made b/w helices
    a3_a7_pt1, a3_a7_pt2, a6_a7_pt1, a6_a7_pt2, a6_a7_pt3 = [],[],[],[],[] #Arrays for the total contacts made b/w sections of helices
    if traj_ns.n_residues > 297: #Only create these arrays if the a7 helix is present
        inter_a3_a7 = np.zeros([num_pairs_a7_a3, time_uncorr]) #Array for the running total of interactions over time between each residue pair for tha a3 and a7 helices
        inter_a6_a7 = np.zeros([num_pairs_a7_a6, time_uncorr]) #Array for the running total of interactions over time between each residue pair for tha a3 and a7 helices
    
    #Loop through all time points
    for i in range(time_uncorr-1):
        #Set all interaction counters to zero for each new time point
        check_a3_a6 = 0
        check_a7_a3 = 0
        check_a7_a3_pt1 = 0
        check_a7_a3_pt2 = 0
        check_a7_a6 = 0
        check_a7_a6_pt1 = 0
        check_a7_a6_pt2 = 0
        check_a7_a6_pt3 = 0
        check_a7_L11 = 0

        #At each time point save the distance between each residue pair of interest
        file_179_191.write(str(dist_ca_a3_WPD[i][47]) + '\n')
        file_185_191.write(str(dist_ca_a3_WPD[i][53]) + '\n')
        file_200_282.write(str(dist_ca_a3_bend[i][70]) + '\n')
        file_151_191.write(str(dist_ca_other[i][0]) + '\n')
        file_264_185.write(str(dist_ca_other[i][1]) + '\n')
        file_178_150.write(str(dist_ca_other[i][2]) + '\n')
        file_81_199.write(str(dist_ca_other[i][3]) + '\n')
        file_117_182.write(str(dist_ca_other[i][4]) + '\n')
        file_117_217.write(str(dist_ca_other[i][5]) + '\n')
        file_192_225.write(str(dist_ca_other[i][6]) + '\n')
        file_182_220.write(str(dist_ca_other[i][7]) + '\n')
        file_187_269.write(str(dist_ca_other[i][8]) + '\n')
        file_178_190.write(str(dist_ca_other[i][9]) + '\n')
        file_152_177.write(str(dist_ca_other[i][10]) + '\n')

        if traj_ns.n_residues > 297: #Only save these distances if the a7 helix is present
            file_200_287.write(str(dist_ca_a7_a3[i][168]) + '\n')
            file_189_295.write(str(dist_ca_a7_a3[i][44]) + '\n')
            file_276_292.write(str(dist_ca_a7_a6[i][149]) + '\n')
            file_280_287.write(str(dist_ca_a7_a6[i][192]) + '\n')
            file_152_297.write(str(dist_ca_a7_L11[i][34]) + '\n')
        
        #Loop through all residue pairs in the a3 and a6 helices
        for j in range(num_pairs_a3_a6): #Determine # of contacts b/w the a3 and a6 helices
            if dist_a3_a6[i][j] <= 0.5: #Count the contact if the residue distance is less than 0.5 nm
                check_a3_a6 += 1
        a3_a6_all.append(check_a3_a6) #Save the total number of residue contacts b/w the a3 and a6 helices to array
        if traj_ns.n_residues > 297: #Only calculate these interactions if the a7 helix is present
            for k in range(num_pairs_a7_a3): #Determine # of contacts b/w the a3 and a7 helices
                if dist_a7_a3[i][k] <= 0.5:
                    check_a7_a3 += 1
                    inter_a3_a7[k][i] = 1 #Record presence of contact for given residue pair(k) and time(i)
                    if k in a3_a7_pt1_ind: #If this residue pair index is in part 1 of the a3 and a7 helices
                        check_a7_a3_pt1 += 1
                    if k in a3_a7_pt2_ind: #If this residue pair index is in part 2 of the a3 and a7 helices
                        check_a7_a3_pt2 += 1
            #Save to array/file for each time point
            a3_a7_all.append(check_a7_a3)
            a3_a7_pt1.append(check_a7_a3_pt1)
            a3_a7_pt2.append(check_a7_a3_pt2)
            for l in range(num_pairs_a7_a6): #Determine # of contacts b/w the a6 and a7 helices
                if dist_a7_a6[i][l] <= 0.5:
                    check_a7_a6 += 1
                    inter_a6_a7[l][i] = 1 #Record presence of contact for given residue pair(l) and time(i)
                    if l in a6_a7_pt1_ind: #If this residue pair index is in part 1 of the a6 and a7 helices
                        check_a7_a6_pt1 += 1
                    if l in a6_a7_pt2_ind: #If this residue pair index is in part 2 of the a6 and a7 helices
                        check_a7_a6_pt2 += 1
                    if l in a6_a7_pt3_ind: #If this residue pair index is in part 3 of the a6 and a7 helices
                        check_a7_a6_pt3 += 1
            #Save to array/file for each time point
            a6_a7_all.append(check_a7_a6)
            a6_a7_pt1.append(check_a7_a6_pt1)
            a6_a7_pt2.append(check_a7_a6_pt2)
            a6_a7_pt3.append(check_a7_a6_pt3)
            for m in range(num_pairs_a7_L11): #Determine # of contacts b/w the L11 loop and a7 helix
                if dist_a7_L11[i][m] <= 0.5:
                    check_a7_L11 += 1
            a7_L11_all.append(check_a7_L11)
    
    np.savetxt('a3_a6_inter.txt', np.asarray(a3_a6_all))
    if traj_ns.n_residues > 297:
        #Save interactions for 2 part helix interactions
        np.savetxt('a3_a7_pt1_tot_inter.txt', np.asarray(a3_a7_pt1))
        np.savetxt('a6_a7_pt1_tot_inter.txt', np.asarray(a6_a7_pt1))
        np.savetxt('a3_a7_pt2_tot_inter.txt', np.asarray(a3_a7_pt2))
        np.savetxt('a6_a7_pt2_tot_inter.txt', np.asarray(a6_a7_pt2))
        np.savetxt('a6_a7_pt3_tot_inter.txt', np.asarray(a6_a7_pt3))
        np.savetxt('a7_a3_inter.txt', np.asarray(a3_a7_all))
        np.savetxt('a7_a6_inter.txt', np.asarray(a6_a7_all))
        np.savetxt('a7_L11_inter.txt', np.asarray(a7_L11_all))

        #Save all individual interactions to file
        file_a3_a7_inters = open('a3_a7_inter_all.txt', 'w')
        for i in range(num_pairs_a7_a3):
            file_a3_a7_inters.write(str(100*sum(inter_a3_a7[i,:])/time_uncorr) + '\n')
        file_a6_a7_inters = open('a6_a7_inter_all.txt', 'w')
        for j in range(num_pairs_a7_a6):
            file_a6_a7_inters.write(str(100*sum(inter_a6_a7[j,:])/time_uncorr) + '\n')
        #Save mean # of interactions over all time points
        file_mean.write('a6-a7 inters mean: ' + str(sum(a6_a7_all)/time_uncorr) + '\n')
        file_mean.write('a3-a7 inters mean: ' + str(sum(a3_a7_all)/time_uncorr))
    file_mean.write('a3-a6 inters mean: ' + str(sum(a3_a6_all)/time_uncorr) + '\n')
    
    #Close all files
    file_200_282.close()
    file_179_191.close()
    file_185_191.close()
    file_151_191.close()
    file_264_185.close()
    file_178_150.close()
    file_mean.close()
    if traj_ns.n_residues > 297: #Only close files that were opened
        file_200_287.close()
        file_276_292.close()
        file_189_295.close()
        file_280_287.close()
        file_152_297.close()

    print('Helix Interaction Analysis Completed')
else:
    print('Helix Interaction Analysis Skipped')

#Ligand Interaction analysis
if lig_check == True: 
    #Compute Ligand location
    pair_a3 = list(product(group_l, group_3))
    pair_a4 = list(product(group_l, group_4))
    pair_a5 = list(product(group_l, group_5))
    pair_a6 = list(product(group_l, group_6))
    pair_bend = list(product(group_l, group_bend))
    if lig == 'both':
        pair2_a3 = list(product(group_l2, group_3))
        pair2_a4 = list(product(group_l2, group_4))
        pair2_a5 = list(product(group_l2, group_5))
        pair2_a6 = list(product(group_l2, group_6))
        pair2_bend = list(product(group_l2, group_bend))

    if traj_ns.n_residues == 300 or traj_ns.n_residues == 299:
        pair_a7 = list(product(group_l, group_7))
        if lig == 'both':
            pair2_a7 = list(product(group_l2, group_7))

        #set up array for total number of contacts with each residue
        tot_pairs = len(pair_a3) + len(pair_a6) + len(pair_bend) + len(pair_a7)
        if lig == 'both':
            tot_pairs2 = len(pair2_a3) + len(pair2_a6) + len(pair2_bend) + len(pair2_a7)

    else:
        #set up array for total number of contacts with each residue
        tot_pairs = len(pair_a3) + len(pair_a6) + len(pair_bend)
        if lig == 'both':
            tot_pairs2 = len(pair2_a3) + len(pair2_a6) + len(pair2_bend)

    #Array for the total number of contacts during the trajectory for each residue in the helices with the ligand
    lig_tot_cont = np.zeros(tot_pairs)
    if lig == 'both':
        lig2_tot_cont = np.zeros(tot_pairs2)

    #compute distances
    [dist_a3_all, pairs] = md.compute_contacts(traj_ns, contacts=pair_a3, scheme='closest', ignore_nonprotein = False, periodic=True, soft_min = False)
    [dist_a4_all, pairs] = md.compute_contacts(traj_ns, contacts=pair_a4, scheme='closest', ignore_nonprotein = False, periodic=True, soft_min = False)
    [dist_a5_all, pairs] = md.compute_contacts(traj_ns, contacts=pair_a5, scheme='closest', ignore_nonprotein = False, periodic=True, soft_min = False)
    [dist_a6_all, pairs] = md.compute_contacts(traj_ns, contacts=pair_a6, scheme='closest', ignore_nonprotein = False, periodic=True, soft_min = False)
    [dist_bend_all, pairs] = md.compute_contacts(traj_ns, contacts=pair_bend, scheme='closest', ignore_nonprotein = False, periodic=True, soft_min = False)
    if traj_ns.n_residues == 300 or traj_ns.n_residues == 299:
        [dist_a7_all, pairs] = md.compute_contacts(traj_ns, contacts=pair_a7, scheme='closest', ignore_nonprotein = False, periodic=True, soft_min = False)
    if lig == 'both':
        [dist2_a3_all, pairs] = md.compute_contacts(traj_ns, contacts=pair2_a3, scheme='closest', ignore_nonprotein = False, periodic=True, soft_min = False)
        [dist2_a4_all, pairs] = md.compute_contacts(traj_ns, contacts=pair2_a4, scheme='closest', ignore_nonprotein = False, periodic=True, soft_min = False)
        [dist2_a5_all, pairs] = md.compute_contacts(traj_ns, contacts=pair2_a5, scheme='closest', ignore_nonprotein = False, periodic=True, soft_min = False)
        [dist2_a6_all, pairs] = md.compute_contacts(traj_ns, contacts=pair2_a6, scheme='closest', ignore_nonprotein = False, periodic=True, soft_min = False)
        [dist2_bend_all, pairs] = md.compute_contacts(traj_ns, contacts=pair2_bend, scheme='closest', ignore_nonprotein = False, periodic=True, soft_min = False)
        if traj_ns.n_residues == 300 or traj_ns.n_residues == 299:
            [dist2_a7_all, pairs] = md.compute_contacts(traj_ns, contacts=pair2_a7, scheme='closest', ignore_nonprotein = False, periodic=True, soft_min = False)

    #Compute time points in which the distance to any residue is less than
    time, num_pairs_a3 = np.shape(dist_a3_all)
    time, num_pairs_a4 = np.shape(dist_a4_all)
    time, num_pairs_a5 = np.shape(dist_a5_all)
    time, num_pairs_a6 = np.shape(dist_a6_all)
    time, num_pairs_bend = np.shape(dist_bend_all)
    if traj_ns.n_residues == 300 or traj_ns.n_residues == 299:
        time, num_pairs_a7 = np.shape(dist_a7_all)

    #Limit to uncorrelated samples
    t_dist = t_full #Determine indices of uncorrelated samples
    time_uncorr = len(t_dist)
    print(time_uncorr)

    #Set open arrays for all samples
    dist_a3 = np.zeros((time_uncorr, num_pairs_a3))
    dist_bend = np.zeros((time_uncorr, num_pairs_bend))
    dist_a4 = np.zeros((time_uncorr, num_pairs_a4))
    dist_a5 = np.zeros((time_uncorr, num_pairs_a5))
    dist_a6 = np.zeros((time_uncorr, num_pairs_a6))
    if traj_ns.n_residues > 297:
        dist_a7 = np.zeros((time_uncorr, num_pairs_a7))
    if lig == 'both':
        dist2_a3 = np.zeros((time_uncorr, num_pairs_a3))
        dist2_bend = np.zeros((time_uncorr, num_pairs_bend))
        dist2_a4 = np.zeros((time_uncorr, num_pairs_a4))
        dist2_a5 = np.zeros((time_uncorr, num_pairs_a5))
        dist2_a6 = np.zeros((time_uncorr, num_pairs_a6))
        if traj_ns.n_residues > 297:
            dist2_a7 = np.zeros((time_uncorr, num_pairs_a7))

    #Set new arrays with uncorrelated samples
    for i in range(num_pairs_a3):
        dist = dist_a3_all[:,i]
        dist_a3[:,i] = uncorr.sort(dist, t_dist)
        if lig == 'both':
            dist2 = dist2_a3_all[:,i]
            dist2_a3[:,i] = uncorr.sort(dist2, t_dist)

    for i in range(num_pairs_a4):
        dist = dist_a4_all[:,i]
        dist_a4[:,i] = uncorr.sort(dist, t_dist)
        if lig == 'both':
            dist2 = dist2_a4_all[:,i]
            dist2_a4[:,i] = uncorr.sort(dist2, t_dist)

    for i in range(num_pairs_a5):
        dist = dist_a5_all[:,i]
        dist_a5[:,i] = uncorr.sort(dist, t_dist)
        if lig == 'both':
            dist2 = dist2_a5_all[:,i]
            dist2_a5[:,i] = uncorr.sort(dist2, t_dist)

    for i in range(num_pairs_a6):
        dist = dist_a6_all[:,i]
        dist_a6[:,i] = uncorr.sort(dist, t_dist)
        if lig == 'both':
            dist2 = dist2_a6_all[:,i]
            dist2_a6[:,i] = uncorr.sort(dist2, t_dist)

    if traj_ns.n_residues > 297:
        for i in range(num_pairs_a7):
            dist = dist_a7_all[:,i]
            dist_a7[:,i] = uncorr.sort(dist, t_dist)
            if lig == 'both':
                dist2 = dist2_a7_all[:,i]
                dist2_a7[:,i] = uncorr.sort(dist2, t_dist)

    #Set array for the crystal structure binding location and two alternatives
    contact_loc1, contact_loc2, contact_loc3, contact_loc4, contact_unb = [],[],[],[],[]
    contact2_loc1, contact2_loc2, contact2_loc3, contact2_loc4, contact2_unb = [],[],[],[],[]

    #Open files for the number of interactions with each protein region at each point in time
    file_a3 = open('a3_inter.txt', 'w')
    file_a4 = open('a4_inter.txt', 'w')
    file_a5 = open('a5_inter.txt', 'w')
    file_a6 = open('a6_inter.txt', 'w')
    if traj_ns.n_residues > 298:
        file_a7 = open('a7_inter.txt', 'w')
    if lig == 'both':
        file2_a3 = open('a3_lig2_inter.txt', 'w')
        file2_a4 = open('a4_lig2_inter.txt', 'w')
        file2_a5 = open('a5_lig2_inter.txt', 'w')
        file2_a6 = open('a6_lig2_inter.txt', 'w')
        if traj_ns.n_residues > 298:
            file2_a7 = open('a7_lig2_inter.txt', 'w')

    #Loop through all frames
    for i in range(time_uncorr):
        #At each time point reset counters for helix interactions
        check_a3 = 0
        check_a4 = 0
        check_a5 = 0
        check_a6 = 0
        check_bend = 0
        check_a7 = 0
        if lig == 'both':
            check_a3 = 0
            check_a4 = 0
            check_a5 = 0
            check_a6 = 0
            check_bend = 0
            check_a7 = 0

        #Index for residues interactions in a3, a6, and a7 + the bend with ligand
        bond = 0
        for j in range(num_pairs_a3): #Determine # of contacts with a3 helix
            if dist_a3[i][j] <= 0.5:
                check_a3 += 1
                lig_tot_cont[bond] += 1
            if lig == 'both' and dist2_a3[i][j] <= 0.5:
                check2_a3 += 1
                lig2_tot_cont[bond] += 1
            bond += 1
        for k in range(num_pairs_a4): #Determine # of contacts with a4 helix
            if dist_a4[i][k] <= 0.5:
                check_a4 += 1
            if lig == 'both' and dist2_a4[i][k] <= 0.5:
                check2_a4 += 1
        for l in range(num_pairs_a5):  #Determine # of contacts with a5 helix
            if dist_a5[i][l] <= 0.5:
                check_a5 += 1
            if lig == 'both' and dist2_a5[i][l] <= 0.5:
                check2_a5 += 1
        for m in range(num_pairs_a6):  #Determine # of contacts with a6 helix
            if dist_a6[i][m] <= 0.5:
                check_a6 += 1
                lig_tot_cont[bond] += 1
            if lig == 'both' and dist2_a6[i][m] <= 0.5:
                check2_a6 += 1
                lig2_tot_cont[bond] += 1
            bond += 1
        for o in range(num_pairs_bend): #Determine the # of contacts between a6 and a7 helices
            if dist_bend[i][o] <= 0.5:
                check_bend += 1
                lig_tot_cont[bond] += 1
            if lig == 'both' and dist2_bend[i][o] <= 0.5:
                check2_bend += 1
                lig2_tot_cont[bond] += 1
            bond += 1
    
        #Output number of interactions to files for a3-a6 at each time point
        file_a3.write(str(check_a3) + '\n')
        file_a4.write(str(check_a4) + '\n')
        file_a5.write(str(check_a5) + '\n')
        file_a6.write(str(check_a6) + '\n')
        if lig == 'both':
            file2_a3.write(str(check2_a3) + '\n')
            file2_a4.write(str(check2_a4) + '\n')
            file2_a5.write(str(check2_a5) + '\n')
            file2_a6.write(str(check2_a6) + '\n')

        
        if traj_ns.n_residues == 300 or traj_ns.n_residues == 299:
            for n in range(num_pairs_a7):#Determine # of contacts with the a7 helix pt1
                if dist_a7[i][n] <= 0.5:
                    check_a7 += 1
                    lig_tot_cont[bond] += 1
                if lig == 'both' and dist2_a7[i][n] <= 0.5:
                    check2_a7 += 1
                    lig2_tot_cont[bond] += 1
                bond += 1
            #Output number of interactions to file for a7
            file_a7.write(str(check_a7) + '\n')
            if lig == 'both':
                file2_a7.write(str(check2_a7) + '\n')
            
            #Determine ligand binding Location
            #Sum total interactions for each binding location
            total_contact_loc1 = check_a3 + check_a6 + check_a7 #Crystal structure binding location
            total_contact_loc2 = check_a6 + check_bend + check_a7
            total_contact_loc3 = check_a4 + check_a5 + check_a6 
            total_contact_loc4 = check_a3 + check_a4 + check_a6
            if lig == 'both':
                total2_contact_loc1 = check2_a3 + check2_a6 + check2_a7 #Crystal structure binding location
                total2_contact_loc2 = check2_a6 + check2_bend + check2_a7
                total2_contact_loc3 = check2_a4 + check2_a5 + check2_a6 
                total2_contact_loc4 = check2_a3 + check2_a4 + check2_a6
           
            #Binding location 1
            if check_a3 >= 1 and check_a6 >= 1 and check_a7 >= 1 and total_contact_loc1 >= 4 and check_a5 == 0 and check_bend == 0:
                contact_loc1.append(1)
            else:
                contact_loc1.append(0)
            if lig == 'both':
                if check2_a3 >= 1 and check2_a6 >= 1 and check2_a7 >= 1 and total2_contact_loc1 >= 4 and check2_a5 == 0 and check2_bend == 0:
                    contact2_loc1.append(1)
                else:
                    contact2_loc1.append(0)

            #Binding Location 2
            if check_bend >= 1 and check_a6 >= 1 and check_a7 >= 1 and total_contact_loc2 >= 4 and check_a4 == 0 and check_a5 == 0:
                contact_loc2.append(1)
            else:
                contact_loc2.append(0)
            if lig == 'both':
                if check2_bend >= 1 and check2_a6 >= 1 and check2_a7 >= 1 and total2_contact_loc2 >= 4 and check2_a4 == 0 and check2_a5 == 0:
                    contact2_loc2.append(1)
                else:
                    contact2_loc2.append(0)

        else:
            #Sum total interactions without the a7 helix
            total_contact_loc1 = check_a3 + check_a6
            total_contact_loc2 = check_a6 + check_bend
            total_contact_loc3 = check_a4 + check_a5 + check_a6 
            total_contact_loc4 = check_a3 + check_a4 + check_a6
            if lig == 'both':
                total2_contact_loc1 = check2_a3 + check2_a6
                total2_contact_loc2 = check2_a6 + check2_bend
                total2_contact_loc3 = check2_a4 + check2_a5 + check2_a6 
                total2_contact_loc4 = check2_a3 + check2_a4 + check2_a6

            #Binding Location 1 w/o a7
            if check_a3 >= 1 and check_a6 >= 1 and total_contact_loc1 >= 5 and check_a5 == 0 and check_bend == 0 and check_a4 == 0:
                contact_loc1.append(1)
            else:
                contact_loc1.append(0)
            if lig == 'both':
                if check2_a3 >= 1 and check2_a6 >= 1 and total2_contact_loc1 >= 5 and check2_a5 == 0 and check2_bend == 0 and check2_a4 == 0:
                    contact2_loc1.append(1)
                else:
                    contact2_loc1.append(0)

            #Binding location 2 w/o a7
            if check_bend >= 1 and check_a6 >= 1 and total_contact_loc2 >= 5 and check_a4 == 0 and check_a5 == 0:
                contact_loc2.append(1)
            else:
                contact_loc2.append(0)
            if lig == 'both':
                if check2_bend >= 1 and check2_a6 >= 1 and total2_contact_loc2 >= 5 and check2_a4 == 0 and check2_a5 == 0:
                    contact2_loc2.append(1)
                else:
                    contact2_loc2.append(0)

        #Binding Location 3
        if check_a4 >= 1 and check_a6 >= 1 and total_contact_loc3 >= 5 and check_a3 == 0 and check_a7 == 0:
            contact_loc3.append(1)
        else:
            contact_loc3.append(0)
        if lig == 'both':
            if check2_a4 >= 1 and check2_a6 >= 1 and total2_contact_loc3 >= 5 and check2_a3 == 0 and check2_a7 == 0:
                contact2_loc3.append(1)
            else:
                contact2_loc3.append(0)

        #Binding Location 4
        if check_a3 >= 1 and check_a4 >= 1 and check_a6 >= 1 and total_contact_loc4 >= 5 and check_a5 == 0 and check_a7 == 0:
            contact_loc4.append(1)
        else:
            contact_loc4.append(0)
        if lig == 'both':
            if check2_a3 >= 1 and check2_a4 >= 1 and check2_a6 >= 1 and total2_contact_loc4 >= 5 and check2_a5 == 0 and check2_a7 == 0:
                contact2_loc4.append(1)
            else:
                contact2_loc4.append(0)

        if check_a3 == 0 and check_a6 == 0 and check_a7 == 0:
            contact_unb.append(1)
        else:
            contact_unb.append(0)
        if lig == 'both':
            if check2_a3 == 0 and check2_a6 == 0 and check2_a7 == 0:
                contact2_unb.append(1)
            else:
                contact2_unb.append(0)
    
    #Print % contact with each bond
    output_bond_cont = open('all_iter_frac.txt', 'w')
    if lig == 'both':
        output_bond_cont2 = open('all_iter_frac_lig2.txt', 'w')

    #Set array for names of residues interacting with ligand
    contacts = np.append(group_3, group_6)
    contacts = np.append(contacts, group_bend)
    contacts = np.append(contacts, group_7)
    for i in range(tot_pairs):
        output_bond_cont.write(str(contacts[i]) + ' ' + str(lig_tot_cont[i]/time_uncorr) + '\n') #Write name of residue and fraction of time in contact fot full trajectory to file
        if lig == 'both':
            output_bond_cont2.write(str(contacts[i]) + ' ' + str(lig2_tot_cont[i]/time_uncorr) + '\n') #Write name of residue and fraction of time in contact fot full trajectory to file

    output_bond_cont.close() #Close File
    output_bond_cont2.close() #Close File

    #Determine the fraction the ligand is in each location
    lig_frac_loc1 = sum(contact_loc1)/len(contact_loc1)
    lig_frac_loc3 = sum(contact_loc3)/len(contact_loc3)
    lig_frac_loc4 = sum(contact_loc4)/len(contact_loc4)
    lig_frac_loc2 = sum(contact_loc2)/len(contact_loc2)
    lig_unb = sum(contact_unb)/len(contact_unb)
    lig_other = 1 - lig_frac_loc1 - lig_frac_loc2 - lig_frac_loc3 - lig_frac_loc4 - lig_unb
    if lig == 'both':
        lig2_frac_loc1 = sum(contact2_loc1)/len(contact2_loc1)
        lig2_frac_loc3 = sum(contact2_loc3)/len(contact2_loc3)
        lig2_frac_loc4 = sum(contact2_loc4)/len(contact2_loc4)
        lig2_frac_loc2 = sum(contact2_loc2)/len(contact2_loc2)
        lig2_unb = sum(contact_unb)/len(contact_unb)
        lig2_other = 1 - lig2_frac_loc1 - lig2_frac_loc2 - lig2_frac_loc3 - lig2_frac_loc4 - lig2_unb

    #Print the percent time in each location
    Loc_frac = open('Lig_bind_frac' + File_base + '.txt', 'w')
    Loc_frac.write('Binding location 1: ' + str(100 * lig_frac_loc1) + '\n')
    Loc_frac.write('Binding location 2: ' + str(100 * lig_frac_loc2) + '\n')
    Loc_frac.write('Binding location 3: ' + str(100 * lig_frac_loc3) + '\n')
    Loc_frac.write('Binding location 4: ' + str(100 * lig_frac_loc4) + '\n')
    Loc_frac.write('Unbound: ' + str(100 * lig_unb) + '\n')
    Loc_frac.write('Other Bound: ' + str(100 * lig_other))
    Loc_frac.close() #Close file
    if lig == 'both':
        Loc_frac = open('Lig2_bind_frac' + File_base + '.txt', 'w')
        Loc_frac.write('Binding location 1: ' + str(100 * lig2_frac_loc1) + '\n')
        Loc_frac.write('Binding location 2: ' + str(100 * lig2_frac_loc2) + '\n')
        Loc_frac.write('Binding location 3: ' + str(100 * lig2_frac_loc3) + '\n')
        Loc_frac.write('Binding location 4: ' + str(100 * lig2_frac_loc4) + '\n')
        Loc_frac.write('Unbound: ' + str(100 * lig2_unb) + '\n')
        Loc_frac.write('Other Bound: ' + str(100 * lig2_other))
        Loc_frac.close() #Close file


    #Compute Simultaneous Ligand Contacts
    if traj_ns.n_residues == 300 or traj_ns.n_residues == 299:
        pairs = num_pairs_a3 + num_pairs_a6 + num_pairs_a7

    else:
        pairs = num_pairs_a3 + num_pairs_a6

    simul_contacts = np.zeros([pairs, pairs])
    simul_contacts2 = np.zeros([pairs, pairs])

    num_1 = 0
    num_2 = 0
    for i in range(num_pairs_a3):
        count = 0
        count2 = 0
        num_2 = 0
        for j in range(num_pairs_a3):
            for t in range(time_uncorr):
                if dist_a3[t][i] <= 0.5 and dist_a3[t][j] <= 0.5:
                    count += 1
                if lig == 'both':
                    if dist2_a3[t][i] <= 0.5 and dist2_a3[t][j] <= 0.5:
                        count2 += 1
            simul_contacts[num_1][num_2] = count/time_uncorr
            if lig == 'both':
                simul_contacts2[num_1][num_2] = count2/time_uncorr

            count = 0
            count2 = 0
            num_2 += 1
        for k in range(num_pairs_a6):
            for t in range(time_uncorr):
                if dist_a3[t][i] <= 0.5 and dist_a6[t][k] <= 0.5:
                    count += 1
                if lig == 'both':
                    if dist2_a3[t][i] <= 0.5 and dist2_a6[t][k] <= 0.5:
                        count2 += 1
            simul_contacts[num_1][num_2] = count/time_uncorr
            if lig == 'both':
                simul_contacts2[num_1][num_2] = count2/time_uncorr

            count = 0
            count2 = 0
            num_2 += 1
        if traj_ns.n_residues == 300 or traj_ns.n_residues == 299:
            for l in range(num_pairs_a7):
                for t in range(time_uncorr):
                    if dist_a3[t][i] <= 0.5 and dist_a7[t][l] <= 0.5:
                        count += 1
                    if lig == 'both':
                        if dist2_a3[t][i] <= 0.5 and dist2_a7[t][l] <= 0.5:
                            count2 += 1
                simul_contacts[num_1][num_2] = count/time_uncorr
                if lig == 'both':
                    simul_contacts2[num_1][num_2] = count2/time_uncorr

                count = 0
                count2 = 0
                num_2 += 1
        num_1 += 1
    for i in range(num_pairs_a6):
        count = 0
        count2 = 0
        num_2 = 0
        for j in range(num_pairs_a3):
            for t in range(time_uncorr):
                if dist_a6[t][i] <= 0.5 and dist_a3[t][j] <= 0.5:
                    count += 1
                if lig == 'both':
                    if dist2_a6[t][i] <= 0.5 and dist2_a3[t][j] <= 0.5:
                        count2 += 1
            simul_contacts[num_1][num_2] = count/time_uncorr
            simul_contacts2[num_1][num_2] = count2/time_uncorr

            count = 0
            count2 = 0
            num_2 += 1
        for k in range(num_pairs_a6):
            for t in range(time_uncorr):
                if dist_a6[t][i] <= 0.5 and dist_a6[t][k] <= 0.5:
                    count += 1
                if lig == 'both':
                    if dist2_a6[t][i] <= 0.5 and dist2_a6[t][k] <= 0.5:
                        count2 += 1
            simul_contacts[num_1][num_2] = count/time_uncorr
            simul_contacts2[num_1][num_2] = count2/time_uncorr

            count = 0
            count2 = 0
            num_2 += 1
        if traj_ns.n_residues == 300 or traj_ns.n_residues == 299:
            for l in range(num_pairs_a7):
                for t in range(time_uncorr):
                    if dist_a6[t][i] <= 0.5 and dist_a7[t][l] <= 0.5:
                        count += 1
                    if lig == 'both':
                        if dist2_a6[t][i] <= 0.5 and dist2_a7[t][l] <= 0.5:
                            count2 += 1
                simul_contacts[num_1][num_2] = count/time_uncorr
                if lig == 'both':
                    simul_contacts2[num_1][num_2] = count2/time_uncorr

                count = 0
                count2 = 0
                num_2 += 1
        num_1 += 1
    if traj_ns.n_residues == 300 or traj_ns.n_residues == 299:
        for i in range(num_pairs_a7):
            count = 0
            count2 = 0
            num_2 = 0
            for j in range(num_pairs_a3):
                for t in range(time_uncorr):
                    if dist_a7[t][i] <= 0.5 and dist_a3[t][j] <= 0.5:
                        count += 1
                    if lig == 'both':
                        if dist2_a7[t][i] <= 0.5 and dist2_a3[t][j] <= 0.5:
                            count2 += 1
                simul_contacts[num_1][num_2] = count/time_uncorr
                if lig == 'both':
                    simul_contacts2[num_1][num_2] = count2/time_uncorr

                count = 0
                count2 = 0
                num_2 += 1
            for k in range(num_pairs_a6):
                for t in range(time_uncorr):
                    if dist_a7[t][i] <= 0.5 and dist_a6[t][k] <= 0.5:
                        count += 1
                simul_contacts[num_1][num_2] = count/time_uncorr
    
                count = 0
                count2 = 0
                num_2 += 1
            for l in range(num_pairs_a7):
                for t in range(time_uncorr):
                    if dist_a7[t][i] <= 0.5 and dist_a7[t][l] <= 0.5:
                        count += 1
                    if lig == 'both':
                        if dist2_a7[t][i] <= 0.5 and dist2_a7[t][l] <= 0.5:
                            count2 += 1
                simul_contacts[num_1][num_2] = count/time_uncorr
                if lig == 'both':
                    simul_contacts2[num_1][num_2] = count2/time_uncorr

                count = 0
                count2 = 0
                num_2 += 1
            num_1 += 1
    np.savetxt('simul_lig_contact_' + File_base + '.txt', simul_contacts)
    if lig == 'both':
        np.savetxt('simul_lig2_contact_' + File_base + '.txt', simul_contacts2)

    #Simultaneour Ligand Contacts for Sections
    lig_cont_sect = np.zeros([9,9])
    if lig == 'both':
        lig2_cont_sect = np.zeros([9,9])

    for t in range(time_uncorr):
        #Residue 186 to 192
        if dist_a3[t][0] <= 0.5 or dist_a3[t][1] <= 0.5 or dist_a3[t][2] <= 0.5 or dist_a3[t][3] <= 0.5 or dist_a3[t][4] <= 0.5 or dist_a3[t][5] <= 0.5 or dist_a3[t][6] <= 0.5:
            a3_sect1 = 1
        else:
            a3_sect1 = 0
        #Residue 193 to 196
        if dist_a3[t][7] <= 0.5 or dist_a3[t][8] <= 0.5 or dist_a3[t][9] <= 0.5 or dist_a3[t][10] <= 0.5:
            a3_sect2 = 1
        else:
            a3_sect2 = 0
        #Residue 197 to 200
        if dist_a3[t][11] <= 0.5 or dist_a3[t][12] <= 0.5 or dist_a3[t][13] <= 0.5 or dist_a3[t][14] <= 0.5:
            a3_sect3 = 1
        else:
            a3_sect3 = 0
        #Residue 264 to 270
        if dist_a6[t][0] <= 0.5 or dist_a6[t][1] <= 0.5 or dist_a6[t][2] <= 0.5 or dist_a6[t][3] <= 0.5 or dist_a6[t][4] <= 0.5 or dist_a6[t][5] <= 0.5 or dist_a6[t][6] <= 0.5: 
            a6_sect1 = 1
        else:
            a6_sect1 = 0
        #Residue 193 to 196
        if dist_a6[t][7] <= 0.5 or dist_a6[t][8] <= 0.5 or dist_a6[t][9] <= 0.5 or dist_a6[t][10] <= 0.5 or dist_a6[t][11] <= 0.5:
            a6_sect2 = 1
        else:
            a6_sect2 = 0
        #Residue 197 to 200
        if dist_a6[t][12] <= 0.5 or dist_a6[t][13] <= 0.5 or dist_a6[t][14] <= 0.5 or dist_a6[t][15] <= 0.5 or dist_a6[t][16] <= 0.5: 
            a6_sect3 = 1
        else:
            a6_sect3 = 0
        if traj_ns.n_residues == 300 or traj_ns.n_residues == 299:
            #Residue 287 to 290
            if dist_a7[t][0] <= 0.5 or dist_a7[t][1] <= 0.5 or dist_a7[t][2] <= 0.5 or dist_a7[t][3] <= 0.5:
                a7_sect1 = 1
            else:
                a7_sect1 = 0
            #Residue 291 to 294
            if dist_a7[t][4] <= 0.5 or dist_a7[t][5] <= 0.5 or dist_a7[t][6] <= 0.5 or dist_a7[t][7] <= 0.5: 
                a7_sect2 = 1
            else:
                a7_sect2 = 0
            #Residue 295 to 298
            if dist_a7[t][8] <= 0.5 or dist_a7[t][9] <= 0.5 or dist_a7[t][10] <= 0.5 or dist_a7[t][11] <= 0.5:
                a7_sect3 = 1
            else:
                a7_sect3 = 0
        #Make matrix for simultaneous contacts
        if a3_sect1 == 1:
            lig_cont_sect[0][0] += 1
        if a3_sect1 == 1 and a3_sect2 == 1:
            lig_cont_sect[0][1] += 1
            lig_cont_sect[1][0] += 1
        if a3_sect1 == 1 and a3_sect3 == 1:
            lig_cont_sect[0][2] += 1
            lig_cont_sect[2][0] += 1
        if a3_sect1 == 1 and a6_sect1 == 1:
            lig_cont_sect[0][3] += 1
            lig_cont_sect[3][0] += 1
        if a3_sect1 == 1 and a6_sect2 == 1:
            lig_cont_sect[0][4] += 1
            lig_cont_sect[4][0] += 1
        if a3_sect1 == 1 and a6_sect3 == 1:
            lig_cont_sect[0][5] += 1
            lig_cont_sect[5][0] += 1
        if traj_ns.n_residues == 300 or traj_ns.n_residues == 299:
            if a3_sect1 == 1 and a7_sect1 == 1:
                lig_cont_sect[0][6] += 1
                lig_cont_sect[6][0] += 1
            if a3_sect1 == 1 and a7_sect2 == 1:
                lig_cont_sect[0][7] += 1
                lig_cont_sect[7][0] += 1
            if a3_sect1 == 1 and a7_sect3 == 1:
                lig_cont_sect[0][8] += 1
                lig_cont_sect[8][0] += 1
        if a3_sect2 == 1:
            lig_cont_sect[1][1] += 1
        if a3_sect2 == 1 and a3_sect3 == 1:
            lig_cont_sect[1][2] += 1
            lig_cont_sect[2][1] += 1
        if a3_sect2 == 1 and a6_sect1 == 1:
            lig_cont_sect[1][3] += 1
            lig_cont_sect[3][1] += 1
        if a3_sect2 == 1 and a6_sect2 == 1:
            lig_cont_sect[1][4] += 1
            lig_cont_sect[4][1] += 1
        if a3_sect2 == 1 and a6_sect3 == 1:
            lig_cont_sect[1][5] += 1
            lig_cont_sect[5][1] += 1
        if traj_ns.n_residues == 300 or traj_ns.n_residues == 299:
            if a3_sect2 == 1 and a7_sect1 == 1:
                lig_cont_sect[1][6] += 1
                lig_cont_sect[6][1] += 1
            if a3_sect2 == 1 and a7_sect2 == 1:
                lig_cont_sect[1][7] += 1
                lig_cont_sect[7][1] += 1
            if a3_sect2 == 1 and a7_sect3 == 1:
                lig_cont_sect[1][8] += 1
                lig_cont_sect[8][1] += 1
        if a3_sect3 == 1:
            lig_cont_sect[2][2] += 1
        if a3_sect3 == 1 and a6_sect1 == 1:
            lig_cont_sect[2][3] += 1
            lig_cont_sect[3][2] += 1
        if a3_sect3 == 1 and a6_sect2 == 1:
            lig_cont_sect[2][4] += 1
            lig_cont_sect[4][2] += 1
        if a3_sect3 == 1 and a6_sect3 == 1:
            lig_cont_sect[2][5] += 1
            lig_cont_sect[5][2] += 1
        if traj_ns.n_residues == 300 or traj_ns.n_residues == 299:
            if a3_sect3 == 1 and a7_sect1 == 1:
                lig_cont_sect[2][6] += 1
                lig_cont_sect[6][2] += 1
            if a3_sect3 == 1 and a7_sect2 == 1:
                lig_cont_sect[2][7] += 1
                lig_cont_sect[7][2] += 1
            if a3_sect3 == 1 and a7_sect3 == 1:
                lig_cont_sect[2][8] += 1
                lig_cont_sect[8][2] += 1
        if a6_sect1 == 1:
            lig_cont_sect[3][3] += 1
        if a6_sect1 == 1 and a6_sect2 == 1:
            lig_cont_sect[3][4] += 1
            lig_cont_sect[4][3] += 1
        if a6_sect1 == 1 and a6_sect3 == 1:
            lig_cont_sect[3][5] += 1
            lig_cont_sect[5][3] += 1
        if traj_ns.n_residues == 300 or traj_ns.n_residues == 299:
            if a6_sect1 == 1 and a7_sect1 == 1:
                lig_cont_sect[3][6] += 1
                lig_cont_sect[6][3] += 1
            if a6_sect1 == 1 and a7_sect2 == 1:
                lig_cont_sect[3][7] += 1
                lig_cont_sect[7][3] += 1
            if a6_sect1 == 1 and a7_sect3 == 1:
                lig_cont_sect[3][8] += 1
                lig_cont_sect[8][3] += 1
        if a6_sect2 == 1:
            lig_cont_sect[4][4] += 1
        if a6_sect2 == 1 and a6_sect3 == 1:
            lig_cont_sect[4][5] += 1
            lig_cont_sect[5][4] += 1
        if traj_ns.n_residues == 300 or traj_ns.n_residues == 299:
            if a6_sect2 == 1 and a7_sect1 == 1:
                lig_cont_sect[4][6] += 1
                lig_cont_sect[6][4] += 1
            if a6_sect2 == 1 and a7_sect2 == 1:
                lig_cont_sect[4][7] += 1
                lig_cont_sect[7][4] += 1
            if a6_sect2 == 1 and a7_sect3 == 1:
                lig_cont_sect[4][8] += 1
                lig_cont_sect[8][4] += 1
        if a6_sect3 == 1:
            lig_cont_sect[5][5] += 1
        if traj_ns.n_residues == 300 or traj_ns.n_residues == 299:
            if a6_sect3 == 1 and a7_sect1 == 1:
                lig_cont_sect[5][6] += 1
                lig_cont_sect[6][5] += 1
            if a6_sect3 == 1 and a7_sect2 == 1:
                lig_cont_sect[5][7] += 1
                lig_cont_sect[7][5] += 1
            if a6_sect3 == 1 and a7_sect3 == 1:
                lig_cont_sect[5][8] += 1
                lig_cont_sect[8][5] += 1
        if traj_ns.n_residues == 300 or traj_ns.n_residues == 299:
            if a7_sect1 == 1:
                lig_cont_sect[6][6] += 1
            if a7_sect1 == 1 and a7_sect2 == 1:
                lig_cont_sect[6][7] += 1
                lig_cont_sect[7][6] += 1
            if a7_sect1 == 1 and a7_sect3 == 1:
                lig_cont_sect[6][8] += 1
                lig_cont_sect[8][6] += 1
            if a7_sect2 == 1:
                lig_cont_sect[7][7] += 1
            if a7_sect2 == 1 and a7_sect3 == 1:
                lig_cont_sect[7][8] += 1
                lig_cont_sect[8][7] += 1
            if a7_sect3 == 1:
                lig_cont_sect[8][8] += 1
        if lig == 'both':
            #Residue 186 to 192
            if dist2_a3[t][0] <= 0.5 or dist2_a3[t][1] <= 0.5 or dist2_a3[t][2] <= 0.5 or dist2_a3[t][3] <= 0.5 or dist2_a3[t][4] <= 0.5 or dist2_a3[t][5] <= 0.5 or dist2_a3[t][6] <= 0.5:
                a3_sect1 = 1
            else:
                a3_sect1 = 0
            #Residue 193 to 196
            if dist2_a3[t][7] <= 0.5 or dist2_a3[t][8] <= 0.5 or dist2_a3[t][9] <= 0.5 or dist2_a3[t][10] <= 0.5:
               a3_sect2 = 1
            else:
                a3_sect2 = 0
            #Residue 197 to 200
            if dist2_a3[t][11] <= 0.5 or dist2_a3[t][12] <= 0.5 or dist2_a3[t][13] <= 0.5 or dist2_a3[t][14] <= 0.5:
                a3_sect3 = 1
            else:
                a3_sect3 = 0
            #Residue 264 to 270
            if dist2_a6[t][0] <= 0.5 or dist2_a6[t][1] <= 0.5 or dist2_a6[t][2] <= 0.5 or dist2_a6[t][3] <= 0.5 or dist2_a6[t][4] <= 0.5 or dist2_a6[t][5] <= 0.5 or dist2_a6[t][6] <= 0.5: 
                a6_sect1 = 1
            else:
                a6_sect1 = 0
            #Residue 193 to 196
            if dist2_a6[t][7] <= 0.5 or dist2_a6[t][8] <= 0.5 or dist2_a6[t][9] <= 0.5 or dist2_a6[t][10] <= 0.5 or dist2_a6[t][11] <= 0.5:
                a6_sect2 = 1
            else:
                a6_sect2 = 0
            #Residue 197 to 200
            if dist2_a6[t][12] <= 0.5 or dist2_a6[t][13] <= 0.5 or dist2_a6[t][14] <= 0.5 or dist2_a6[t][15] <= 0.5 or dist2_a6[t][16] <= 0.5: 
                a6_sect3 = 1
            else:
                a6_sect3 = 0
            if traj_ns.n_residues == 300 or traj_ns.n_residues == 299:
                #Residue 287 to 290
                if dist2_a7[t][0] <= 0.5 or dist2_a7[t][1] <= 0.5 or dist2_a7[t][2] <= 0.5 or dist2_a7[t][3] <= 0.5:
                    a7_sect1 = 1
                else:
                    a7_sect1 = 0
                #Residue 291 to 294
                if dist2_a7[t][4] <= 0.5 or dist2_a7[t][5] <= 0.5 or dist2_a7[t][6] <= 0.5 or dist2_a7[t][7] <= 0.5: 
                    a7_sect2 = 1
                else:
                    a7_sect2 = 0
                #Residue 295 to 298
                if dist2_a7[t][8] <= 0.5 or dist2_a7[t][9] <= 0.5 or dist2_a7[t][10] <= 0.5 or dist2_a7[t][11] <= 0.5:
                    a7_sect3 = 1
                else:
                    a7_sect3 = 0
            #Make matrix for simultaneous contacts
            if a3_sect1 == 1:
                lig2_cont_sect[0][0] += 1
            if a3_sect1 == 1 and a3_sect2 == 1:
                lig2_cont_sect[0][1] += 1
                lig2_cont_sect[1][0] += 1
            if a3_sect1 == 1 and a3_sect3 == 1:
                lig2_cont_sect[0][2] += 1
                lig2_cont_sect[2][0] += 1
            if a3_sect1 == 1 and a6_sect1 == 1:
                lig2_cont_sect[0][3] += 1
                lig2_cont_sect[3][0] += 1
            if a3_sect1 == 1 and a6_sect2 == 1:
                lig2_cont_sect[0][4] += 1
                lig2_cont_sect[4][0] += 1
            if a3_sect1 == 1 and a6_sect3 == 1:
                lig2_cont_sect[0][5] += 1
                lig2_cont_sect[5][0] += 1
            if traj_ns.n_residues == 300 or traj_ns.n_residues == 299:
                if a3_sect1 == 1 and a7_sect1 == 1:
                    lig2_cont_sect[0][6] += 1
                    lig2_cont_sect[6][0] += 1
                if a3_sect1 == 1 and a7_sect2 == 1:
                    lig2_cont_sect[0][7] += 1
                    lig2_cont_sect[7][0] += 1
                if a3_sect1 == 1 and a7_sect3 == 1:
                    lig2_cont_sect[0][8] += 1
                    lig2_cont_sect[8][0] += 1
            if a3_sect2 == 1:
                lig2_cont_sect[1][1] += 1
            if a3_sect2 == 1 and a3_sect3 == 1:
                lig2_cont_sect[1][2] += 1
                lig2_cont_sect[2][1] += 1
            if a3_sect2 == 1 and a6_sect1 == 1:
                lig2_cont_sect[1][3] += 1
                lig2_cont_sect[3][1] += 1
            if a3_sect2 == 1 and a6_sect2 == 1:
                lig2_cont_sect[1][4] += 1
                lig2_cont_sect[4][1] += 1
            if a3_sect2 == 1 and a6_sect3 == 1:
                lig2_cont_sect[1][5] += 1
                lig2_cont_sect[5][1] += 1
            if traj_ns.n_residues == 300 or traj_ns.n_residues == 299:
                if a3_sect2 == 1 and a7_sect1 == 1:
                    lig2_cont_sect[1][6] += 1
                    lig2_cont_sect[6][1] += 1
                if a3_sect2 == 1 and a7_sect2 == 1:
                    lig2_cont_sect[1][7] += 1
                    lig2_cont_sect[7][1] += 1
                if a3_sect2 == 1 and a7_sect3 == 1:
                    lig2_cont_sect[1][8] += 1
                    lig2_cont_sect[8][1] += 1
            if a3_sect3 == 1:
                lig2_cont_sect[2][2] += 1
            if a3_sect3 == 1 and a6_sect1 == 1:
                lig2_cont_sect[2][3] += 1
                lig2_cont_sect[3][2] += 1
            if a3_sect3 == 1 and a6_sect2 == 1:
                lig2_cont_sect[2][4] += 1
                lig2_cont_sect[4][2] += 1
            if a3_sect3 == 1 and a6_sect3 == 1:
                lig2_cont_sect[2][5] += 1
                lig2_cont_sect[5][2] += 1
            if traj_ns.n_residues == 300 or traj_ns.n_residues == 299:
                if a3_sect3 == 1 and a7_sect1 == 1:
                    lig2_cont_sect[2][6] += 1
                    lig2_cont_sect[6][2] += 1
                if a3_sect3 == 1 and a7_sect2 == 1:
                    lig2_cont_sect[2][7] += 1
                    lig2_cont_sect[7][2] += 1
                if a3_sect3 == 1 and a7_sect3 == 1:
                    lig2_cont_sect[2][8] += 1
                    lig2_cont_sect[8][2] += 1
            if a6_sect1 == 1:
                lig2_cont_sect[3][3] += 1
            if a6_sect1 == 1 and a6_sect2 == 1:
                lig2_cont_sect[3][4] += 1
                lig2_cont_sect[4][3] += 1
            if a6_sect1 == 1 and a6_sect3 == 1:
                lig2_cont_sect[3][5] += 1
                lig2_cont_sect[5][3] += 1
            if traj_ns.n_residues == 300 or traj_ns.n_residues == 299:
                if a6_sect1 == 1 and a7_sect1 == 1:
                    lig2_cont_sect[3][6] += 1
                    lig2_cont_sect[6][3] += 1
                if a6_sect1 == 1 and a7_sect2 == 1:
                    lig2_cont_sect[3][7] += 1
                    lig2_cont_sect[7][3] += 1
                if a6_sect1 == 1 and a7_sect3 == 1:
                    lig2_cont_sect[3][8] += 1
                    lig2_cont_sect[8][3] += 1
            if a6_sect2 == 1:
                lig2_cont_sect[4][4] += 1
            if a6_sect2 == 1 and a6_sect3 == 1:
                lig2_cont_sect[4][5] += 1
                lig2_cont_sect[5][4] += 1
            if traj_ns.n_residues == 300 or traj_ns.n_residues == 299:
                if a6_sect2 == 1 and a7_sect1 == 1:
                    lig2_cont_sect[4][6] += 1
                    lig2_cont_sect[6][4] += 1
                if a6_sect2 == 1 and a7_sect2 == 1:
                    lig2_cont_sect[4][7] += 1
                    lig2_cont_sect[7][4] += 1
                if a6_sect2 == 1 and a7_sect3 == 1:
                    lig2_cont_sect[4][8] += 1
                    lig2_cont_sect[8][4] += 1
            if a6_sect3 == 1:
                lig2_cont_sect[5][5] += 1
            if traj_ns.n_residues == 300 or traj_ns.n_residues == 299:
                if a6_sect3 == 1 and a7_sect1 == 1:
                    lig2_cont_sect[5][6] += 1
                    lig2_cont_sect[6][5] += 1
                if a6_sect3 == 1 and a7_sect2 == 1:
                    lig2_cont_sect[5][7] += 1
                    lig2_cont_sect[7][5] += 1
                if a6_sect3 == 1 and a7_sect3 == 1:
                    lig2_cont_sect[5][8] += 1
                    lig2_cont_sect[8][5] += 1
            if traj_ns.n_residues == 300 or traj_ns.n_residues == 299:
                if a7_sect1 == 1:
                    lig2_cont_sect[6][6] += 1
                if a7_sect1 == 1 and a7_sect2 == 1:
                    lig2_cont_sect[6][7] += 1
                    lig2_cont_sect[7][6] += 1
                if a7_sect1 == 1 and a7_sect3 == 1:
                    lig2_cont_sect[6][8] += 1
                    lig2_cont_sect[8][6] += 1
                if a7_sect2 == 1:
                    lig2_cont_sect[7][7] += 1
                if a7_sect2 == 1 and a7_sect3 == 1:
                    lig2_cont_sect[7][8] += 1
                    lig2_cont_sect[8][7] += 1
                if a7_sect3 == 1:
                    lig2_cont_sect[8][8] += 1
    lig_cont_sect_per = 100 * (lig_cont_sect/time_uncorr)
    if lig == 'both':
        lig2_cont_sect_per = 100 * (lig2_cont_sect/time_uncorr)
        np.savetxt('simul_lig2_contact_sect' + File_base + '.txt', lig2_cont_sect_per)

    np.savetxt('simul_lig_contact_sect' + File_base + '.txt', lig_cont_sect_per)
    print('Ligand Interaction Analysis Complete')
else:
    print('Ligand Interaction Analysis Skipped')

if pca_ck == True:
    #Principle Component Analysis
    pca = PCA(n_components=10) #initialize
    traj_bb.superpose(traj_bb, 0) #align trajectory
    reduced_cartesian = pca.fit_transform(traj_bb.xyz.reshape(traj_bb.n_frames, traj_bb.n_atoms * 3))
    per_var = np.round(pca.explained_variance_ratio_ * 100, decimals =2)
    accum_per_var = [ i for i in [ np . sum ( per_var [: j ]) for j in range (1 ,11)]]

    #Plot PCA
    fig = plt.figure()
    plt.scatter(reduced_cartesian[:, 0], reduced_cartesian[:,1], marker='x', c=traj.time)
    plt.xlabel('PC1')
    plt.ylabel('PC2')
    plt.title('Cartesian coordinate PCA for '+File_base)
    cbar = plt.colorbar()
    cbar.set_label('Time [ps]')
    fig.savefig('PCA_'+File_base+'.png')
    plt.close(fig)

    fig.ax = plt.figure(figsize=(10,5.5))

    #Scree Plot
    labels = [ 'PC' + str(x) for x in range(1 , len(per_var) +1) ]
    ax1 = plt.subplot (1 ,2 ,1)
    ax2 = ax1.twinx ()
    ax1.bar( x = range(1 , len(per_var) +1) , height = per_var , tick_label = labels , alpha = 0.85)
    ax2.plot(range(1 , len(per_var) +1), accum_per_var, color = 'r' , marker = 'o')
    ax2.grid( True )
    xlocs , xlabs = plt.xticks()
    ax1.set_ylabel ( 'Percentage of explained variance (%)', size = '12')
    ax2.set_ylabel ( 'Accumulated explained variance (%)', size = '12')
    ax1.set_xlabel ( 'Principal Components' , size = '12')
    ax1.set_ylim ([0 , max(per_var)*1.1 ] )
    ax2.set_ylim ([0 , max(accum_per_var)*1.1 ] )
    plt.title ( 'Scree Plot' , size = '14')
    plt.grid ( True )
    plt.savefig('Scree_'+File_base+'.png')
    print('PCA Completed')

if pca_ck == False:
    print('PCA Skipped')


