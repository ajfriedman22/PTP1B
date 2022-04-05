#!/ usr / bin / env python

#Load MDTraj trajectory and process to output one full trajectory, one with no solvent, one with all protein residues, and one with protein bb atoms only
def mdtraj_load(File_traj, File_gro):
    #Import required packages
    import mdtraj as md

    #Load trajectories
    traj = md.load(File_traj, top=File_gro)
    top = traj.topology

    #Process Trajectory
    traj_bb = traj.atom_slice(top.select('backbone')) #Backbond atoms of PTP1B only
    traj_prot = traj.atom_slice(top.select('protein')) #Select only atoms in the protein
    traj_ns = traj.remove_solvent() #Remove solvent from the trajectory leaving only protein (and ligand if applicable)

    #Detemine if first residue is missing
    if traj_prot.n_residues == 299:
        miss_first = False
    else:
        miss_first = True

    if miss_first == False:
        traj_a7 = traj.atom_slice(top.select('286 <= resid and resid <= 294')) #Select only atoms in the a7 helix
    else:
        traj_a7 = traj.atom_slice(top.select('285 <= resid and resid <= 293')) #Select only atoms in the a7 helix

    return traj, traj_bb, traj_prot, traj_ns, traj_a7, miss_first

#Function to set arrays for all PTP1B protein sections
def set_sect(miss_first, lig):
    import numpy as np
    #MDtraj numbers residues starting with zero. Add an offset to ensure that the same residues are compared between trajectories
    if miss_first == True:
        offset = 1
    else:
        offset = 0
    group_WPD = np.linspace(176-offset, 184-offset, num = 9) #residues in the WPD loop
    group_3 = np.linspace(185-offset, 199-offset, num = 15) #residues in the a3 helix
    group_4 = np.linspace(220-offset, 237-offset, num = 18) #residues in the a4 helix
    group_5 = np.linspace(244-offset, 251-offset, num = 8) #residues in the a5 helix
    group_6 = np.linspace(263-offset, 280-offset, num = 18) #residues in the a6 helix
    group_bend = np.linspace(281-offset, 285-offset, num = 5) #residues in the bend b/w the a6 and a7 helices
    group_7 = np.linspace(286-offset, 297-offset, num = 12) #residues in the a7 helix
    group_L11 = np.linspace(149-offset, 152-offset, num = 4) #residues in the L11 loop
    pair_other = [[150-offset, 190-offset], [263-offset, 184-offset], [149-offset, 177-offset], [80-offset, 198-offset], [116-offset, 181-offset], [116-offset, 216-offset], [191-offset, 224-offset], 
            [181-offset, 219-offset], [186-offset, 268-offset], [177-offset, 189-offset], [151-offset, 176-offset]] #Residue pairs for distances not included above
    
    #Set the index for the ligand if present in trajectory
    if lig != 'none':
        if traj_prot.n_residues == 299: #w/ a7 helix
            if lig == 'both':
                group_l = [299]
                group_l2 = [300]
            else:
                group_l = [299]
        elif traj_prot.n_residues == 287: #w/o a7 helix
            if lig == 'both':
                group_l = [287]
                group_l2 = [288]
            else:
                group_l = [287]
        if lig == 'both':
            return group_WPD, group_3, group_4, group_6, group_7, group_L11, pair_other, group_l, group_l2
        else:
            return group_WPD, group_3, group_4, group_6, group_7, group_L11, pair_other, group_l
    else:
        return group_WPD, group_3, group_4, group_6, group_7, group_L11, pair_other

#Function computes the % of time there are simultaneous contacts between the ligand and residues in pairs A with those in pairs A, B, and C 
def compute_simul_comtacts(pairs_A, pairs_B, pairs_C, time_uncorr, dist_A, dist_B, dist_C, num_1, simul_contacts):
    for i in range(pairs_A):
        count = 0
        count2 = 0
        num_2 = 0
        for j in range(pairs_A):
            for t in range(time_uncorr):
                if dist_A[t][i] <= 0.5 and dist_A[t][j] <= 0.5:
                    count += 1
            simul_contacts[num_1][num_2] = count/time_uncorr

            count = 0
            count2 = 0
            num_2 += 1
        for k in range(pairs_B):
            for t in range(time_uncorr):
                if dist_A[t][i] <= 0.5 and dist_B[t][k] <= 0.5:
                    count += 1
            simul_contacts[num_1][num_2] = count/time_uncorr

            count = 0
            count2 = 0
            num_2 += 1
        if a7 == True:
            for l in range(pairs_C):
                for t in range(time_uncorr):
                    if dist_A[t][i] <= 0.5 and dist_C[t][l] <= 0.5:
                        count += 1
                simul_contacts[num_1][num_2] = count/time_uncorr

                count = 0
                count2 = 0
                num_2 += 1
        num_1 += 1
    return simul_contacts, num1

#Function to determine if there are any contacts formed within the specified residue section
def sect_contact(dist, t, low, high):
    dist_sect = dist[t][low:high]
    if min(dist_sect) <= 0.5:
        return 1
    else:
        return 0

