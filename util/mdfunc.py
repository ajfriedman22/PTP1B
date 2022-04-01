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

    return traj, traj_bb, traj_prot, traj_ns
