# Allosteric Inhibition of PTP1B by a Nonpolar Terpenoid
## Input Files
This repository contains the MDP files with all settings used to run the simulations included within the paper. The depository also contains the TPR files necessary to replicated the production simulations used for the paper. All data files are held within the data diectory. 

## Analysis Scripts
All analysis used for this paper can be regenerated using the files containded within the analysis_scripts directory. A guide to the files found within this direcory can be found below:

traj_analysis.py: This is the primary analysis script used to generate a majority of the data presented within the paper. This script takes a GROMACS trajectory and gro file as well as options specifying the analysis which should be done. The "-h" flag can be used when running the script for syntax and input options. This script will compute RMSD and RMSF from a specified reference for the full protein and important sections for PTP1B, utilize DSSP to identify the secondary structure of all residues within the $\alpha$7 helix of PTP1B, Identify all h-bonds formed between residues in the protein for over 60% of the trajectory, identify all h-bonds formed between the protein and a specificed ligand over 10% of the trajectory, identify residue contacts fromed by a ligand with a protein, identify binding conformation for AD, compute the COM RMSD for the ligand, identify interactions between the residues within the helical triad of the protein, and run PCA.

equil_deter.py: Equilibration determination for a GROMACS trajectory based on stabilization of protein backbone RMSD.

pdb_analysis.py: Analysis of the confrmation of the WPD loop and radius of gyration for the ligand from a PDB file.

hbond_Apo_analysis.py: Determination of the % of the trajectory a h-bond was fromed from an input list.

Hbond_corr.py: Determination of the % of the trajectory that a given h-bond was formed as well as the % of a trajectory that two h-bonds were formed simultaneously from a given input list.

Hbond_sect.py: Determine the number of h-bonds formed between section A and secion B of a protein as defined by the input.

Hbonds_note.txt: Hbonds taken from previous literature that were shown to correlate to various PTP1B behaviors (ie. WPD loop conformation or inhibtion).

Hbonds_note_atom.txt: Atom indices for the Hbonds specified above.
