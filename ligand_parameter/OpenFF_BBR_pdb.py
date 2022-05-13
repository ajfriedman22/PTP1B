# -*- coding: utf-8 -*-
"""
Created on Mon Jan 25 10:14:00 2021

@author: Anika
"""
#Loading OpenFF Forcefield and parmed
from openff.toolkit.typing.engines.smirnoff.forcefield import ForceField
from simtk.openmm.app import PDBFile
import parmed
from simtk import unit
import numpy as np

parsley = ForceField('openff_unconstrained-1.3.0.offxml')

#Loading up an OpenFF Topology from amorphadiene's SMILES
from openff.toolkit.topology import Molecule, Topology

BBR = Molecule.from_smiles('CCC1=C(C2=C(O1)C=C(C=C2)S(=O)(=O)NC3=CC=C(C=C3)S(=O)(=O)NC4=NC=CS4)C(=O)C5=CC(=C(C(=C5)Br)O)Br', allow_undefined_stereo=True)
top_ff = Topology.from_molecules(1*[BBR])
#Import PDB file
pdbfile = PDBFile('BBR_processed.pdb')

#Generate topology
top_mm = top_ff.to_openmm()

#Set box vectors for OpenFF topology
top_ff.box_vectors = np.array([4, 4, 4]) * unit.nanometer

#Create an OpenMM System
openmm_sys = parsley.create_openmm_system(top_ff)

# Convert OpenMM System to a ParmEd structure.
parmed_structure = parmed.openmm.load_topology(top_mm, openmm_sys, pdbfile.positions)

# Export GROMACS files.
parmed_structure.save('BBR_OpenFF_pdb.top', overwrite=True)
parmed_structure.save('BBR_OpenFF_pdb.gro', overwrite=True)

#Validate Conversion
from simtk import openmm
for force in openmm_sys.getForces():
        if isinstance(force, openmm.NonbondedForce):
                    break
print(force.getCutoffDistance())
print(force.getUseSwitchingFunction())
print(force.getNonbondedMethod() == openmm.NonbondedForce.PME)
