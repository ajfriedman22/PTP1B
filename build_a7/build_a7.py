#!/usr/bin/env python

"""
Build the missing residues based on input structure (input.pdb) and template structure (temp.pdb)
"""

import argparse
from modeller import *
from modeller.automodel import *

#Create the parser
parser = argparse.ArgumentParser()

#Input argument for PDB file
parser.add_argument('-i', type=str, required=True)
parser.add_argument('-r', type=str, required=True)
parser.add_argument('-a', type=str, required=True)

#Parse the argument
args = parser.parse_args()

#Set input and template files
i = args.i
r = args.r
align = args.a

log.verbose()
env = Environ()

#Directories for input atom files
env.io.atom_files_directory = ['.', '../atom_files']

#Model in only missing residues from a7 helix
class MyModel(AutoModel):
	def special_restraints(self, aln):
		rsr = self.restraints
		at = self.atoms
		#Residue 288-299 should be an alpha helix
		rsr.add(secondary_structure.Alpha(self.residue_range('288:A', '299:A')))
		return Selection(self.residue_range('288:A', '297:A'))

a = MyModel(env, alnfile = align, knowns = i, sequence = r)
a.starting_model=1
a.ending_model=1

a.make()

