#!/usr/bin/env python

"""
Get Alignment file from two PDB files
"""

from modeller import *
from modeller.automodel import *

#input pdb name
code = "PTP1B_a7_helix"
code_ref = "1sug_clean"

#Set MODELLER Environment
env = Environ()
env.io.atom_file_directory = "../atom_files/"

#Create Empty Alignment and model
aln = Alignment(env)

#Read full input atom file
mdl = Model(env, file=code)

#Add the model sequence to the alignment
aln.append_model(mdl, atom_files=code, align_codes=code)

#Read reference atom file
mdl_ref = Model(env, file=code_ref)
aln.append_model(mdl_ref, align_codes=code_ref, atom_files=code_ref)

#Allign them by structure
aln.salign(gap_penalties_3d=(0.0, 2.0))

#Check alignment for suitability for modeling
aln.check()

aln.write(file="PTP1B-a7-align.ali")

#Model in only missing residues from a7 helix
class MyModel(AutoModel):
	def select_atoms(self):
		return Selection(self.residue_range('280:A', '298:A'))

a = MyModel(env, alnfile = "PTP1B-a7-align.ali", knowns = code_ref, sequence = code)
a.starting_model=1
a.ending_model=5

a.make()


