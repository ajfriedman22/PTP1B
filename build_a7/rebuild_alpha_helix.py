#!/usr/bin/env python

"""
Build the missing residues based on input structure (input.pdb) and template structure (temp.pdb)
"""

import argparse
from modeller import *
from modeller.optimizers import ConjugateGradients

log.verbose()
env = Environ()
env.libs.topology.read('${LIB}/top_heav.lib')
env.libs.parameters.read('${LIB}/par.lib')
env.io.atom_files_directory = ['./']

code = 'PTP1B_fill'
##Read full input atom file
mdl = Model(env)
aln = Alignment(env)
aln.append(file='PTP1B_fill.ali', align_codes=code)

#Build an extended chain model
mdl.read(file=code)
mdl.generate_topology(aln[code])

#Make steriochemistry restraints
allatoms = Selection(mdl)
mdl.restraints.make(allatoms, restraint_type='STEREO', spline_on_site=False)

#Constrain residues 288-299
mdl.restraints.add(secondary_structure.Alpha(mdl.residue_range('288:A', '299:A')))

#Get an optimized structure with CG
cg = ConjugateGradients()
cg.optimize(mdl.residue_range('288:A', '299:A'), max_iterations=100)
m.write(file='PTP1B_a7_helix.pdb')

