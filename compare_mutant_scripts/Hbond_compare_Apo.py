#Import Necessary Packages
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as lines
from scipy import stats
from itertools import product
import seaborn as sns

#Read names of H-bonds
WT_name = open('../../../Apo_dis/analysis/Hbonds_alt_Apo.txt', 'r').readlines()
F196A_name = open('../../F196A/Apo/analysis/Hbonds_F196A_Apo.txt', 'r').readlines()
L192F_name = open('../../L192F/Apo/analysis/Hbonds_L192F_Apo.txt', 'r').readlines()
L195F_name = open('../../L195F/Apo/analysis/Hbonds_L195F_Apo.txt', 'r').readlines()
F280Y_name = open('../../F280Y/Apo/analysis/Hbonds_F280Y_Apo.txt', 'r').readlines()
E276F_name = open('../../E276F/Apo/analysis/Hbonds_E276F_Apo.txt', 'r').readlines()
V287T_name = open('../../V287T/Apo/analysis/Hbonds_V287T_Apo.txt', 'r').readlines()

#Read names of H-bonds
WT_atom = open('../../../Apo_dis/analysis/Hbonds_atom_alt_Apo.txt', 'r').readlines()
F196A_atom = open('../../F196A/Apo/analysis/Hbonds_atom_F196A_Apo.txt', 'r').readlines()
L192F_atom = open('../../L192F/Apo/analysis/Hbonds_atom_L192F_Apo.txt', 'r').readlines()
L195F_atom = open('../../L195F/Apo/analysis/Hbonds_atom_L195F_Apo.txt', 'r').readlines()
F280Y_atom = open('../../F280Y/Apo/analysis/Hbonds_atom_F280Y_Apo.txt', 'r').readlines()
E276F_atom = open('../../E276F/Apo/analysis/Hbonds_atom_E276F_Apo.txt', 'r').readlines()
V287T_atom = open('../../V287T/Apo/analysis/Hbonds_atom_V287T_Apo.txt', 'r').readlines()

#Determine bonds common to all
file_all_bonds = open('Hbonds_common_all.txt', 'w')
common_bonds = []
for i in WT_name:
    if i in F196A_name and i in L192F_name and i in L195F_name and i in F280Y_name and i in E276F_name and i in V287T_name:
        file_all_bonds.write(str(i))
        common_bonds.append(i)

#Bonds common to those mutants with activity equal to or less than but not those greater than WT activity
file_allo_network = open('hbonds_possible_allo_network.txt', 'w')
file_allo_network_atom = open('hbonds_possible_allo_network_atom.txt', 'w')
n = 0
for i in WT_name:
    if i in L192F_name and i in F280Y_name and i not in F196A_name and i not in L195F_name and i not in E276F and i not in V287T_name:
        file_allo_network.write(str(i))
        file_allo_network_atom.write(str(WT_atom[n]) + '\n')
    n += 1

n = 0
file_allo_network_strict = open('Hbonds_possible_allo_network_strict.txt', 'w')
file_allo_network_strict_atom = open('Hbonds_possible_allo_network_strict_atom.txt', 'w')
for i in F280Y_name:
    if i not in F196A_name and i not in L195F_name and i not in V287T_name:
        file_allo_network_strict.write(str(i))
        file_allo_network_strict_atom.write(str(F280Y_atom[n]) + '\n')
    n += 1

#Bonds common to those mutants with activity greater than WT but not those less than or equal to WT activity
file_allo_break_network = open('Hbonds_possible_allo_break_network.txt', 'w')
file_allo_break_network_atom = open('Hbonds_possible_allo_break_network_atom.txt', 'w')
n = 0
for i in V287T_name:
    if i not in L192F_name and i not in F280Y_name and i not in WT_name and i in F196A_name and i in L195F_name and i in E276F:
        file_allo_break_network.write(str(i))
        file_allo_break_network_atom.write(str(V287T_atom[n]) + '\n')
    n += 1

file_allo_break_network_strict = open('Hbonds_possible_allo_break_network_strict.txt', 'w')
file_allo_break_network_strict_atom = open('Hbonds_possible_allo_break_network_strict_atom.txt', 'w')
n = 0
for i in V287T_name:
    if i not in F280Y_name and i in F196A_name and i in L195F_name:
        file_allo_break_network_strict.write(str(i))
        file_allo_break_network_strict_atom.write(str(V287T_atom[n]) + '\n')
    n += 1

#Determine all bonds which are present in at least one trajectory but not in all of them
uncommon_bonds = []
file_uncommon = open('Hbonds_uncommon_Apo.txt', 'w')
file_uncommon_atom = open('Hbonds_uncommon_atom_Apo.txt', 'w')
n = 0
for i in WT_name:
    if i not in common_bonds:
        uncommon_bonds.append(i)
        file_uncommon.write(i)
        file_uncommon_atom.write(WT_atom[n])
    n += 1
n = 0
for i in F196A_name:
    if i not in common_bonds and i not in uncommon_bonds:
        uncommon_bonds.append(i)
        file_uncommon.write(i)
        file_uncommon_atom.write(F196A_atom[n])
    n += 1
n = 0
for i in L192F_name:
    if i not in common_bonds and i not in uncommon_bonds:
        uncommon_bonds.append(i)
        file_uncommon.write(i)
        file_uncommon_atom.write(L192F_atom[n])
    n += 1
n = 0
for i in L195F_name:
    if i not in common_bonds and i not in uncommon_bonds:
        uncommon_bonds.append(i)
        file_uncommon.write(i)
        file_uncommon_atom.write(L195F_atom[n])
    n += 1
n = 0
for i in E276F_name:
    if i not in common_bonds and i not in uncommon_bonds:
        uncommon_bonds.append(i)
        file_uncommon.write(i)
        file_uncommon_atom.write(E276F_atom[n])
    n += 1
n = 0
for i in F280Y_name:
    if i not in common_bonds and i not in uncommon_bonds:
        uncommon_bonds.append(i)
        file_uncommon.write(i)
        file_uncommon_atom.write(F280Y_atom[n])
    n += 1
n = 0
for i in V287T_name:
    if i not in common_bonds and i not in uncommon_bonds:
        uncommon_bonds.append(i)
        file_uncommon.write(i)
        file_uncommon_atom.write(V287T_atom[n])
    n += 1

