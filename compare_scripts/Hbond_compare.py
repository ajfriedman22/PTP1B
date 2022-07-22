#!/ usr / bin / env python

#Import packages
from matplotlib_venn import venn2, venn3, venn2_circles
from matplotlib import pyplot as plt
import numpy as np

#Open files listing H_bonds in read-only format
list_a7 = open('../../rebuild_a7/analysis/Hbonds_a7.txt','r').readlines() #WPD open, a7 ordered, no AD
list_a7_AD = open('../../AD_rebuild_a7/analysis/Hbonds_a7_AD.txt', 'r').readlines() #WPD open, a7 ordered, AD present
list_a7_BBR = open('../../BBR_a7/analysis/Hbonds_BBR_a7.txt', 'r').readlines() #WPD open, a7 ordered, BBR present
list_Apo = open('../../Apo/analysis/Hbonds_Apo.txt','r').readlines() #WPD open, no a7, no AD
list_Apo_dis7 = open('../../rebuild_a7_high/config7/analysis/Hbonds_dis7.txt','r').readlines() #WPD open, no a7, no AD
list_Apo_dis9 = open('../../rebuild_a7_high/config9/analysis/Hbonds_dis9.txt','r').readlines() #WPD open, no a7, no AD
list_Apo_dis11 = open('../../rebuild_a7_high/config11/analysis/Hbonds_dis11.txt','r').readlines() #WPD open, no a7, no AD
list_Apo_dis_alt = open('../../Apo_dis/analysis/Hbonds_alt.txt','r').readlines() #WPD open, no a7, no AD
list_AD = open('../../AD/analysis/Hbonds_AD.txt','r').readlines() #WPD open, no a7, AD present
list_1sug = open('../../Apo_1SUG/analysis/1sug/Hbonds_1sug.txt','r').readlines() #WPD closed, a7 ordered, no AD
list_1sug2 = open('../../Apo_1SUG/analysis/1sug2/Hbonds_1sug2.txt','r').readlines() #WPD closed, a7 ordered, no AD
list_1sug3 = open('../../Apo_1SUG/analysis/1sug3/Hbonds_1sug3.txt','r').readlines() #WPD closed/open, a7 ordered, no AD
list_1sug_AD = open('../../1sug_AD/analysis/Hbonds_1sug_AD.txt','r').readlines() #WPD closed/open, a7 ordered, AD present
list_1sug_BBR = open('../../BBR_1sug/analysis/Hbonds_1sug_BBR.txt', 'r').readlines() #WPD closed, a7 ardered, BBR present
list_1sug_na7 = open('../../1sug_no_a7/analysis/Hbonds_1sug_naz7.txt','r').readlines() #WPD closed, no a7, no AD
list_1sug_na7_AD = open('../../1sug_no_a7_AD/analysis/Hbonds_1sug_na7_AD.txt','r').readlines() #WPD closed, no a7, no AD
list_1sug_dis = open('../../1sug_dis/analysis/config11/Hbonds_1sug_dis11.txt', 'r').readlines() #WPD closed, disordered a7, no AD
list_1sug_dis11_AD = open('../../1sug_dis_AD/analysis/config11/Hbonds_1sug_AD_dis11.txt', 'r').readlines() #WPD closed/open, disordered a7, AD present
list_1sug_dis11_2AD = open('../../1sug_dis_AD/analysis/config11_2/Hbonds_1sug_AD_dis11_2.txt', 'r').readlines() #WPD closed/open, disordered a7, AD present
list_1sug_alt_AD = open('../../1sug_dis_AD/analysis/config_alt/Hbonds_alt.txt', 'r').readlines() #WPD closed/open, disordered a7, AD present
list_1sug_alt2_AD = open('../../1sug_dis_AD/analysis/config_alt2/Hbonds_alt2.txt', 'r').readlines() #WPD closed/open, disordered a7, AD present

#Open file for atom numbers all bonds
list_atom_a7 = open('../../rebuild_a7/analysis/Hbonds_atom_a7.txt','r').readlines() #WPD open, a7 ordered, no AD
list_atom_a7_AD = open('../../AD_rebuild_a7/analysis/Hbonds_atom_a7_AD.txt', 'r').readlines() #WPD open, a7 ordered, AD present
list_atom_a7_BBR = open('../../BBR_a7/analysis/Hbonds_atom_BBR_a7.txt', 'r').readlines() #WPD open, a7 ordered, BBR present
list_atom_Apo = open('../../Apo/analysis/Hbonds_atom_Apo.txt','r').readlines() #WPD open, no a7, no AD
list_atom_Apo_dis7 = open('../../rebuild_a7_high/config7/analysis/Hbonds_atom_dis7.txt','r').readlines() #WPD open, no a7, no AD
list_atom_Apo_dis9 = open('../../rebuild_a7_high/config9/analysis/Hbonds_atom_dis9.txt','r').readlines() #WPD open, no a7, no AD
list_atom_Apo_dis11 = open('../../rebuild_a7_high/config11/analysis/Hbonds_atom_dis11.txt','r').readlines() #WPD open, no a7, no AD
list_atom_Apo_dis_alt = open('../../Apo_dis/analysis/Hbonds_atom_alt.txt','r').readlines() #WPD open, no a7, no AD
list_atom_AD = open('../../AD/analysis/Hbonds_atom_AD.txt','r').readlines() #WPD open, no a7, AD present
list_atom_1sug = open('../../Apo_1SUG/analysis/1sug/Hbonds_atom_1sug.txt','r').readlines() #WPD closed, a7 ordered, no AD
list_atom_1sug2 = open('../../Apo_1SUG/analysis/1sug2/Hbonds_atom_1sug2.txt','r').readlines() #WPD closed, a7 ordered, no AD
list_atom_1sug3 = open('../../Apo_1SUG/analysis/1sug3/Hbonds_atom_1sug3.txt','r').readlines() #WPD closed/open, a7 ordered, no AD
list_atom_1sug_AD = open('../../1sug_AD/analysis/Hbonds_atom_1sug_AD.txt','r').readlines() #WPD closed/open, a7 ordered, AD present
list_atom_1sug_BBR = open('../../BBR_1sug/analysis/Hbonds_atom_1sug_BBR.txt', 'r').readlines() #WPD closed, a7 ardered, BBR present
list_atom_1sug_na7 = open('../../1sug_no_a7/analysis/Hbonds_atom_1sug_naz7.txt','r').readlines() #WPD closed, no a7, no AD
list_atom_1sug_na7_AD = open('../../1sug_no_a7_AD/analysis/Hbonds_atom_1sug_na7_AD.txt','r').readlines() #WPD closed, no a7, no AD
list_atom_1sug_dis = open('../../1sug_dis/analysis/config11/Hbonds_atom_1sug_dis11.txt', 'r').readlines() #WPD closed, disordered a7, no AD
list_atom_1sug_dis11_AD = open('../../1sug_dis_AD/analysis/config11/Hbonds_atom_1sug_AD_dis11.txt', 'r').readlines() #WPD closed/open, disordered a7, AD present
list_atom_1sug_dis11_2AD = open('../../1sug_dis_AD/analysis/config11_2/Hbonds_atom_1sug_AD_dis11_2.txt', 'r').readlines() #WPD closed/open, disordered a7, AD present
list_atom_1sug_alt_AD = open('../../1sug_dis_AD/analysis/config_alt/Hbonds_atom_alt.txt', 'r').readlines() #WPD closed/open, disordered a7, AD present
list_atom_1sug_alt2_AD = open('../../1sug_dis_AD/analysis/config_alt2/Hbonds_atom_alt2.txt', 'r').readlines() #WPD closed/open, disordered a7, AD present

#Open output file or create if it doesn't exsist
output_uncomm = open('Hbond_uncommon.txt', 'w')
print('Input Files Open')

#Empty list for common hbonds
Apo_common, AD_common, Common_all, Apo_open, Apo_closed, AD_open, AD_closed, AD_open_crys= [],[],[],[],[],[],[],[]
AD_break, AD_form = [],[]

#Write all lines in common with all Apo files
for i in list_a7:
    if i in list_Apo and i in list_1sug and i in list_1sug_na7 and i in list_1sug_dis and i in list_1sug2 and i in list_1sug3 and i in list_Apo_dis7 and i in list_Apo_dis9 and i in list_Apo_dis11 and i in list_Apo_dis_alt:
        Apo_common.append(i)

#Write all lines in common with all AD files
for i in list_a7_AD:
    if i in list_AD and i in list_1sug_AD and i in list_1sug_na7_AD and i in list_1sug_dis11_AD and i in list_1sug_dis11_2AD and i in list_1sug_alt_AD and i in list_1sug_alt2_AD:
        AD_common.append(i)
#Write bonds that are common between all simms
for i in list_a7:
    if i in list_a7_AD and i in list_Apo and i in list_AD and i in list_1sug and i in list_1sug2 and i in list_1sug3 and i in list_1sug_AD and i in list_1sug_na7 and i in list_1sug_na7_AD and i in list_1sug_dis and i in list_1sug_dis11_AD and i in list_1sug_dis11_2AD and i in list_1sug_alt_AD and i in list_1sug_alt2_AD and i in list_a7_BBR and i in list_1sug_BBR and i in list_Apo_dis7 and i in list_Apo_dis9 and i in list_Apo_dis11 and i in list_Apo_dis_alt:
        Common_all.append(i)

#Write bonds in common betwen Apo open files
for i in list_a7:
    if i in list_Apo and i in list_1sug3 and i in list_Apo_dis7 and i in list_Apo_dis9 and i in list_Apo_dis11 and i in list_Apo_dis_alt:
        Apo_open.append(i)

#Write bonds in common between Apo closed files
for i in list_1sug:
    if i in list_1sug_dis and i in list_1sug2:
        Apo_closed.append(i)

#Write bonds in common b/w AD open files
for i in list_a7_AD:
    if i in list_AD and i in list_1sug_AD and i in list_1sug_dis11_AD and i in list_1sug_alt_AD and i in list_1sug_alt2_AD:
        AD_open.append(i)

#Write bonds in common b/w AD bound to crys + open
for i in list_1sug_dis11_AD:
    if i in list_1sug_alt_AD and i in list_1sug_dis11_2AD and i in list_1sug_alt2_AD:
        AD_open_crys.append(i)

#Write bonds in common b/w AD closed files
for i in list_1sug_na7_AD:
    AD_closed.append(i)

#Write all bonds broken by binding of AD
for i in Apo_common:
    if i not in AD_common:
        AD_break.append(i)

#Write all bonds formed by binding of AD
for i in AD_common:
    if i not in Apo_common:
        AD_form.append(i)

#Create VennDiagram b/w Apo and AD bound structures
plt.figure()
vd2=venn2([set(Apo_common),set(AD_common)], 
        set_labels = ('Apo Structures', 'AD Bound Structures'),
        set_colors = ('yellow', 'blue'))
venn2_circles([set(Apo_common), set(AD_common)], 
        linestyle='--', 
        linewidth=2, 
        color='black')
for text in vd2.set_labels:  #change label size
     text.set_fontsize(16);
     for text in vd2.subset_labels:  #change number size
          text.set_fontsize(16)
plt.title('Change in H-bond Network with AD Binding\n',
        fontsize=16,
        fontweight = 'bold')
plt.savefig('Venn_AD_binding.png')

#Create VennDiagram b/w open Apo and open AD bound structures
plt.figure()
vd2=venn2([set(Apo_open),set(AD_open_crys)], 
        set_labels = ('Open Apo Structures', 'Open AD Bound Crys Structures'),
        set_colors = ('yellow', 'blue'))
venn2_circles([set(Apo_open), set(AD_open_crys)], 
        linestyle='--', 
        linewidth=2, 
        color='black')
for text in vd2.set_labels:  #change label size
     text.set_fontsize(16);
     for text in vd2.subset_labels:  #change number size
          text.set_fontsize(16)
plt.title('Change in H-bond Network with AD Binding\n',
        fontsize=16,
        fontweight = 'bold')
plt.savefig('Venn_AD_open_binding.png')

#Create VennDiagram b/w open and closed AD bound structures
plt.figure()
vd2=venn2([set(AD_closed),set(AD_open)], 
        set_labels = ('Closed AD Bound', 'Open AD Bound'),
        set_colors = ('yellow', 'blue'))
venn2_circles([set(AD_closed), set(AD_open)], 
        linestyle='--', 
        linewidth=2, 
        color='black')
for text in vd2.set_labels:  #change label size
     text.set_fontsize(16);
     for text in vd2.subset_labels:  #change number size
          text.set_fontsize(16)
plt.title('Change in H-bond Network b/w AD opening WPD\n',
        fontsize=16,
        fontweight = 'bold')
plt.savefig('Venn_AD_open_close.png')

#Create VennDiagram b/w open and closed Apo structures
plt.figure()
vd2=venn2([set(Apo_closed),set(Apo_open)], 
        set_labels = ('Closed Apo', 'Open Apo'),
        set_colors = ('yellow', 'blue'))
venn2_circles([set(Apo_closed), set(Apo_open)], 
        linestyle='--', 
        linewidth=2, 
        color='black')
for text in vd2.set_labels:  #change label size
     text.set_fontsize(16);
     for text in vd2.subset_labels:  #change number size
          text.set_fontsize(16)
plt.title('Change in H-bond Network w/ WPD opening\n',
        fontsize=16,
        fontweight = 'bold')
plt.savefig('Venn_Apo_WPD.png')

print('Figures Completed')

#Determine Atom numbers on broken and formed bonds
File_atom_uncomm = open('Atom_uncomm.txt', 'w')
File_atom_uncomm_1sug = open('Atom_uncomm_1sug.txt', 'w')
#Determine list of uncommon bonds from all structures
uncomm = []
num = 0
for i in list_a7:
    if i not in Common_all:
        uncomm.append(i)
        File_atom_uncomm.write(list_atom_a7[num])
        atoms = list_atom_a7[num].split()
        File_atom_uncomm_1sug.write(str(float(atoms[0])-17) + ' ' + str(float(atoms[1])-17) + ' ' + str(float(atoms[2])-17) + '\n') #Change atom index for when the first residues is missing
        output_uncomm.write(i)
    num += 1
num = 0
for i in list_a7_AD:
    if i not in Common_all and i not in uncomm:
        uncomm.append(i)
        File_atom_uncomm.write(list_atom_a7_AD[num])
        atoms = list_atom_a7_AD[num].split()
        File_atom_uncomm_1sug.write(str(float(atoms[0])-17) + ' ' + str(float(atoms[1])-17) + ' ' + str(float(atoms[2])-17) + '\n')
        output_uncomm.write(i)
num = 0
for i in list_a7_BBR:
    if i not in Common_all and i not in uncomm:
        uncomm.append(i)
        File_atom_uncomm.write(list_atom_a7_BBR[num])
        atoms = list_atom_a7_BBR[num].split()
        File_atom_uncomm_1sug.write(str(float(atoms[0])-17) + ' ' + str(float(atoms[1])-17) + ' ' + str(float(atoms[2])-17) + '\n')
        output_uncomm.write(i)
num = 0
for i in list_Apo:
    if i not in Common_all and i not in uncomm:
        uncomm.append(i)
        File_atom_uncomm.write(list_atom_Apo[num])
        atoms = list_atom_Apo[num].split()
        File_atom_uncomm_1sug.write(str(float(atoms[0])-17) + ' ' + str(float(atoms[1])-17) + ' ' + str(float(atoms[2])-17) + '\n')
        output_uncomm.write(i)
    num += 1
num = 0
for i in list_Apo_dis7:
    if i not in Common_all and i not in uncomm:
        uncomm.append(i)
        File_atom_uncomm.write(list_atom_Apo_dis7[num])
        atoms = list_atom_Apo_dis7[num].split()
        File_atom_uncomm_1sug.write(str(float(atoms[0])-17) + ' ' + str(float(atoms[1])-17) + ' ' + str(float(atoms[2])-17) + '\n')
        output_uncomm.write(i)
    num += 1
num = 0
for i in list_Apo_dis9:
    if i not in Common_all and i not in uncomm:
        uncomm.append(i)
        File_atom_uncomm.write(list_atom_Apo_dis9[num])
        atoms = list_atom_Apo_dis9[num].split()
        File_atom_uncomm_1sug.write(str(float(atoms[0])-17) + ' ' + str(float(atoms[1])-17) + ' ' + str(float(atoms[2])-17) + '\n')
        output_uncomm.write(i)
    num += 1
num = 0
for i in list_Apo_dis11:
    if i not in Common_all and i not in uncomm:
        uncomm.append(i)
        File_atom_uncomm.write(list_atom_Apo_dis11[num])
        atoms = list_atom_Apo_dis11[num].split()
        File_atom_uncomm_1sug.write(str(float(atoms[0])-17) + ' ' + str(float(atoms[1])-17) + ' ' + str(float(atoms[2])-17) + '\n')
        output_uncomm.write(i)
    num += 1
num = 0
for i in list_Apo_dis_alt:
    if i not in Common_all and i not in uncomm:
        uncomm.append(i)
        File_atom_uncomm.write(list_atom_Apo_dis_alt[num])
        atoms = list_atom_Apo_dis_alt[num].split()
        File_atom_uncomm_1sug.write(str(float(atoms[0])-17) + ' ' + str(float(atoms[1])-17) + ' ' + str(float(atoms[2])-17) + '\n')
        output_uncomm.write(i)
    num += 1
num = 0
for i in list_AD:
    if i not in Common_all and i not in uncomm:
        uncomm.append(i)
        File_atom_uncomm.write(list_atom_AD[num])
        atoms = list_atom_AD[num].split()
        File_atom_uncomm_1sug.write(str(float(atoms[0])-17) + ' ' + str(float(atoms[1])-17) + ' ' + str(float(atoms[2])-17) + '\n')
        output_uncomm.write(i)
    num += 1
num = 0
for i in list_1sug:
    if i not in Common_all and i not in uncomm:
        uncomm.append(i)
        File_atom_uncomm_1sug.write(list_atom_1sug[num])
        atoms = list_atom_1sug[num].split()
        File_atom_uncomm.write(str(float(atoms[0])+17) + ' ' + str(float(atoms[1])+17) + ' ' + str(float(atoms[2])+17) + '\n')
        output_uncomm.write(i)
    num += 1
num = 0
for i in list_1sug2:
    if i not in Common_all and i not in uncomm:
        uncomm.append(i)
        File_atom_uncomm_1sug.write(list_atom_1sug2[num])
        atoms = list_atom_1sug2[num].split()
        File_atom_uncomm.write(str(float(atoms[0])+17) + ' ' + str(float(atoms[1])+17) + ' ' + str(float(atoms[2])+17) + '\n')
        output_uncomm.write(i)
    num += 1
num = 0
for i in list_1sug3:
    if i not in Common_all and i not in uncomm:
        uncomm.append(i)
        File_atom_uncomm_1sug.write(list_atom_1sug3[num])
        atoms = list_atom_1sug3[num].split()
        File_atom_uncomm.write(str(float(atoms[0])+17) + ' ' + str(float(atoms[1])+17) + ' ' + str(float(atoms[2])+17) + '\n')
        output_uncomm.write(i)
    num += 1
num = 0
for i in list_1sug_AD:
    if i not in Common_all and i not in uncomm:
        uncomm.append(i)
        File_atom_uncomm_1sug.write(list_atom_1sug_AD[num])
        atoms = list_atom_1sug_AD[num].split()
        File_atom_uncomm.write(str(float(atoms[0])+17) + ' ' + str(float(atoms[1])+17) + ' ' + str(float(atoms[2])+17) + '\n')
        output_uncomm.write(i)
    num += 1
num = 0
for i in list_1sug_na7:
    if i not in Common_all and i not in uncomm:
        uncomm.append(i)
        File_atom_uncomm_1sug.write(list_atom_1sug_na7[num])
        atoms = list_atom_1sug_na7[num].split()
        File_atom_uncomm.write(str(float(atoms[0])+17) + ' ' + str(float(atoms[1])+17) + ' ' + str(float(atoms[2])+17) + '\n')
        output_uncomm.write(i)
    num += 1
num = 0
for i in list_1sug_na7_AD:
    if i not in Common_all and i not in uncomm:
        uncomm.append(i)
        File_atom_uncomm_1sug.write(list_atom_1sug_na7_AD[num])
        atoms = list_atom_1sug_na7_AD[num].split()
        File_atom_uncomm.write(str(float(atoms[0])+17) + ' ' + str(float(atoms[1])+17) + ' ' + str(float(atoms[2])+17) + '\n')
        output_uncomm.write(i)
    num += 1
num = 0
for i in list_1sug_dis:
    if i not in Common_all and i not in uncomm:
        uncomm.append(i)
        File_atom_uncomm_1sug.write(list_atom_1sug_dis[num])
        atoms = list_atom_1sug_dis[num].split()
        File_atom_uncomm.write(str(float(atoms[0])+17) + ' ' + str(float(atoms[1])+17) + ' ' + str(float(atoms[2])+17) + '\n')
        output_uncomm.write(i)
    num += 1
num = 0
for i in list_1sug_dis11_AD:
    if i not in Common_all and i not in uncomm:
        uncomm.append(i)
        File_atom_uncomm_1sug.write(list_atom_1sug_dis11_AD[num])
        atoms = list_atom_1sug_dis11_AD[num].split()
        File_atom_uncomm.write(str(float(atoms[0])+17) + ' ' + str(float(atoms[1])+17) + ' ' + str(float(atoms[2])+17) + '\n')
        output_uncomm.write(i)
    num += 1
num = 0
for i in list_1sug_dis11_2AD:
    if i not in Common_all and i not in uncomm:
        uncomm.append(i)
        File_atom_uncomm_1sug.write(list_atom_1sug_dis11_2AD[num])
        atoms = list_atom_1sug_dis11_2AD[num].split()
        File_atom_uncomm.write(str(float(atoms[0])+17) + ' ' + str(float(atoms[1])+17) + ' ' + str(float(atoms[2])+17) + '\n')
        output_uncomm.write(i)
    num += 1
num = 0
for i in list_1sug_alt_AD:
    if i not in Common_all and i not in uncomm:
        uncomm.append(i)
        File_atom_uncomm_1sug.write(list_atom_1sug_alt_AD[num])
        atoms = list_atom_1sug_alt_AD[num].split()
        File_atom_uncomm.write(str(float(atoms[0])+17) + ' ' + str(float(atoms[1])+17) + ' ' + str(float(atoms[2])+17) + '\n')
        output_uncomm.write(i)
    num += 1
num = 0
for i in list_1sug_alt2_AD:
    if i not in Common_all and i not in uncomm:
        uncomm.append(i)
        File_atom_uncomm_1sug.write(list_atom_1sug_alt2_AD[num])
        atoms = list_atom_1sug_alt2_AD[num].split()
        File_atom_uncomm.write(str(float(atoms[0])+17) + ' ' + str(float(atoms[1])+17) + ' ' + str(float(atoms[2])+17) + '\n')
        output_uncomm.write(i)
    num += 1
num = 0
for i in list_1sug_BBR:
    if i not in Common_all and i not in uncomm:
        uncomm.append(i)
        File_atom_uncomm_1sug.write(list_atom_1sug_BBR[num])
        atoms = list_atom_1sug_BBR[num].split()
        File_atom_uncomm.write(str(float(atoms[0])+17) + ' ' + str(float(atoms[1])+17) + ' ' + str(float(atoms[2])+17) + '\n')
        output_uncomm.write(i)
    num += 1
