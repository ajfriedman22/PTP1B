from matplotlib import pyplot as plt
from matplotlib_venn import venn2, venn2_circles
from matplotlib_venn import venn3, venn3_circles

#AD contacts 
AD = ['MET282', 'PRO241', 'ALA278', 'PHE196', 'MET235', 'MET1', 'ILE281', 'PHE280', 'LEU232', 'LEU192', 'GLY277', 'ASP236'] #65% alt2 25% alt3
AD_a7 = ['MET3', 'LYS239', 'MET235', 'PHE280', 'ILE281', 'VAL274', 'LEU232', 'ALA278', 'ASP236', 'PRO241', 'GLY277'] #33% alt2 65% alt3
AD_dis7 = ['LEU294', 'ARG45', 'MET133', 'ASN42', 'LYS120', 'MET114', 'SER151', 'ALA122', 'PRO89', 'SER190', 'LEU88', 'GLN123', 'ASN90', 'PHE135', 
        'LEU119', 'HIS296', 'TYR152', 'GLN78', 'GLU186', 'CYS121', 'ARG43'] #Unstable bonding
AD_dis9 = ['ASP240', 'ILE281', 'ARG105', 'MET235', 'PRO241', 'ASP236', 'MET282', 'LYS103', 'LYS239', 'ALA278', 'ARG169'] #86% alt2
AD_dis11 = ['PHE196', 'SER190', 'LYS279', 'VAL287', 'TRP291', 'ASP289', 'SER187', 'GLN288', 'SER295', 'GLY283', 'ALA189', 'TYR152', 'PHE280', 'GLN290'] #83% alt1
AD_1sug = ['ASP236', 'GLY277', 'PHE280', 'ALA278', 'LYS239', 'MET3', 'MET235', 'PRO241', 'ILE281', 'LEU232'] #67% alt2 33% alt3
AD_1sug_na7 = ['LEU232', 'ASP236', 'PHE196', 'LEU192', 'ALA189', 'LEU195', 'PHE280', 'ARG199', 'VAL287', 'ILE281', 'GLY277', 'ASN193'] #86% alt3
AD_1sug_dis7 = ['PHE135', 'ARG268', 'PRO188', 'LEU272', 'PHE269', 'MET133', 'PRO89', 'ALA189', 'ASN42', 'GLN123', 'GLU276', 'LEU192', 'CYS92', 'ALA273'] #27% Crys binding loc
AD_1sug_dis9 = []
AD_1sug_dis11 = ['LEU232', 'ASP289', 'PHE280', 'PHE196', 'GLY277', 'TRP291', 'ASN193', 'ASP298', 'LEU192', 'GLU297', 'SER286', 'ALA189', 'GLN290'] #95% crys
#AD_1sug_dis11_2 =
AD_1sug_alt = ['PHE280', 'GLU297', 'TYR152', 'PHE196', 'LYS292', 'HIS296', 'GLN288', 'SER295', 'TRP291', 'ASP289', 'ALA189', 'LYS197', 'ASN193'] #100% alt1
#AD_1sug_alt2 =

#BBR contacts
bbr_a7 = ['SER295', 'PHE280', 'GLN288', 'LYS279', 'VAL287', 'LYS292', 'TRP291', 'ALA189', 'GLU293', 'PHE196', 'SER285', 'GLU276', 'ASN193', 'LEU192',
        'GLN290', 'LEU294']
bbr_1sug = ['PRO188', 'LEU294', 'GLU186', 'GLN290', 'ALA189', 'PHE280', 'GLU276', 'GLU293', 'SER295', 'SER187', 'VAL287', 'SER286', 'LYS279', 'TRP291']
bbr_dis9 = ['GLN288', 'TRP291', 'GLU200', 'PHE196', 'LYS292', 'GLN290', 'ARG199', 'PHE280', 'LYS279', 'LEU294', 'ASN193', 'ILE281', 'ASP289', 'GLU293']
bbr_1sug_dis7 = ['LYS292', 'GLN290', 'ARG268', 'VAL287', 'LEU272', 'SER187', 'PHE280', 'GLU293', 'ILE275', 'LYS279', 'GLN288', 'PHE7', 'TRP291',
        'ALA189', 'GLU276', 'PRO188']

#Output folder for unique bonds to open conformations
file_uniq = open('AD_interact.txt', 'w')
file_uniq_bbr = open('BBR_open.txt', 'w')
file_AD_BBR = open('AD_BBR_common.txt', 'w')
file_AD_no_BBR = open('AD_BBR_uncommon.txt', 'w')
crys, alt1, alt2, alt3 = [],[],[],[]
open_bbr_all, close_bbr_all, open_bbr_uniq = [],[],[]

#AD crystal binding site
for i in AD_1sug_dis11:
    crys.append(i)

#AD alt 1
file_uniq.write('Alt 1' + '\n')
for i in AD_1sug_alt:
    alt1.append(i)
    if i not in AD_dis11:
        file_uniq.write(i + '\n')

#AD alt 2
file_uniq.write('Alt 2' + '\n')
for i in AD_a7:
    if i in AD_1sug:
        alt2.append(i)
        if i not in AD_dis9:
            file_uniq.write(i + '\n')
            alt2.append(i)

#AD alt 3
file_uniq.write('Alt 3' + '\n')
for i in AD:
    alt3.append(i)
    if i not in AD_1sug_na7:
        file_uniq.write(i + '\n')

#Compare crystal to alt 1
plt.figure()
vd2=venn2([set(crys),set(alt1)], 
        set_labels = ('Binding Loc 1', 'Binding Loc 2'))
venn2_circles([set(crys), set(alt1)],
        linestyle='--', 
        linewidth=2, 
        color='black')
for text in vd2.set_labels:  #change label size
     text.set_fontsize(16);
plt.title('Change AD Interactions by binding site\n',
        fontsize=16,
        fontweight = 'bold')
plt.savefig('Venn_AD_crys_alt1_bind.png')

#Compare crystal to alt 2/3
plt.figure()
vd3=venn3([set(crys),set(alt2), set(alt3)], 
        set_labels = ('Binding Loc 1', 'Binding Loc 3', 'Binding Loc 4'))
venn3_circles([set(crys), set(alt2), set(alt3)],
        linestyle='--', 
        linewidth=2, 
        color='black')
for text in vd3.set_labels:  #change label size
     text.set_fontsize(16);
plt.title('Change AD Interactions by binding site\n',
        fontsize=16,
        fontweight = 'bold')
plt.savefig('Venn_AD_crys_alt2-3_bind.png')

#BBR open
for i in bbr_a7:
    open_bbr_all.append(i)

for i in bbr_dis9:
    if i not in open_bbr_all:
        open_bbr_all.append(i)

#BBR closed
for i in bbr_1sug:
    close_bbr_all.append(i)

for i in bbr_1sug_dis7:
    if i not in close_bbr_all:
        close_bbr_all.append(i)

#Compare BBR open and closed
plt.figure()
vd2=venn2([set(open_bbr_all),set(close_bbr_all)], 
        set_labels = ('WPD Open', 'WPD Closed'))
venn2_circles([set(open_bbr_all), set(close_bbr_all)],
        linestyle='--', 
        linewidth=2, 
        color='black')
for text in vd2.set_labels:  #change label size
     text.set_fontsize(16);
plt.title('Change BBR Interactions by Open or Closed Bound States\n',
        fontsize=16,
        fontweight = 'bold')
plt.savefig('Venn_BBR_bind.png')

#Print bonds unique to the open conformation
for i in open_bbr_all:
    if i not in close_bbr_all:
        file_uniq_bbr.write(i + '\n')
        open_bbr_uniq.append(i)

#Compare BBR to AD binding
file_AD_BBR.write('Crystal Binding Site\n')
for i in crys:
    if i in open_bbr_uniq:
        file_AD_BBR.write(i + '\n')
file_AD_BBR.write('Alt1 Binding Site\n')
for i in alt1:
    if i in open_bbr_uniq:
        file_AD_BBR.write(i + '\n')

file_AD_no_BBR.write('BBR Binding Site\n')
for i in open_bbr_uniq:
    if i not in crys and i not in alt1:
        file_AD_no_BBR.write(i + '\n')

file_AD_no_BBR.write('Crystal Binding Site\n')
for i in crys:
    if i not in open_bbr_all:
        file_AD_no_BBR.write(i + '\n')
file_AD_no_BBR.write('Alt1 Binding Site\n')
for i in alt1:
    if i not in open_bbr_all:
        file_AD_no_BBR.write(i + '\n')
plt.figure()
vd3=venn3([set(crys),set(alt1), set(open_bbr_all)], 
        set_labels = ('AD Loc 1', 'AD Loc 2', 'BBR'))
venn3_circles([set(crys), set(alt1), set(open_bbr_all)],
        linestyle='--', 
        linewidth=2, 
        color='black')
for text in vd3.set_labels:  #change label size
     text.set_fontsize(16);
plt.title('Difference in Interactions b/w BBR and AD\n',
        fontsize=16,
        fontweight = 'bold')
plt.savefig('Venn_AD_BBR_bind.png')

