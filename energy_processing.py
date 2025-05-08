import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as stats #statistical tests

obl_list=pd.read_csv('obl_pdbs',header=None)
non_obl_list=pd.read_csv('non_obl_pdbs',header=None)

#Directly PPCheck run
#=============
#OBLIGATORY
#=============

obl_pdb=[]
obl_energy=[]#total energy of pdbs
obl_interface_res_count=[]#number of interface residues
os.chdir('/home/sidhant/Sidhanta/work/multidomain/energy/PPCheck/obl_results/')
for i in range(len(obl_list)):
    pdb_dir=obl_list.iloc[i,0]
    #print(pdb_dir)
    os.chdir(pdb_dir)
    obl_pdb.append(pdb_dir)
    obl_energy.append(os.popen("grep ^Total A-B_results.dat | grep -oP '\-*\d+\.\d+'").read().split('\n')[0])
    obl_interface_res_count.append(os.popen("grep ^Number A-B_results.dat | grep -oP '\d+'").read().split('\n')[0])
    os.chdir('..')

#==============
#NON-OBLIGATORY
#==============
#modified on 22/4/22 to not include 5IK4 as it has positive energy even after energy minimization.

non_obl_pdb=[]
non_obl_energy=[] # non-obl energies
non_obl_interface_res_count=[] #non-obl interface res count
os.chdir("/home/sidhant/Sidhanta/work/multidomain/energy/PPCheck/non_obl_results/")
for j in range(len(non_obl_list)):
    Npdb_dir=non_obl_list.iloc[j,0]
    #print(Npdb_dir)
    if Npdb_dir=='5ik4':
        continue
    else:
        os.chdir(Npdb_dir)
        non_obl_pdb.append(Npdb_dir)
        non_obl_energy.append(os.popen("grep ^Total A-B_results.dat | grep -oP '\-*\d+\.\d+'").read().split('\n')[0])
        non_obl_interface_res_count.append(os.popen("grep ^Number A-B_results.dat | grep -oP '\d+'").read().split('\n')[0])
        os.chdir('..')

os.chdir("/home/sidhant/Sidhanta/work/multidomain/energy/PPCheck/random/")

obl_energy_float=[float(k) for k in obl_energy]
#obl_energy_float[57]= -152.122 #editing the energy of 1zrl manually. refer nb 22/4/22. replaced energy min with original.
obl_interface_res_count_float=[float(k) for k in obl_interface_res_count]

non_obl_energy_float=[float(k) for k in non_obl_energy]
non_obl_e_rand=[non_obl_energy_float[i] for i in np.random.randint(0,153,size=109)]#randomizing

#non_obl_energy_float[124]= -234.3631 #editing the energy of 3wvt manually. refer nb 22/4/22. file replaced
non_obl_interface_res_count_float=[float(k) for k in non_obl_interface_res_count]

obl_results=pd.DataFrame(list(zip(obl_pdb,obl_energy_float,obl_interface_res_count_float)),columns=['PDB','Energy','Interface_res_count'])
non_obl_results=pd.DataFrame(list(zip(non_obl_pdb,non_obl_energy_float,non_obl_interface_res_count_float)),columns=['PDB','Energy','Interface_res_count'])

only_interface_count=list((obl_results.Interface_res_count,non_obl_results.Interface_res_count))# interface res count


#not needed as foldx repair gives little higher energy. see plots.
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#After_repair by foldx RepairPDB
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#=============
#OBLIGATORY
#=============
"""
repair_obl_pdb=[]
repair_obl_energy=[]#total energy of pdbs
repair_obl_interface_res_count=[]#number of interface residues
os.chdir('../after_repair/obl_results/')
for i in range(len(obl_list)):
    pdb_dir=obl_list.iloc[i,0]
    #print(pdb_dir)
    os.chdir(pdb_dir)
    repair_obl_pdb.append(pdb_dir)
    repair_obl_energy.append(os.popen("grep ^Total A-B_results.dat | grep -oP '\-*\d+\.\d+'").read().split('\n')[0])
    repair_obl_interface_res_count.append(os.popen("grep ^Number A-B_results.dat | grep -oP '\d+'").read().split('\n')[0])
    os.chdir('..')

#==============
#NON-OBLIGATORY
#==============

repair_non_obl_pdb=[]
repair_non_obl_energy=[] # non-obl energies
repair_non_obl_interface_res_count=[] #non-obl interface res count
os.chdir("../non_obl_results/")
for j in range(len(non_obl_list)):
    Npdb_dir=non_obl_list.iloc[j,0]
    #print(Npdb_dir)
    os.chdir(Npdb_dir)
    repair_non_obl_pdb.append(Npdb_dir)
    repair_non_obl_energy.append(os.popen("grep ^Total A-B_results.dat | grep -oP '\-*\d+\.\d+'").read().split('\n')[0])
    repair_non_obl_interface_res_count.append(os.popen("grep ^Number A-B_results.dat | grep -oP '\d+'").read().split('\n')[0])
    os.chdir('..')


repair_obl_energy_float=[float(k) for k in repair_obl_energy]
repair_obl_interface_res_count_float=[float(k) for k in repair_obl_interface_res_count]

repair_non_obl_energy_float=[float(k) for k in repair_non_obl_energy]
repair_non_obl_interface_res_count_float=[float(k) for k in repair_non_obl_interface_res_count]

repair_obl_results=pd.DataFrame(list(zip(repair_obl_pdb,repair_obl_energy_float,repair_obl_interface_res_count_float)),columns=['PDB','Energy','Interface_res_count'])
repair_non_obl_results=pd.DataFrame(list(zip(repair_non_obl_pdb,repair_non_obl_energy_float,repair_non_obl_interface_res_count_float)),columns=['PDB','Energy','Interface_res_count'])
"""


print('p-value: ',stats.mannwhitneyu(obl_energy_float,non_obl_e_rand)[1])
print('ks-test: ',stats.ks_2samp(obl_energy_float,non_obl_e_rand))
#------------------------------------------
#plotting
#------------------------------------------
#energy plot
energies=list((obl_energy_float,non_obl_e_rand))#obl and non-obl energies in a list
ax=sns.violinplot(data=energies)
ax.set_xticklabels(['Obligatory','Non-obligatory'],weight='bold')
plt.ylabel('Energy in KJ/mol',weight='bold')
#plt.savefig('energy_rand5.png',format='png',dpi=1300)
plt.show()
#only energy is randomized

"""
#swarmplot of number of interface residues
ax1=sns.swarmplot(data=only_interface_count)
ax1.set_xticklabels(['Obligatory','Non-obligatory'],weight='bold')
plt.ylabel('Number of interfacial residues',fontweight='bold')
plt.yticks(np.arange(0,336,25))
plt.show()

repair_energies=list((repair_obl_energy_float,repair_non_obl_energy_float))
ax1=sns.violinplot(data=repair_energies)
ax1.set_xticklabels(['Obligatory','Non-obligatory'],weight='bold')
plt.ylabel('Energy in KJ/mol',weight='bold')
plt.title('after repair')
plt.show()
"""
