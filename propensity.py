#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 27 13:06:26 2023

@author: sidhant
"""
import os
#from interface_residues_pymol import interfaceResidues
import pandas as pd
#from pymol import cmd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from Bio.PDB import *
#from scipy.stats import mannwhitneyu #mann-whitney test
import scipy.stats as stats #ks 2 sample test
import collections


parent_dir= '/home/sidhant/Sidhanta/work/multidomain/energy/PPCheck/'

data=pd.read_csv('obl_pdbs', header =None)
data1=pd.read_csv('non_obl_pdbs', header= None)

#`````````````````````````````````````````````````````````````````````
def resnames(intrf):#function to find the residue names from residue numbers using bioPDB
    temp=[]
    for z in range(len(chains)):
        for x in chains[z].get_residues():
            residue_id=x.get_id()
            if residue_id[0]==' ':
                res_num=residue_id[1]
                if res_num in intrf:
                    temp.append(x.get_resname())
    return temp                
#```````````````````````````````````````````````````````````````````````

#OBLIGATORY--------interface_residues--------------------------------------               
o_intf_res_numbers=[]#list of interface residue numbers
o_intf_dict=dict()#pdb_id and interfacial residue names in dict
i=0
for pdb in data.iloc[:,0]:
    #os.chdir(parent_dir+'obl_pdb_files/'+pdb)
    pdb_file=parent_dir+'obl_pdb_files/'+pdb+'_edit.pdb'
    cmd.load(pdb_file,pdb)
    cmd.disable('all')
    cmd.enable(pdb)
    #cmd.get_chains()
    a=interfaceResidues(pdb,cA="c. A", cB="c. B", cutoff=1, selName="foundIntrf")
    #interface_area_obl.append(float(a))
    #obl_pdb.append(pdb)
    i+=1
    print(i,'done1',pdb)
    #print([int(item[1]) for item in a])#interfacial residue numbers
    o_intf_res_numbers.append([int(item[1]) for item in a])
    
    #------ BioPDB
    parser=PDBParser(QUIET=True)
    structure=parser.get_structure(pdb, pdb_file)
    chains=list(structure.get_chains())
    lis=resnames([int(item[1]) for item in a])
    #print(lis)
    o_intf_dict[pdb]=lis
    #print(o_intf_dict)
    
    
    
#NON-OBLIGATORY-----------interface_residues------------------------------------------
n_intf_res_numbers=[]#list of interface residue numbers
n_intf_dict=dict()#pdb_id and interfacial residue names in dict
j=0
for pdb in data1.iloc[:,0]:
    #os.chdir(parent_dir+'obl_pdb_files/'+pdb)
    pdb_file=parent_dir+'non_obl_pdb_files/'+pdb+'_edit.pdb'
    cmd.load(pdb_file,pdb)
    cmd.disable('all')
    cmd.enable(pdb)
    #cmd.get_chains()
    a=interfaceResidues(pdb,cA="c. A", cB="c. B", cutoff=1, selName="foundIntrf")
    #interface_area_obl.append(float(a))
    #obl_pdb.append(pdb)
    j+=1
    print(j,'done2',pdb)
    #print([int(item[1]) for item in a])#interfacial residue numbers
    n_intf_res_numbers.append([int(item[1]) for item in a])
    
    #------ BioPDB
    parser=PDBParser(QUIET=True)
    structure=parser.get_structure(pdb, pdb_file)
    chains=list(structure.get_chains())
    lis=resnames([int(item[1]) for item in a])
    #print(lis)
    n_intf_dict[pdb]=lis
    #print(n_intf_dict)




o_intf_residues=sum(list(o_intf_dict.values()),[])#making list of lists to list
n_intf_residues=sum(list(n_intf_dict.values()),[])

#collctions.Counter(o_intf_residues)#gives 21 a.a. types. 2GNX has 2 'UNK' residues. 
o_intf_residues.remove('UNK')
o_intf_residues.remove('UNK')#as this list has two 'UNK' residues.

residues=sorted(collections.Counter(o_intf_residues).keys())
o_res_freq=collections.Counter(o_intf_residues)#frequency
n_res_freq=collections.Counter(n_intf_residues)


#for obligatory
o_propensity=[]
for i in residues:
    num1=o_res_freq[i]
    num2=o_res_freq[i]+n_res_freq[i]
    num3=len(o_intf_residues)
    num4=len(o_intf_residues)+len(n_intf_residues)
    propensity=(num1/num2)/(num3/num4)
    o_propensity.append("%.2f" %propensity)
    
#for non-obligatory
n_propensity=[]
for i in residues:
    num1=n_res_freq[i]
    num2=o_res_freq[i]+n_res_freq[i]
    num3=len(n_intf_residues)
    num4=len(o_intf_residues)+len(n_intf_residues)
    propensity=(num1/num2)/(num3/num4)
    n_propensity.append("%.2f" %propensity)
    
propensity_df=pd.DataFrame({'Residues':residues,'O-propensity':o_propensity,'N-propensity':n_propensity})
propensity_df = propensity_df.astype({"O-propensity":float, "N-propensity": float})
print(propensity_df)

"""
#plotting!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
propensity_df['O-propensity']=pd.to_numeric(propensity_df['O-propensity'])
propensity_df['N-propensity']=pd.to_numeric(propensity_df['N-propensity'])
sns.set_style("darkgrid")
plt.figure(figsize=(14.0, 8.0))
x_axis=np.arange(0,20)
#obl_std=[np.std(df_obl.iloc[:,i]) for i in range(1,10)]
#non_obl_std=[np.std(df_non_obl.iloc[:,i]) for i in range(1,10)]
plt.bar(x_axis-0.2,propensity_df['O-propensity'], width=0.4, label='Permanent',edgecolor='black',color='deepskyblue')
plt.bar(x_axis+0.2,propensity_df['N-propensity'], width=0.4,label='Transient',edgecolor='black',color='limegreen')
plt.hlines(1, -1, 20,color='black')
plt.xticks(x_axis,propensity_df.iloc[:,0])
plt.xticks(rotation= 45,fontweight='bold')
#plt.yticks(np.arange(0.8,1.4,0.05),np.arange(0.8,1.4,0.05))
plt.ylabel('Normalized propensity',fontweight='bold',size=14)
plt.legend()
plt.show()
#plt.savefig('propensity.png',format='png',dpi=1300)
"""

#````````````````````````````````````````````````
#revision work for small interfaces
data3=pd.read_csv('small_intf_o_pdbs',header=None)
data4=pd.read_csv('small_intf_n_pdbs',header=None)
data3_list=list(data3.iloc[:,0])
data4_list=list(data4.iloc[:,0])
small_o_intf=[]#interfacial residues of proteins having small interface
small_n_intf=[]
for i in o_intf_dict.keys():
    if i in data3_list:
        small_o_intf.append(o_intf_dict[i])

for j in n_intf_dict.keys():
    if j in data4_list:
        small_n_intf.append(n_intf_dict[j])
small_o_intf=sum(small_o_intf,[])#making list of lists to list. All interfacial residues
small_n_intf=sum(small_n_intf,[])
small_o_res_freq=collections.Counter(small_o_intf)#frequency
small_n_res_freq=collections.Counter(small_n_intf)

small_o_propensity=[]
for i in residues:
    num21=small_o_res_freq[i]
    num22=small_o_res_freq[i]+small_n_res_freq[i]
    num23=len(small_o_intf)
    num24=len(small_o_intf)+len(small_n_intf)
    propensity=(num21/num22)/(num23/num24)
    small_o_propensity.append("%.2f" %propensity)
    
#for non-obligatory
small_n_propensity=[]
for i in residues:
    num21=small_n_res_freq[i]
    num22=small_o_res_freq[i]+small_n_res_freq[i]
    num23=len(small_n_intf)
    num24=len(small_o_intf)+len(small_n_intf)
    propensity=(num21/num22)/(num23/num24)
    small_n_propensity.append("%.2f" %propensity)
    
small_propensity_df=pd.DataFrame({'Residues':residues,'O-propensity':small_o_propensity,'N-propensity':small_n_propensity})
small_propensity_df = small_propensity_df.astype({"O-propensity":float, "N-propensity": float})
print(small_propensity_df)

"""
#plotting!!!!!!!!!!!!!!!!!!!!!!!!!!!!11
small_propensity_df['O-propensity']=pd.to_numeric(small_propensity_df['O-propensity'])
small_propensity_df['N-propensity']=pd.to_numeric(small_propensity_df['N-propensity'])
sns.set_style("darkgrid")
plt.figure(figsize=(14.0, 8.0))
x_axis=np.arange(0,20)
#obl_std=[np.std(df_obl.iloc[:,i]) for i in range(1,10)]
#non_obl_std=[np.std(df_non_obl.iloc[:,i]) for i in range(1,10)]
plt.bar(x_axis-0.2,small_propensity_df['O-propensity'], width=0.4, label='Permanent',edgecolor='black',color='deepskyblue')
plt.bar(x_axis+0.2,small_propensity_df['N-propensity'], width=0.4,label='Transient',edgecolor='black',color='limegreen')
plt.hlines(1, -1, 20,color='black')
plt.xticks(x_axis,small_propensity_df.iloc[:,0])
plt.xticks(rotation= 45,fontweight='bold')
#plt.yticks(np.arange(0.8,1.4,0.05),np.arange(0.8,1.4,0.05))
plt.ylabel('Normalized propensity',fontweight='bold',size=14)
plt.legend()
plt.show()
plt.savefig('small_propensity.png',format='png',dpi=1300)
"""

#~~~~~~~~~~~~~~CHECKING CYS IN SMALL INTERFACES
subdict={x:n_intf_dict[x] for x in data4_list}#small interface sub-dictionary
temp=[]
for i in subdict.keys():
    for j in subdict[i]:
        if j=='CYS':
            temp.append(i)
cys=collections.Counter(temp)

#delta-propensity
delta_df=pd.DataFrame({'Residues':residues,'O-propensity':propensity_df['O-propensity']-small_propensity_df['O-propensity'],'N-propensity':propensity_df['N-propensity']-small_propensity_df['N-propensity']})

#plotting
sns.set_style("darkgrid")
plt.figure(figsize=(14.0, 8.0))
x_axis=np.arange(0,20)
#obl_std=[np.std(df_obl.iloc[:,i]) for i in range(1,10)]
#non_obl_std=[np.std(df_non_obl.iloc[:,i]) for i in range(1,10)]
plt.bar(x_axis-0.2,delta_df['O-propensity'], width=0.4, label='Permanent',edgecolor='black',color='deepskyblue')
plt.bar(x_axis+0.2,delta_df['N-propensity'], width=0.4,label='Transient',edgecolor='black',color='limegreen')
#plt.hlines(1, -1, 20,color='black')
plt.xticks(x_axis,small_propensity_df.iloc[:,0])
plt.xticks(rotation= 45,fontweight='bold')
#plt.yticks(np.arange(0.8,1.4,0.05),np.arange(0.8,1.4,0.05))
plt.ylabel(r'$\Delta$ Normalized propensity',fontweight='bold',size=14)
plt.legend()
plt.show()
#plt.savefig('delta_propensity.png',format='png',dpi=1300)