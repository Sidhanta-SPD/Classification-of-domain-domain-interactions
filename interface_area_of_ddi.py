#date: 10/05/21
#this program is a modificatiton of interface residue script of pymol
#from https://pymolwiki.org/index.php/InterfaceResidues#More_Complex_Example
#also available in ~/Sidhanta/programs/
#AIM: FIND INTERFACE AREA OF COMPLEX. HERE INTERFACE AREA OF DOMAIN DOMAIN INTERACTION
#removed 1LDM on 27/4/23.
#finds domain interfaces larger than 4000 sq Ang


import os
from interface_surface_area import interfaceResidues_area
import pandas as pd
#from pymol import cmd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from scipy.stats import mannwhitneyu #mann-whitney test
import scipy.stats as stats #ks 2 sample test

parent_dir= '/home/sidhant/Sidhanta/work/multidomain/energy/PPCheck/'

data=pd.read_csv('obl_pdbs', header =None)
data1=pd.read_csv('non_obl_pdbs', header= None)


interface_area_obl=[] #interface areas of obl in A^2
obl_pdb=[]#pdbs
interface_area_non_obl=[] #interface areas of non_obl in A^2
n_obl_pdb=[]#pdbs
i=0
j=0

for pdb in data.iloc[:,0]:
    #os.chdir(parent_dir+'obl_pdb_files/'+pdb)
    pdb_file=parent_dir+'obl_pdb_files/'+pdb+'_edit.pdb'
    cmd.load(pdb_file,pdb)
    cmd.disable('all')
    cmd.enable(pdb)
    #cmd.get_chains()
    a=interfaceResidues_area(pdb,cA="c. A", cB="c. B", cutoff=1, selName="foundIntrf")
    interface_area_obl.append(float(a))
    obl_pdb.append(pdb)
    i+=1
    print('done',i)

for pdb in data1.iloc[:,0]:
    #os.chdir(parent_dir+'non_obl_pdb_files/'+pdb)
    pdb_file=parent_dir+'non_obl_pdb_files/'+pdb+'_edit.pdb'
    cmd.load(pdb_file,pdb)
    cmd.disable('all')
    cmd.enable(pdb)
    #cmd.get_chains()
    a=interfaceResidues_area(pdb,cA="c. A", cB="c. B", cutoff=1, selName="foundIntrf")
    interface_area_non_obl.append(float(a))
    n_obl_pdb.append(pdb)#appending pdb ids for dict creation
    j+=1
    print('done1',j)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#plotting
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
obl_area_dict=dict(zip(obl_pdb,interface_area_obl_sasa))#dict of pdb and their interface
n_obl_area_dict=dict(zip(n_obl_pdb,interface_area_non_obl_sasa))

large_intf_obl=dict((k,v) for k, v in obl_area_dict.items() if v >=4000) #proteins having domain interface larger than or equal to 4000 sq A.
large_intf_n_obl=dict((k,v) for k, v in n_obl_area_dict.items() if v >=4000)
small_intf_obl=dict((k,v) for k,v in obl_area_dict.items() if v <= 1500)#proteins having domain interfaces smaller than or equal to 1500 sq A.
small_intf_n_obl=dict((k,v) for k,v in n_obl_area_dict.items() if v <= 1500)


#df=pd.DataFrame({'Obl':pd.Series(interface_area_obl),'Non_obl':pd.Series(interface_area_non_obl)})
data10=list([interface_area_obl_sasa,interface_area_non_obl_sasa])
#plt.figure(figsize=(13,10),dpi=80)
#sns.set_style("darkgrid")
fig,ax=plt.subplots()
#box= plt.boxplot(data10, showmeans=True,patch_artist=True)

#sns.violinplot(data=df, scale='width', inner='quartile')
#colors=['bisque', 'lightsteelblue']
#sns.boxplot(data=data10, notch=False, showmeans=True,palette=colors,meanprops={"markerfacecolor":"black"})
sns.boxplot(data=data10, notch=False, showmeans=True,meanprops={"markerfacecolor":"black"},palette=['deepskyblue','limegreen'])
sns.swarmplot(data=data10,color='blue',alpha=0.4,size=3)
plt.yticks(np.arange(0,max(interface_area_obl_sasa)+2, step=1000))
#ax.violinplot(data10, showmeans=True, showmedians=False)
ax.set_xticks([0,1])
xticklabels=['Permanent', 'Transient']
ax.set_xticklabels(xticklabels,fontweight='bold')
ax.set_xlabel('DDI types', fontweight='bold')
ax.set_ylabel(r'interface area in $\AA^2$', fontweight='bold')
#ax.set_title('Domain-domain interface areas', fontweight='bold')

"""
colors=['bisque', 'lightsteelblue']
for patch, color in zip(box['boxes'], colors):
    patch.set_facecolor(color)
"""
plt.show()
#plt.savefig('intf_area_sasa_swarm_manuscipt_no_grid.png', format='png',dpi=1300)

result=mannwhitneyu(interface_area_obl_sasa, interface_area_non_obl_sasa)
ks=stats.ks_2samp(interface_area_obl_sasa, interface_area_non_obl_sasa)

#with open('obl_large_interface','w') as filehandle:
#    filehandle.writelines('%s\t%s\n' %(key,value) for (key,value) in large_intf_obl.items())

#with open('n_obl_large_interface','w') as filehandle:
#    filehandle.writelines('%s\t%s\n' %(key,value) for (key,value) in large_intf_n_obl.items()) 

#with open('obl_small_interface', 'w') as filehandle:
#    filehandle.writelines('%s\t%s\n' %(key,value) for (key, value) in small_intf_obl.items())
#with open('n_obl_small_interface', 'w') as filehandle:
#    filehandle.writelines('%s\t%s\n' %(key,value) for (key, value) in small_intf_n_obl.items())
