import os
import pandas as pd
pd.options.mode.chained_assignment = None
#from pymol import cmd
import numpy as np
#from interface_res_pymol import interfaceResidues
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu #mann-whitney test
import scipy.stats as stats

#obl_file=pd.read_csv('obl_pdbs_nError', header=None)#obl pdb list
#non_obl_file=pd.read_csv('non_obl_pdbs_nError',header=None)

#~~~~~~~~~~~~~~~~~~~~~interfacial residues

parent_dir= '/home/sidhant/Sidhanta/work/multidomain/energy/PPCheck/'
trial_pdb=parent_dir+'obl_pdb_files/1md8_edit.pdb'
cmd.load(trial_pdb,"pdb")
cmd.disable('all')
cmd.enable("pdb")
intf= interfaceResidues("pdb",cA="c. A", cB="c. B", cutoff=10, selName="foundIntrf")#2nd entry of tuple is the residue number
int_res=[intf[i][1] for i in range(len(intf))]#interface residue numbers





interface_res=[]#list of interface residue number
interface_res_dict={}#dict of pdb id and interface residue number
conserved_interface_obl=[]#list of len/count of conserved interfacial residues in proteins
conserved_int_dict={}#dictionary of pdb and its conserved interfacial residues
#current_dir='/home/sidhant/Sidhanta/work/multidomain/evolution/conservation/'
j=0 #for counting

for pdb in obl_file.iloc[:,0]:
    pdb_low=pdb[0:4].lower()
    #interfacial residues----------------
    os.chdir('/home/sidhant/Sidhanta/work/multidomain/energy/PPCheck/obl_pdb_files')
    #os.chdir(parent_dir)
    trial_pdb=pdb_low + '_edit.pdb'
    cmd.load(trial_pdb, pdb)
    cmd.disable('all')
    cmd.enable(pdb)
    print(pdb_low)
    intf=interfaceResidues(pdb,cA="c. A", cB="c. B", cutoff=1,selName="found") #2nd entry of tuple is the residue number
    #print(intf)
    int_res=[intf[i][1] for i in range(len(intf))]#interface residue numbers
    key=pdb_low
    interface_res_dict[key]=int_res
    #print(intf)
    #cmd.quit()
    interface_res.append(int_res)
    j+=1
    print('done obl',j)
    
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #!!!!!!!!!!!!!!!!!conserved residues from consurf
    os.chdir('/home/sidhant/Sidhanta/work/multidomain/evolution/conservation/obl_consurf-db/'+pdb+'_ConSurf_DB_Outputs')
    consurf_file=pd.read_csv(pdb_low+'_consurf_summary',header=None,sep='\t+',engine='python')
    consurf_file.drop(0,axis=0,inplace=True) #dropping the headers
    consurf_file.reset_index(level=0,drop=True,inplace=True)
    consurf_file.columns=['POS','SEQ','3LATOM','SCORE','COLOR','CONFIDENCE_INTERVAL','CONFIDENCE_INTERVAL_COLORS','MSA_DATA','RESIDUE_VARIETY']
    for i in range(len(consurf_file)):
        consurf_file.loc[i]['COLOR']=consurf_file.loc[i]['COLOR'][2]#dealing with entries having * in color
    consurf_file['COLOR']=pd.to_numeric(consurf_file['COLOR'])
    consurf_file1=consurf_file[consurf_file['COLOR']>7]#filtering if color > 7 as these residues are thought to be conserved
    consurf_file1.reset_index(level=0,drop=True,inplace=True)
    consurf_file1['3LATOM']=consurf_file1['3LATOM'].astype('string')
    consurf_file1['3LATOM']=consurf_file1['3LATOM'].str.strip() #strip whitespace
    consurf_file2=consurf_file1[consurf_file1['3LATOM']!='-']#i don't want such conserved residues which is not available in structure
    consurf_file2.reset_index(level=0,drop=True,inplace=True)
    conserved_res=[]#conserved residues of the protein;whose color value > 7;residue should be present in structure
    for i in range(len(consurf_file2)):
        conserved_res.append(consurf_file2.loc[i]['3LATOM'].split(':')[0][3:])
    
    
    con_int=[x for x in int_res if x in conserved_res]
    print(len(con_int))
    conserved_interface_obl.append(len(con_int))
    conserved_int_dict[key]=len(con_int)
    print('conserved interfaces:',con_int)
    print('~~~~~~')
    
#~~~~~~~~~~~~~~~~~~~~~~~~NON-OBLIGATORY~~~~~~~~~~~~~~~~~~~~

interface_res_n=[]
interface_res_dict_n={}
conserved_interface_obl_n=[]
conserved_int_dict_n={}#dictionary of pdb and its conserved interfacial residues
j=0
for pdb1 in non_obl_file.iloc[:,0]:
    pdb_low=pdb1[0:4].lower()
    #interfacial residues----------------
    os.chdir('/home/sidhant/Sidhanta/work/multidomain/energy/PPCheck/non_obl_pdb_files')
    #os.chdir(parent_dir)
    trial_pdb1=pdb_low + '_edit.pdb'
    cmd.load(trial_pdb1, pdb1)
    cmd.disable('all')
    cmd.enable(pdb1)
    print(pdb_low)
    intf1=interfaceResidues(pdb1,cA="c. A", cB="c. B", cutoff=1,selName="found") #2nd entry of tuple is the residue number
    #print(intf)
    int_res1=[intf1[i][1] for i in range(len(intf1))]#interface residue numbers
    key=pdb_low
    interface_res_dict_n[key]=int_res1
    #print(intf)
    #cmd.quit()
    interface_res_n.append(int_res1)
    j+=1
    print('done non-obl',j)
    
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #!!!!!!!!!!!!!!!!!conserved residues from consurf
    os.chdir('/home/sidhant/Sidhanta/work/multidomain/evolution/conservation/non_obl_consurf-db/'+pdb1+'_ConSurf_DB_Outputs')
    consurf_file=pd.read_csv(pdb_low+'_consurf_summary',header=None,sep='\t+',engine='python')
    consurf_file.drop(0,axis=0,inplace=True) #dropping the headers
    consurf_file.reset_index(level=0,drop=True,inplace=True)
    consurf_file.columns=['POS','SEQ','3LATOM','SCORE','COLOR','CONFIDENCE_INTERVAL','CONFIDENCE_INTERVAL_COLORS','MSA_DATA','RESIDUE_VARIETY']
    for i in range(len(consurf_file)):
        consurf_file.loc[i]['COLOR']=consurf_file.loc[i]['COLOR'][2]#dealing with entries having * in color
    consurf_file['COLOR']=pd.to_numeric(consurf_file['COLOR'])
    consurf_file1=consurf_file[consurf_file['COLOR']>7]#filtering if color > 7 as these residues are thought to be conserved
    consurf_file1.reset_index(level=0,drop=True,inplace=True)
    consurf_file1['3LATOM']=consurf_file1['3LATOM'].astype('string')
    consurf_file1['3LATOM']=consurf_file1['3LATOM'].str.strip() #strip whitespace
    consurf_file2=consurf_file1[consurf_file1['3LATOM']!='-']#i don't want such conserved residues which is not available in structure
    consurf_file2.reset_index(level=0,drop=True,inplace=True)
    conserved_res=[]#conserved residues of the protein;whose color value > 7;residue should be present in structure
    for i in range(len(consurf_file2)):
        conserved_res.append(consurf_file2.loc[i]['3LATOM'].split(':')[0][3:])
    
    
    con_int=[x for x in int_res1 if x in conserved_res]
    print(len(con_int))
    conserved_interface_obl_n.append(len(con_int))
    conserved_int_dict_n[key]=len(con_int)
    #print('conserved interfaces:',con_int)
    print('~~~~~~')

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


conserved_interface_obl_n_rand=[conserved_interface_obl_n[i] for i in np.random.randint(0,147,size=102)]

result=mannwhitneyu(conserved_interface_obl, conserved_interface_obl_n_rand)
print('mann-whitney p-value: ',result)
result1=stats.ks_2samp(conserved_interface_obl, conserved_interface_obl_n_rand)
print('ks-2sample p-value: ', result1)
#++++++++++++++++++++++++++PLOTTING
data10=list([conserved_interface_obl,conserved_interface_obl_n_rand])
fig,ax=plt.subplots()
#box=plt.boxplot(data10, showmeans=True,patch_artist=True)
colors=['bisque', 'lightsteelblue']
sns.boxplot(data=data10,showmeans=True,palette=colors,meanprops={"markerfacecolor":"black"})
sns.swarmplot(data=data10,color='blue',size=3,alpha=0.4)
ax.set_xticks([0,1])
xticklabels=['Obligatory','Non-obligatory']
ax.set_xticklabels(xticklabels,fontweight='bold')
ax.set_ylabel('Number of conserved interfacial residues',fontweight='bold')
plt.savefig('cons_rand3.png',format='png',dpi=1300)
plt.show()
