import os
import pandas as pd
import re
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import scipy.stats as stats #statistical tests

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#this is for getting number of interactions for the interacting one
def count_int(line):
    j=line+3
    while bool(re.search(r'\t',int_file.iloc[j,0])):#looping interaction lines using presence of \t
        j+=1
        if j==len(int_file):break
    df=int_file.iloc[line+3:j,0].str.split('\t',expand=True)
    if df.shape[1]<10:
        df=df.drop_duplicates(subset=[0,3])#changed 4 to 3
        df[0]=pd.to_numeric(df[0])
        df[3]=pd.to_numeric(df[3])
        int_dist=abs(df[0]-df[3])
        count=(int_dist>6).sum()#at least 7 residues apart 
    elif df.shape[1]>10:
        #df=df.drop_duplicates(subset=[0,4])
        df=df[df[8]!='2']#removing rows having MO=2
        df[0]=pd.to_numeric(df[0])
        df[4]=pd.to_numeric(df[4])
        int_dist=abs(df[0]-df[4])
        count=(int_dist>6).sum()
    return(count)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


obl_file=pd.read_csv('obl_pdbs',header=None)#list of obl_pdbs wc 109. removed 1LDM.pdb. refer nb 8/2/23
non_obl_file=pd.read_csv('non_obl_pdbs',header=None)#list of non-obl pdbs. wc 154

#||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#OBLIGATORY

pdb_id=[]#ids
hydrophobic=[]#hydrophobic
ionic=[]#ionic
aro_aro=[]#aromatic-aromatic
aro_sul=[]#aromatic-sulphur
cat_pi=[]#cation-pi
disul=[]#disulphide
mc_mc=[]#main chain-main chain
mc_sc=[]#main chain-side chain
sc_sc=[]#side chain-side chain

for o_pdb in obl_file.iloc[:,0]:
    o_pdb=o_pdb.lower()
#o_pdb='1rhs'
    #print(o_pdb)
    pdb_id.append(o_pdb)
    #print('-----')
    os.chdir('/home/sidhant/Sidhanta/work/multidomain/interaction/interaction_files/'+o_pdb)
    int_file=pd.read_csv('inter_'+o_pdb+'_edit',header=None)#interaction file
    int_file.columns=['all']

    headlines=[]#index of all the headlines of interactions (E.g. Protein-protein hydrophobic....)
    for line_no in range(len(int_file)):
        if bool(re.search(r'\s--+', int_file.iloc[line_no,0])): headlines.append(line_no-1)#for the interactions
        elif bool(re.search(r'PROTEIN-PROTEIN',int_file.iloc[line_no,0])): headlines.append(line_no)#for no interactions

    #print(headlines)

    for i in range(len(headlines)):
        if i==0:#hydrophobic
           # print('hydrophobic')
            if re.search(r'PROTEIN-PROTEIN',int_file.iloc[headlines[i],0]):
                hydrophobic.append(0)
            elif re.search(r'\s--+', int_file.iloc[headlines[i]+1,0]):
                hydrophobic.append(count_int(headlines[i]))
            
        if i==1:#ionic
            #print('ionic')
            if re.search(r'PROTEIN-PROTEIN',int_file.iloc[headlines[i],0]):
                ionic.append(0)
            elif re.search(r'\s--+', int_file.iloc[headlines[i]+1,0]):
                ionic.append(count_int(headlines[i]))
            
        if i==2:#aromatic-aromatic
            #print('aromatic-aromatic')
            if re.search(r'PROTEIN-PROTEIN',int_file.iloc[headlines[i],0]):
                aro_aro.append(0)
            elif re.search(r'\s--+', int_file.iloc[headlines[i]+1,0]):
                aro_aro.append(count_int(headlines[i]))
            
        if i==3:#aromatic-sulphur
            #print('aromatic-sulphur')
            if re.search(r'PROTEIN-PROTEIN',int_file.iloc[headlines[i],0]):
                aro_sul.append(0)
            elif re.search(r'\s--+', int_file.iloc[headlines[i]+1,0]):
                aro_sul.append(count_int(headlines[i]))
            
        if i==4:#cation-pi
            #print('cation-pi')
            if re.search(r'PROTEIN-PROTEIN',int_file.iloc[headlines[i],0]):
                cat_pi.append(0)
            elif re.search(r'\s--+', int_file.iloc[headlines[i]+1,0]):
                cat_pi.append(count_int(headlines[i]))
            
        if i==5:#disulphide
            #print('disulphide')
            if re.search(r'\s--+', int_file.iloc[headlines[i]+1,0]):
                disul.append(count_int(headlines[i]))
            elif re.search(r'PROTEIN-PROTEIN',int_file.iloc[headlines[i],0]):
                disul.append(0)
        if i==6:#main chain-main chain
            #print('main chain-main chain')
            if re.search(r'PROTEIN-PROTEIN',int_file.iloc[headlines[i],0]):
                mc_mc.append(0)
            elif re.search(r'\s--+', int_file.iloc[headlines[i]+1,0]):
                mc_mc.append(count_int(headlines[i]))
            
        if i==7:#main chain-side chain
            #print('main chain-side chain')
            if re.search(r'PROTEIN-PROTEIN',int_file.iloc[headlines[i],0]):
                mc_sc.append(0)
            elif re.search(r'\s--+', int_file.iloc[headlines[i]+1,0]):
                mc_sc.append(count_int(headlines[i]))
            
        if i==8:#side chain-side chain
            #print('side chain-side chain')
            if re.search(r'PROTEIN-PROTEIN',int_file.iloc[headlines[i],0]):
                sc_sc.append(0)
            elif re.search(r'\s--+', int_file.iloc[headlines[i]+1,0]):
                sc_sc.append(count_int(headlines[i]))
            
    #os.chdir('..')

df_obl=pd.DataFrame({'pdb':pdb_id,'hphbic':hydrophobic,'ionic':ionic,'aro-aro':aro_aro,'aro-s':aro_sul,'cat-pi':cat_pi,'disulf':disul,'mc-mc':mc_mc,'mc-sc':mc_sc,'sc-sc':sc_sc})


#||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#NON-OBLIGATORY

n_pdb_id=[]#ids
n_hydrophobic=[]#hydrophobic
n_ionic=[]#ionic
n_aro_aro=[]#aromatic-aromatic
n_aro_sul=[]#aromatic-sulphur
n_cat_pi=[]#cation-pi
n_disul=[]#disulphide
n_mc_mc=[]#main chain-main chain
n_mc_sc=[]#main chain-side chain
n_sc_sc=[]#side chain-side chain

for o_pdb in non_obl_file.iloc[:,0]:
    o_pdb=o_pdb.lower()
#o_pdb='1rhs'
    #print(o_pdb)
    n_pdb_id.append(o_pdb)
    #print('-----')
    os.chdir('/home/sidhant/Sidhanta/work/multidomain/interaction/interaction_files/'+o_pdb)
    int_file=pd.read_csv('inter_'+o_pdb+'_edit',header=None)#interaction file
    int_file.columns=['all']

    headlines=[]#index of all the headlines of interactions (E.g. Protein-protein hydrophobic....)
    for line_no in range(len(int_file)):
        if bool(re.search(r'\s--+', int_file.iloc[line_no,0])): headlines.append(line_no-1)#for the interactions. As it contains heading followed by -------
        elif bool(re.search(r'PROTEIN-PROTEIN',int_file.iloc[line_no,0])): headlines.append(line_no)#for no interactions

    #print(headlines)

    for i in range(len(headlines)):
        if i==0:#hydrophobic
            #print('hydrophobic')
            if re.search(r'PROTEIN-PROTEIN',int_file.iloc[headlines[i],0]):
                n_hydrophobic.append(0)
            elif re.search(r'\s--+', int_file.iloc[headlines[i]+1,0]):
                n_hydrophobic.append(count_int(headlines[i]))
            
        if i==1:#ionic
            #print('ionic')
            if re.search(r'PROTEIN-PROTEIN',int_file.iloc[headlines[i],0]):
                n_ionic.append(0)
            elif re.search(r'\s--+', int_file.iloc[headlines[i]+1,0]):
                n_ionic.append(count_int(headlines[i]))
            
        if i==2:#aromatic-aromatic
            #print('aromatic-aromatic')
            if re.search(r'PROTEIN-PROTEIN',int_file.iloc[headlines[i],0]):
                n_aro_aro.append(0)
            elif re.search(r'\s--+', int_file.iloc[headlines[i]+1,0]):
                n_aro_aro.append(count_int(headlines[i]))
            
        if i==3:#aromatic-sulphur
            #print('aromatic-sulphur')
            if re.search(r'PROTEIN-PROTEIN',int_file.iloc[headlines[i],0]):
                n_aro_sul.append(0)
            elif re.search(r'\s--+', int_file.iloc[headlines[i]+1,0]):
                n_aro_sul.append(count_int(headlines[i]))
            
        if i==4:#cation-pi
            #print('cation-pi')
            if re.search(r'PROTEIN-PROTEIN',int_file.iloc[headlines[i],0]):
                n_cat_pi.append(0)
            elif re.search(r'\s--+', int_file.iloc[headlines[i]+1,0]):
                n_cat_pi.append(count_int(headlines[i]))
            
        if i==5:#disulphide
            #print('disulphide')
            if re.search(r'\s--+', int_file.iloc[headlines[i]+1,0]):
                n_disul.append(count_int(headlines[i]))
            elif re.search(r'PROTEIN-PROTEIN',int_file.iloc[headlines[i],0]):
                n_disul.append(0)
        if i==6:#main chain-main chain
            #print('main chain-main chain')
            if re.search(r'PROTEIN-PROTEIN',int_file.iloc[headlines[i],0]):
                n_mc_mc.append(0)
            elif re.search(r'\s--+', int_file.iloc[headlines[i]+1,0]):
                n_mc_mc.append(count_int(headlines[i]))
            
        if i==7:#main chain-side chain
            #print('main chain-side chain')
            if re.search(r'PROTEIN-PROTEIN',int_file.iloc[headlines[i],0]):
                n_mc_sc.append(0)
            elif re.search(r'\s--+', int_file.iloc[headlines[i]+1,0]):
                n_mc_sc.append(count_int(headlines[i]))
            
        if i==8:#side chain-side chain
            #print('side chain-side chain')
            if re.search(r'PROTEIN-PROTEIN',int_file.iloc[headlines[i],0]):
                n_sc_sc.append(0)
            elif re.search(r'\s--+', int_file.iloc[headlines[i]+1,0]):
                n_sc_sc.append(count_int(headlines[i]))
            
    #os.chdir('..')

df_non_obl=pd.DataFrame({'pdb':n_pdb_id,'hphbic':n_hydrophobic,'ionic':n_ionic,'aro-aro':n_aro_aro,'aro-s':n_aro_sul,'cat-pi':n_cat_pi,'disulf':n_disul,'mc-mc':n_mc_mc,'mc-sc':n_mc_sc,'sc-sc':n_sc_sc})



df_obl['total_int']=df_obl['hphbic']+df_obl['ionic']+df_obl['aro-aro']+df_obl['aro-s']+df_obl['cat-pi']+df_obl['disulf']+df_obl['mc-mc']+df_obl['mc-sc']+df_obl['sc-sc']
df_non_obl['total_int']=df_non_obl['hphbic']+df_non_obl['ionic']+df_non_obl['aro-aro']+df_non_obl['aro-s']+df_non_obl['cat-pi']+df_non_obl['disulf']+df_non_obl['mc-mc']+df_non_obl['mc-sc']+df_non_obl['sc-sc']
df_obl['total_int']=pd.to_numeric(df_obl['total_int'])
df_non_obl['total_int']=pd.to_numeric(df_non_obl['total_int'])

#print(df_obl)
#print(df_non_obl)

"""
#verification mode
#remove comments when running
df_obl.to_csv('/home/sidhant/Sidhanta/work/multidomain/interaction2.0/obligatory_interactions', index=False)
df_non_obl.to_csv('/home/sidhant/Sidhanta/work/multidomain/interaction2.0/non_obligatory_interactions',index=False)
"""

os.chdir('/home/sidhant/Sidhanta/work/multidomain/interaction2.0/normalized_random')

#df_obl['hphbic'].describe()

#^^^^^^^^ NORMALITY CHECK ^^^^^^^^^^
#refer nb 15/6/22
#mostly ks test is thought to have lower statistical power
for i in range(1,10):
    if stats.shapiro(df_obl.iloc[:,i])[1]> 0.05: print (df_obl.columns[i],' is normally distributed')
print('---')
for i in range(1,10):
    if stats.shapiro(df_non_obl.iloc[:,i])[1]> 0.05: print (df_non_obl.columns[i],' is normally distributed')
#it doesn't print anything. Hence these are non-normal data




# @@@@@@@@@@@@@@@@@ processing @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

#*************************************
df_non_obl_rand=df_non_obl.sample(n=109)#take random rows of non-obl. normalizing the no. of proteins
#*************************************
#normalization-----------------------------------------------------------------------
#normalizing to total numbr of interaction in the protein
df_obl.iloc[:,1:10]=df_obl.iloc[:,1:10].div(df_obl.total_int,axis=0)
df_non_obl_rand.iloc[:,1:10]=df_non_obl_rand.iloc[:,1:10].div(df_non_obl_rand.total_int,axis=0)
df_non_obl.iloc[:,1:10]=df_non_obl.iloc[:,1:10].div(df_non_obl.total_int,axis=0)#normalizing 



#as Himani and Tarun suggested: merge sidechain h-bonds and remove disulphide
df_obl['sch']=df_obl['mc-sc']+df_obl['sc-sc'] #mc-sc and sc-sc together
df_obl=df_obl[['pdb', 'hphbic', 'ionic', 'aro-aro', 'aro-s', 'cat-pi', 'mc-mc', 'sch', 'total_int']]
df_non_obl_rand['sch']=df_non_obl_rand['mc-sc']+df_non_obl_rand['sc-sc']
df_non_obl_rand=df_non_obl_rand[['pdb', 'hphbic', 'ionic', 'aro-aro', 'aro-s', 'cat-pi', 'mc-mc', 'sch', 'total_int']]
df_non_obl['sch']=df_non_obl['mc-sc']+df_non_obl['sc-sc']
df_non_obl=df_non_obl[['pdb', 'hphbic', 'ionic', 'aro-aro', 'aro-s', 'cat-pi', 'mc-mc', 'sch', 'total_int']]



print(df_obl)
print(df_non_obl_rand)


#^^^^^^^^ statistical significance test ^^^^^^^^^^^
for i in range(1,8):
    for j in range(1,8):#7 interaction types
        if i == j:
            p_value= stats.mannwhitneyu(df_obl.iloc[:,i],df_non_obl_rand.iloc[:,j])[1]
            print(df_obl.columns[i],':',p_value)#mann-whitney test for each interaction type for obl and non-obl.
            if p_value < 0.005: print(df_obl.columns[i],' is statistically significant')

"""
#plotting @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@2
#@@@@@@@@@@@@@@@ pie plot @@@@@@@@@@@@@
#+++++use df_non_obl for normalized original values (154) and df_non_obl_rand for normalized random (109)
obl_avg=[(df_obl[df_obl.columns[i]].sum()/len(df_obl)).round(3) for i in range(1,len(df_obl.columns)-1)]
#non_obl_avg=[(df_non_obl_rand[df_non_obl_rand.columns[i]].sum()/len(df_non_obl_rand)).round(3) for i in range(1,len(df_non_obl_rand.columns)-1)]
non_obl_avg=[(df_non_obl[df_non_obl.columns[i]].sum()/len(df_non_obl)).round(3) for i in range(1,len(df_non_obl.columns)-1)]

labels=['Hydrophobic','Ionic','Aromatic-aromatic','Aromatic-sulphur','Cation-pi','Main chain-main chain \n H-bonds','Side chain \n H-bonds']

#its not a great idea to plot average values of each interactions for pie plot
obl_sum=[(df_obl[df_obl.columns[i]].sum()).round(3) for i in range(1,len(df_obl.columns)-1)]
#non_obl_sum=[(df_non_obl_rand[df_non_obl_rand.columns[i]].sum()).round(3) for i in range(1,len(df_non_obl_rand.columns)-1)]
non_obl_sum=[(df_non_obl[df_non_obl.columns[i]].sum()).round(3) for i in range(1,len(df_non_obl.columns)-1)]

#explode = [0.1,0,0,0,0,0,0,0,0]
plt.figure(figsize=(14.0, 8.0))
plt.pie(obl_sum,labels=labels,autopct='%.0f%%',textprops={'fontsize':16})
plt.show()
#plt.savefig("pie_obl_final.png", format="png", dpi=1300)
plt.pie(non_obl_sum, labels=labels, autopct='%.0f%%',textprops={'fontsize':16})
plt.show()
#plt.savefig("pie_n_obl_final.png",format="png",dpi=1300)

#@@@@@@@@@@@@@@@@@@@@@ pairplot @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
df_obl['type']=['Obligatory' for i in range(len(df_obl))]
df_non_obl['type']=['Non-obligatory' for i in range(len(df_non_obl))]
#as the number of non-obl entries are more, the frequrency/peak of non-obl plots are siturated higher. So trying to take random 110 non-obl entries. 
df_non_obl_rand=df_non_obl.sample(n=109)#take random rows
#df_both=df_obl.append(df_non_obl)#for taking all non-obl entries
df_both=df_obl.append(df_non_obl_rand)#taking random 110 non-obl entries
sns.pairplot(df_both, hue='type')
plt.show()

#for h-bonds only~~~~~~~~~~~~~~~~~~~~
df_obl_h=df_obl.loc[:,('mc-mc','mc-sc','sc-sc')]
df_obl_h['type']=['Obligatory' for i in range(len(df_obl_h))]
df_nobl_h=df_non_obl_rand.loc[:,('mc-mc','mc-sc','sc-sc')]
df_nobl_h['type']=['Non-obligatory' for i in range(len(df_nobl_h))]
df_both_h=df_obl_h.append(df_nobl_h)
sns.pairplot(df_both_h,hue='type')
plt.show()
sns.scatterplot(x=df_both_h['mc-sc'],y=df_both_h['sc-sc'],hue=df_both_h.type)
plt.show()
plt.savefig('mc-sc_sc-sc.png',format='png',dpi=1300)

#barplot of types of interactions @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
x_axis=np.arange(0,7)
plt.figure(figsize=(8,6))
#obl_std=[np.std(df_obl.iloc[:,i]) for i in range(1,10)]
#non_obl_std=[np.std(df_non_obl.iloc[:,i]) for i in range(1,10)]
plt.bar(x_axis-0.2,obl_avg, width=0.4, label='Obligatory',edgecolor='black')
plt.bar(x_axis+0.2,non_obl_avg, width=0.4,label='Non-obligatory',edgecolor='black')
plt.xticks(x_axis,labels)
plt.xticks(rotation= 45,ha='right')
plt.ylabel('Average no. of interactions')
plt.legend()
plt.tight_layout()
plt.show()
#plt.savefig("bar_final.png",format="png",dpi=1300)

#histplot of total interactions @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
bins=[0,15,30,45,60,75,90,105,120,135,150]
plt.hist([df_obl.total_int,df_non_obl_rand.total_int], bins=bins,edgecolor='black')
plt.xticks([x*15 for x in range(11)],[x*15 for x in range(11)])
plt.ylabel('Count of proteins',fontweight='bold')
plt.xlabel('Range of number of interactions',fontweight='bold')
plt.legend(['Obligatory','Non-obligatory'])
plt.show()
#plt.savefig('tot_int_rand3.png',format='png',dpi=1300)
"""
