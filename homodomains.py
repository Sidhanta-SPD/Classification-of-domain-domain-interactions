import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

obl_cla_file=pd.read_csv('obl_cla_file',sep='\t',header=None)
non_obl_cla_file=pd.read_csv('non_obl_cla_file',sep='\t',header=None)
obl_file=pd.read_csv('obl_pdbs',header=None)
non_obl_file=pd.read_csv('non_obl_pdbs',header=None)


#````````````````````````````````````````````````````````````````
o_sccs=[]#all sccs ids.wc 220
for i in range(len(obl_file)):
    #print(i)
    for j in range(len(obl_cla_file)):     
        if obl_file.iloc[i,0].lower()== obl_cla_file.iloc[j,1]:
            #print(j)
            o_sccs.append(obl_cla_file.iloc[j,3])
            #print(obl_cla_file.iloc[j,3])
n_sccs=[]#wc308
for i in range(len(non_obl_file)):
    for j in range(len(non_obl_cla_file)):     
        if non_obl_file.iloc[i,0].lower()== non_obl_cla_file.iloc[j,1]:
            #print(j)
            n_sccs.append(non_obl_cla_file.iloc[j,3])
#`````````````````````````````````````````````````````````````````    

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
def uniq_homodomain(listt):#function to return a list of uniq sccs. uniq till superfamily level.
    for i in range(len(listt)):
        for j in range(len(listt)):
            if i!=j:
                if listt[i].split('.')[0:3]==listt[j].split('.')[0:3]:
                    listt[j]=listt[i]#if sccs is same till superfamily, updating the j to i.
    x=list(sorted(set(listt),key=listt.index))#sorting and preserving the order according to original index
    return(x)
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

#o_uniq_supfam=uniq_homodomain(o_sccs)#uniq sccs till superfamily
n_uniq_supfam=uniq_homodomain(n_sccs)#uniq_sccs till superfamily in non-obl.wc 149
uniq_supfam=[i for i in uniq_homodomain(o_sccs)]#uniq sccs till superfamily in both. At first obl uniq sccs were appended (wc 110) to have a order.
for j in n_uniq_supfam:#appedning uniq non-obl sccs which are not there in uniq_supfam or obl uniq sccs.
    if j in uniq_supfam:
        continue
    else:
        uniq_supfam.append(j)
uniq_supfam=uniq_homodomain(uniq_supfam)#229 uniq sccs till supfam in both obl and non-obl
#uniq_sccs=list(set(all_sccs))#uniq sccs ids in whole dataset
#uniq_supfam=uniq_homodomain(uniq_sccs)#uniq till superfamily level in all dataset.maps are according to this list index


#trying to order obl and non-obl supfamily maps in ascending order.
ordered_uniq_supfam=[]
for i in range(len(obl_file)):
    for j in range(len(obl_cla_file)):
        for k in range(len(uniq_supfam)):
            if obl_file.iloc[i,0].lower()==obl_cla_file.iloc[j,1] and obl_cla_file.iloc[j,3].split('.')[0:3]==uniq_supfam[k].split('.')[0:3]:
                ordered_uniq_supfam.append(uniq_supfam[k])

for i in range(len(non_obl_file)):
    for j in range(len(non_obl_cla_file)):
        for k in range(len(uniq_supfam)):
            if non_obl_file.iloc[i,0].lower()==non_obl_cla_file.iloc[j,1] and non_obl_cla_file.iloc[j,3].split('.')[0:3]==uniq_supfam[k].split('.')[0:3]:
                if uniq_supfam[k] in ordered_uniq_supfam:
                    continue
                else:
                    ordered_uniq_supfam.append(uniq_supfam[k])
        

o_pdb_id=[]
o_sccs1=[]
o_sccs2=[]
o_map1=[]
o_map2=[]
for i in range(len(obl_file)):
    temp=[]
    for j in range(len(obl_cla_file)):
        if obl_file.iloc[i,0].lower()==obl_cla_file.iloc[j,1]:
            for k in range(len(uniq_supfam)):
               if obl_cla_file.iloc[j,3].split('.')[0:3]==uniq_supfam[k].split('.')[0:3]:
                   temp.append(obl_cla_file.iloc[j,3])
                   temp.append(k)
    print(i,obl_file.iloc[i,0],temp)
    o_pdb_id.append(obl_file.iloc[i,0])
    o_sccs1.append(temp[0])
    o_sccs2.append(temp[2])
    o_map1.append(temp[1])
    o_map2.append(temp[3])

#```````````````````````````````````````````````````````````````````````````````````````````````

n_pdb_id=[]
n_sccs1=[]
n_sccs2=[]
n_map1=[]
n_map2=[]
for i in range(len(non_obl_file)):
    temp=[]
    for j in range(len(non_obl_cla_file)):
        if non_obl_file.iloc[i,0].lower()==non_obl_cla_file.iloc[j,1]:
            for k in range(len(uniq_supfam)):
                if non_obl_cla_file.iloc[j,3].split('.')[0:3]==uniq_supfam[k].split('.')[0:3]:
                    temp.append(non_obl_cla_file.iloc[j,3])
                    temp.append(k)
    print(i,non_obl_file.iloc[i,0],temp)
    n_pdb_id.append(non_obl_file.iloc[i,0])
    n_sccs1.append(temp[0])
    n_sccs2.append(temp[2])
    n_map1.append(temp[1])
    n_map2.append(temp[3])



obl_df=pd.DataFrame(list(zip(o_pdb_id,o_sccs1,o_sccs2,o_map1,o_map2)),columns=['PDB','SCCS1','SCCS2','MAP1','MAP2'])
obl_df['type']=['Permanent' for i in range(len(obl_df))]
n_obl_df=pd.DataFrame(list(zip(n_pdb_id,n_sccs1,n_sccs2,n_map1,n_map2)),columns=['PDB','SCCS1','SCCS2','MAP1','MAP2'])
n_obl_df['type']=['Transient' for i in range(len(n_obl_df))]
n_obl_df_rand=n_obl_df.sample(n=110)
df_both=obl_df.append(n_obl_df)
df_both['Homodomain']=df_both.apply(lambda row: row.MAP1 - row.MAP2==0, axis=1)
df_both['Homodomain'] = df_both['Homodomain'].map({True: 'yes', False: 'no'})
sns.set_style('darkgrid')
#sns.scatterplot(data=df_both,x=df_both['MAP1'],y=df_both['MAP2'],hue=df_both.type,style=df_both.Homodomain)
df_hom=df_both[df_both.Homodomain=='yes']
sns.scatterplot(data=df_hom,x=df_hom['MAP1'],y=df_hom['MAP2'],hue=df_hom.type,size=1)
plt.xlabel('Domain 1')
plt.ylabel('Domain 2')
#plt.savefig('new_homodomain.png',format='png',dpi=1300)
plt.show()
