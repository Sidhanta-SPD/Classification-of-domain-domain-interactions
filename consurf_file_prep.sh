#!/bin/bash

for i in $(less non_obl_pdbs_nError);
do
	echo $i
	pdb_strip=${i:0:4} #taking just the pdb id
	#echo ${pdb_strip,,}
	cd non_obl_consurf-db/${i}_ConSurf_DB_Outputs
	#ls *_consurf_summary.txt
	sed '/^\-/d;/^\t/d;/^$/d;/^[\*a-z]/d' *_consurf_summary.txt  > ${pdb_strip,,}_consurf_summary
	sed -i '2d' ${pdb_strip,,}_consurf_summary
	#head *_consurf_summary
	cd ../..
done	
