#!/bin/bash

while read protein;
do
	#echo "$protein"
	IFS='.'
	read -a pdb_id <<< $protein
	prot_pdb_id=${pdb_id[0]} #pdb_code
	prot_pdb_id_lower=${prot_pdb_id,,} # lowered pdb code. useful for calling split pdb files
	dir=${prot_pdb_id,,} #lowered pdb code. this is the directory names
	#echo "$dir"

	cd ${dir}/NOXclass/python_scripts/
	python predict.py -p 0 -b ${prot_pdb_id_lower}_edit_a.pdb -c ${prot_pdb_id_lower}_edit_b.pdb -e 1_1_1_0_0_0 -m 1 -f 0 > nox_out_${dir}
	cd ../../..


done < ~/Sidhanta/work/multidomain/noxclass_trial/noxclass_inputs/nox_input_pdbs/two_domain_proteins_updated

