# Comparative analysis of permanent and transient domain–domain interactions in multi-domain proteins


This repo contains files related to the above article. [Link](https://onlinelibrary.wiley.com/doi/10.1002/prot.26581)

## Dataset creation
**final\_monomers\_refined** --- list of monomeric two domain PDB ids  
**uniprots** --- uniprot ids of the above PDBs  
**fetch\_uniprot.sh** --- sh to fetch uniprot sequences  
`bash fetch_uniprot.sh`  
**database.fsa** --- all the fetched uniprot sequences are stored as a local database  
_To run cd-hit and cluster:_  
`~/softwares/cd-hit/./cd-hit -i database.fsa -o clusters -c 0.4 -n 2 -sc 1`  

The classification/prediction was performed using [NOXClass](https://doi.org/10.1186/1471-2105-7-27)  
**noxclass\_run.sh** --- sh to run the tool in a batch for several domain pairs  

## Programs related to analysis of different properties of domain interaction types
**interaction23.py** --- Py to calculate total number of interacions between domain pairs. Also filters interacting domain pairs from non-interacting pairs using 5-5 domain interaction rule. This also finds different types of interactions using in-house software [PIC](https://doi.org/10.1186/1471-2105-7-27)  
**interface\_surface\_area.py** --- Pymol script for calculating interface surface area using SASA. Contains minor modifications for desired output.  
**interface\_area\_of\_ddi.py** --- Py to calculate domain interface areas and plot  
**energy\_processing.py** --- Py to extract calculated interaction energies and process for analysis and ploting. [PPCheck](https://pmc.ncbi.nlm.nih.gov/articles/PMC4578551/) was used to calculate energies.  
**propensity.py** --- Py to calculate amino acid propensities at the domain interfaces.  
**random\_non\_obl.ipynb** --- Jupyter notebook for calculating residue cross-correlations using ANM in Prody.  
**consurf\_file\_prep.sh** --- sh to post-process consurf-db outputs for further analysis  
**conserved\_interfacial.py** --- Py to calculated proportion of conserved interfacial residues in two domain interaction types.  
**homodomains.py** --- Py to calculate number of homodomain containing two-domain proteins  
**folds.py** --- Py to find folds that are populated in domain interaction types.  
