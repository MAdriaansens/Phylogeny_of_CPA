This is the output directory 
In this directroy are the scripts which generate the output files and some of the output files used/generated earlier on.
It contains .tsv file of all the antiporters present as well as all the types and the .tsv files

**Supplementary**
here we have the scripts in which we generate some of our supplementary figures and do statistical comparison of our antiporters presence/absence.


**write_fl_CPA_seq**

In write_fl_CPA_seq is a script which writes a .faa and json file containing data and seq data of all CPAs we identified.
It also generated a CPA_fl_taxa.json file. This file is used later on.

Input = fasta file of all CPAs as awell as metadata file
Output = json file containing all seq and taxa info of each protein id.

**Types_CPA**

In types_CPA is a script which types of CPA are present in which organisms.
Input: CPA_fl.json, red-based clade .tsv file, as well as mmseqs.
Output: .tsv file of where each species has the count of each CPA type present in its representative proteome. 

First this file just goes trough the red-interval and makes a dictionary where each protein rep in the trees is given its clade_id.
Second it matches the representatives to the clusteroids and other proteins they represent.
Third we load the .json dictionary.

We then match each entry in the euk_dic to its representative and the clade that representative belongs too. 
We only match or include the 11 main clades of CPA. 
