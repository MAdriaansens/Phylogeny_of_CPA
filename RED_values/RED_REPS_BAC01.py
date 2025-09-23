import os
import sys
import json

tsv = sys.argv[1]
output = tsv.split('_chunk_')[1].split('.tsv')[0]
Allseq_dict = {}
    #make a dictionairy for representatives and what their GTDB_id is

with open('{}'.format(tsv) ,'r') as Proteins:
   next(Proteins, None)
   for protein in Proteins:
       protein_id = protein.split('\t')[0]
       protein_tax = protein.split('\t')[2].replace(' ', '_')
       
       Allseq_dict[protein_id] = protein_tax
       
print(len(set(list(Allseq_dict.keys()))))

with open('/nesi/nobackup/uc04105/new_databases_May/final_tree_set/sequences/MMseq/All_Bac_seq_GTDB226_tax_{}.json'.format(output), 'w') as outfile:
    json.dump(Allseq_dict, outfile)
