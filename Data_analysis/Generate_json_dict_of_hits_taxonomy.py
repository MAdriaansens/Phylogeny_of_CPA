#due to the structure of all databases this code is mostly universal for Archaea, Bacteria and Eukarya
#HOWEVER one edit for Eukarya is needed as the downstream application between Prokarya and Eukarya is a bit different

import sys


all_sequences = sys.argv[1] # '/nesi/nobackup/uc04105/results/hmmalign/rerun_pipeline_euk/complete_euk.fasta'
#protein tsv including taxonomic and species data
tsv = sys.argv[2] #'/nesi/nobackup/uc04105/database/Euk_db_7April/Euk_db_7April_protein.tsv'

from Bio import SeqIO
import time

id_list = []
for record in SeqIO.parse(all_sequences, 'fasta'):
    id_list.append(record.id)

id_list.sort()
id_list = list(set(id_list))
length = len(id_list)

dict_hits = {}
with open(tsv, 'r') as B:
    for i in B:
        #id number protein is i.spli('\t')[0] in all cases
        name = i.split('\t')[0]
        if name in id_list:
            description = []
            
            description.append(i.split('\t')[-2]) #species id
            
            description.append(i.split('\t')[-3]) #gtdb taxonomy

            #in the case of eukarya the first description should be kept the same and now returns the taxonomy
            #only edit for eukarya is to replace [-3] with [1] to get the species name
            dict_hits[name] = description
B.close()

import json
outfile = sys.argv[3]

with open("{}.json".format(outfile), "w") as out:
    json.dump(dict_hits, out)
