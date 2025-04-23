#due to the structure of all databases this code is universal for Archaea, Bacteria and Eukarya

import sys


all_sequences = sys.argv[1] # '/nesi/nobackup/uc04105/results/hmmalign/rerun_pipeline_euk/complete_euk.fasta'

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
        name = i.split('\t')[0]
        if name in id_list:
            description = []
            description.append(i.split('\t')[-2])
            description.append(i.split('\t')[-3])
            dict_hits[name] = description
B.close()

import json
outfile = sys.argv[3]

with open("{}.json".format(outfile), "w") as out:
    json.dump(dict_hits, out)
