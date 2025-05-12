Dir = '/nesi/nobackup/uc04105/new_databases_May/GTDB_226'
Archaea = 'Archaea_GTDB226_protein_May92025.tsv'

import json
Data_dic = {}

with open('{}/{}'.format(Dir, Archaea), 'r') as A:
    next(A, None)
    for line in A:
        Alist = []
        seq = line.split('\t')[-1].split('\n')[0]
        tax = line.split('\t')[2]
        protein_id = line.split('\t')[0]
        Alist.append(seq)
        Alist.append(tax)
        Data_dic[protein_id] = Alist

outfile = open('{}/Archaea_GTDB_protein_May92025.json'.format(Dir), 'a')
json.dump(Data_dic, outfile, indent=6)
