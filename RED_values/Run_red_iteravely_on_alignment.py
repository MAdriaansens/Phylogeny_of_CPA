from redvals import RedTree

red_trees = RedTree("decorated_trees/bac120_r220_decorated.pkl", "decorated_trees/ar53_r220_decorated.pkl")
from Bio import SeqIO

import json
#with open(diction) as json_file:
       # data = json.load(json_file)
#redvalid_dict = {}
#for i in data:
 #   output = data[i]
 #   GTDB_id = output[0]
 #   taxonomy = output[1]
 #   redval_id = red_trees.get_redvals_id(GTDB_id)
 #   important = [GTDB_id, redval_id, taxonomy]
 #   redvalid_dict[i] = important




alignment = '/nesi/nobackup/uc04105/database/GTDB_220/ar53_ssu_reps_r220_autoaligned.fna'
tree = '/nesi/nobackup/uc04105/database/GTDB_220/ar53_ssu_reps_r220_FT_GTRCAT.treefile'
a = []

for record in SeqIO.parse(alignment, 'fasta'):
    a.append(str(record.id))

res = [(a[i], a[j]) for i in range(len(a)) for j in range(i + 1, len(a))]



with open('red_of_all_pairs_Archaea_ssu_r220.tsv', 'a') as R:
    header = 'pair' + '\t' + 'RED' + '\n'
    R.write(header)
    for pair in res:
        x = red_trees.get_redvals_id(pair[0])
        #if redvalid_dict needs to be used: redvalid_dict[pair[0]][1]

        y = = red_trees.get_redvals_id(pair[1])
         #if redvalid_dict needs to be used: redvalid_dict[pair[1]][1]
        RED = red_trees.dist_between_nodes(x, y)[0]
        line = pair[0] + '*' + pair[1] + '\t' + str(RED) + '\n'
        R.write(line)
