import os
from Bio import Phylo
from io import StringIO

from Bio import SeqIO

alignment = '/nesi/nobackup/uc04105/database/GTDB_220/ar53_ssu_reps_r220_autoaligned.fna'
tree = '/nesi/nobackup/uc04105/database/GTDB_220/ar53_ssu_reps_r220_FT_GTRCAT.treefile'
a = []

for record in SeqIO.parse(alignment, 'fasta'):
    a.append(str(record.id))

res = [(a[i], a[j]) for i in range(len(a)) for j in range(i + 1, len(a))]


TreeLoaded = Phylo.read(tree, 'newick')


pair_phylo_distance = {}
import json

for pair in res:
    i = pair[0]
    j = pair[1]
    rename = i + '*' + j
    distance = str(TreeLoaded.distance(i, j))
    pair_phylo_distance[rename] = distance

with open('/nesi/nobackup/uc04105/redvals/Phylo_distance_ssu_Archaea.json', 'w') as fp:
    json.dump(pair_phylo_distance, fp)
